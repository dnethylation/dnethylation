#!/usr/bin/env python
"""
STEP3C_results.py — figures and Spearman correlations for deconvolution output.

INPUT:
    deconvolution_results.tsv               (from STEP3C_deconvolution.py)

OUTPUT:
    step3c_spearman_table.tsv               (supplementary)
    step3c_fig1_all_tissues.png             (Figure 12)
    step3c_fig2_top_celltypes.png           (Figure 13)
    step3c_fig3_granulocyte_scatter.png     (Figure 14)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

plt.rcParams.update({
    "font.family": "Arial", "font.size": 11,
    "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.labelsize": 12, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
})

# Colours
AGE_COLORS = {"Young": "#E69F00", "Middle": "#009E73", "Old": "#0072B2"}
CB_TEAL = "#009E73"; CB_GREY = "#999999"; CB_HIGHLIGHT = "#D55E00"

# ── LOAD ──────────────────────────────────────────────────────────────────────
print("Loading deconvolution_results.tsv...")
df = pd.read_csv("deconvolution_results.tsv", sep="\t")

def assign_group(age):
    if age <= 40: return "Young"
    elif age <= 70: return "Middle"
    else: return "Old"

df["age_group"] = df["age"].apply(assign_group)
df["age_numeric"] = df["age"].astype(float)

META = ["sample","age","age_group","age_numeric"]
cell_cols = [c for c in df.columns if c not in META]
GROUP_ORDER = ["Young","Middle","Old"]

print(f"  Samples: {len(df)}, Cell types: {len(cell_cols)}")

# ── SPEARMAN CORRELATIONS ─────────────────────────────────────────────────────
print("\n── Computing Spearman correlations ──────────────────────────")
sp_rows = []
for ct in cell_cols:
    rho, p = stats.spearmanr(df["age_numeric"], df[ct])
    sp_rows.append({"cell_type": ct, "rho": round(rho, 3), "pval": p,
                     "direction": "increases" if rho > 0 else "decreases",
                     "significant": p < 0.05})

sp_df = pd.DataFrame(sp_rows).sort_values("pval").reset_index(drop=True)
sp_df.to_csv("step3c_spearman_table.tsv", sep="\t", index=False)
n_sig = sp_df["significant"].sum()
print(f"  Significant (p<0.05): {n_sig}/{len(sp_df)}")
print(sp_df.head(10).to_string(index=False))
print(f"  Saved: step3c_spearman_table.tsv")

# ══════════════════════════════════════════════════════════════════════════════
# FIG 1: All tissues — horizontal bar ranked by -log10(p)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Figure 1: All tissues ────────────────────────────────────")

plot_df = sp_df.copy()
plot_df["log10p"] = -np.log10(plot_df["pval"].clip(lower=1e-10))
plot_df = plot_df[plot_df["log10p"] > 0.01]
plot_df = plot_df.sort_values("log10p", ascending=True).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(10, max(8, len(plot_df)*0.28)))

colors = []
for _, r in plot_df.iterrows():
    if r["cell_type"] == "Blood-Granul":
        colors.append(CB_HIGHLIGHT)
    elif r["pval"] < 0.05 and r["rho"] > 0:
        colors.append(CB_HIGHLIGHT)
    elif r["pval"] < 0.05 and r["rho"] < 0:
        colors.append("#0072B2")
    else:
        colors.append(CB_GREY)

ax.barh(plot_df["cell_type"], plot_df["log10p"], color=colors, edgecolor="none", height=0.7)

for _, r in plot_df[plot_df["log10p"] > -np.log10(0.2)].iterrows():
    ax.text(r["log10p"]+0.04, r["cell_type"], f"{r['pval']:.2e}",
            va="center", fontsize=8.5, color="#333")

ax.axvline(-np.log10(0.05), color="#2c3e50", linestyle="--", linewidth=1.5)
ax.text(-np.log10(0.05)+0.04, 1, "p = 0.05", color="#2c3e50", fontsize=9)

legend_handles = [
    mpatches.Patch(color=CB_HIGHLIGHT, label="Increases with age (p<0.05)"),
    mpatches.Patch(color="#0072B2", label="Decreases with age (p<0.05)"),
    mpatches.Patch(color=CB_GREY, label="Not significant"),
]
ax.legend(handles=legend_handles, fontsize=9, loc="lower right")
ax.set_xlabel("-log10(p)  [Spearman, continuous age]")
ax.set_title(f"Age-related changes in cfDNA tissue composition\n"
             f"(Loyfer 2023 atlas, {len(cell_cols)} cell types, n=35)")
plt.tight_layout()
plt.savefig("step3c_fig1_all_tissues.png"); plt.close()
print("  Saved: step3c_fig1_all_tissues.png")

# ══════════════════════════════════════════════════════════════════════════════
# FIG 2: Top cell types grouped bar by age
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Figure 2: Top cell types ─────────────────────────────────")

TOP_N = 8
top_cells = df[cell_cols].mean().sort_values(ascending=False).head(TOP_N).index.tolist()

grp_mean = df.groupby("age_group")[top_cells].mean().reindex(
    [g for g in GROUP_ORDER if g in df["age_group"].unique()])
grp_sem = df.groupby("age_group")[top_cells].sem().reindex(
    [g for g in GROUP_ORDER if g in df["age_group"].unique()])

x = np.arange(len(top_cells))
n_grp = len(grp_mean)
width = 0.22
offsets = np.linspace(-(n_grp-1)/2, (n_grp-1)/2, n_grp) * width

fig, ax = plt.subplots(figsize=(14, 6))

for i, (grp, row) in enumerate(grp_mean.iterrows()):
    vals = row.values * 100
    sems = grp_sem.loc[grp].values * 100
    bars = ax.bar(x+offsets[i], vals, width,
                  label=f"{grp}", color=AGE_COLORS.get(grp,"#aaa"),
                  edgecolor="white", linewidth=0.8, zorder=3)
    ax.errorbar(x+offsets[i], vals, yerr=sems, fmt="none",
                color="black", capsize=3, linewidth=1, zorder=4)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.3,
                f"{val:.1f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")

# Highlight granulocyte
if "Blood-Granul" in top_cells:
    idx = top_cells.index("Blood-Granul")
    ax.axvspan(idx-0.42, idx+0.42, alpha=0.06, color=CB_HIGHLIGHT, zorder=0)
    granu_p = sp_df[sp_df["cell_type"]=="Blood-Granul"]["pval"].values[0]
    ax.text(idx, grp_mean.values.max()*100*1.2,
            f"p = {granu_p:.2e}", ha="center", fontsize=9,
            color=CB_HIGHLIGHT, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(top_cells, rotation=25, ha="right", fontsize=10)
ax.set_ylabel("Proportion of cfDNA (%)")
ax.set_title(f"Top {TOP_N} cell types in cfDNA by age group")
ax.legend(fontsize=10, title="Age Group")
ax.set_ylim(0, grp_mean.values.max()*100*1.35)
plt.tight_layout()
plt.savefig("step3c_fig2_top_celltypes.png"); plt.close()
print("  Saved: step3c_fig2_top_celltypes.png")

# ══════════════════════════════════════════════════════════════════════════════
# FIG 3: Granulocyte scatter
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Figure 3: Granulocyte scatter ────────────────────────────")

CT = "Blood-Granul"
rho, p_sp = stats.spearmanr(df["age_numeric"], df[CT])
slope, intercept, r, pv, se = stats.linregress(df["age_numeric"], df[CT])

fig, ax = plt.subplots(figsize=(8, 6))

for grp in GROUP_ORDER:
    sub = df[df["age_group"]==grp]
    if len(sub)==0: continue
    ax.scatter(sub["age_numeric"], sub[CT]*100, color=AGE_COLORS[grp],
               s=90, alpha=0.85, label=f"{grp} (n={len(sub)})",
               edgecolors="white", linewidths=0.8, zorder=4)

x_line = np.linspace(df["age_numeric"].min()-3, df["age_numeric"].max()+3, 200)
y_line = (slope*x_line + intercept)*100
ax.plot(x_line, y_line, color="#2c3e50", linewidth=2, linestyle="--",
        label="Linear fit", zorder=3)

# 95% CI
n = len(df); x_m = df["age_numeric"].mean()
sxx = np.sum((df["age_numeric"]-x_m)**2)
ci = 1.96*se*np.sqrt(1/n+(x_line-x_m)**2/sxx)
ax.fill_between(x_line, (y_line-ci*100), (y_line+ci*100), alpha=0.12, color="#2c3e50")

young_m = df[df["age_group"]=="Young"][CT].mean()*100
old_m = df[df["age_group"]=="Old"][CT].mean()*100
fold = old_m/young_m

ax.annotate(
    f"Spearman rho = {rho:.3f}\np = {p_sp:.2e}\n"
    f"Young: {young_m:.1f}%\nOld: {old_m:.1f}%\nFold: {fold:.1f}x",
    xy=(0.05,0.97), xycoords="axes fraction", va="top", ha="left", fontsize=10,
    bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#ccc", alpha=0.95))

ax.set_xlabel("Age (years)")
ax.set_ylabel("Blood-Granulocyte proportion of cfDNA (%)")
ax.set_title(f"Granulocyte cfDNA fraction increases {fold:.1f}x with age\n"
             f"(Loyfer 2023 atlas, n=35)")
ax.legend(fontsize=10)
ax.grid(alpha=0.2); ax.set_axisbelow(True)
plt.savefig("step3c_fig3_granulocyte_scatter.png"); plt.close()
print("  Saved: step3c_fig3_granulocyte_scatter.png")

print("\n── ALL OUTPUTS ─────────────────────────────────────────────")
print("  step3c_fig1_all_tissues.png")
print("  step3c_fig2_top_celltypes.png")
print("  step3c_fig3_granulocyte_scatter.png")
print("  step3c_spearman_table.tsv")
