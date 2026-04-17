#!/usr/bin/env python
"""
STEP4B_chromhmm.py — ChromHMM Epigenomic Annotation


Annotates concordant vs discordant CpGs with chromatin states
from Roadmap Epigenomics 15-state ChromHMM model.

Runs on TWO input sets:
  (A) CGI CpGs only (cgi_concordant.bed, cgi_discordant.bed)
  (B) ALL significant CpGs (all_concordant.bed, all_discordant.bed)

TISSUES (Roadmap Epigenomics IDs):
  E062  Peripheral blood mononuclear cells (PBMC)
  E030  Primary neutrophils (CD15+) — i.e. granulocytes
  E029  Primary monocytes (CD14+)
  E034  Primary T cells
  E067  Brain angular gyrus
  E066  Liver
  E075  Colonic mucosa


INPUT:  BED files from STEP4_prep.py
OUTPUT: Per-set summary tables, figures, PCDHGA annotation
"""

import os, subprocess, gzip, shutil, urllib.request
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ══════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════

ROADMAP_BASE = ("https://egg2.wustl.edu/roadmap/data/byFileType/"
                "chromhmmSegmentations/ChmmModels/coreMarks/"
                "jointModel/final/")

TISSUES = {
    "Blood_PBMC":       "E062",
    "Neutrophil_CD15":  "E030",
    "Monocyte_CD14":    "E029",
    "T_cell":           "E034",
    "Brain_AG":         "E067",
    "Liver":            "E066",
    "Colon":            "E075",
}

STATE_GROUPS = {
    "1_TssA": "Active promoter", "2_TssAFlnk": "Active promoter",
    "3_TxFlnk": "Tran- scription", "4_Tx": "Tran- scription",
    "5_TxWk": "Tran- scription",
    "6_EnhG": "Enhancer", "7_Enh": "Enhancer",
    "8_ZNF/Rpts": "Other", "9_Het": "Hetero- chromatin",
    "10_TssBiv": "Bivalent", "11_BivFlnk": "Bivalent",
    "12_EnhBiv": "Bivalent",
    "13_ReprPC": "Polycomb repressed", "14_ReprPCWk": "Polycomb repressed",
    "15_Quies": "Quiescent",
}

GROUP_ORDER = ["Active promoter", "Enhancer", "Tran- scription",
               "Bivalent", "Polycomb repressed", "Quiescent", "Hetero- chromatin",
                "Other"]

# Okabe-Ito
OI = {"concordant": "#009E73", "discordant": "#E69F00"}

# ══════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════

def download_chromhmm(eid, tissue_name, outdir="chromhmm_beds"):
    os.makedirs(outdir, exist_ok=True)
    fname = f"{eid}_15_coreMarks_mnemonics.bed.gz"
    url = ROADMAP_BASE + fname
    gz_path = os.path.join(outdir, fname)
    bed_path = gz_path.replace(".gz", "")
    if os.path.exists(bed_path):
        print(f"  {tissue_name} ({eid}): cached")
        return bed_path
    print(f"  {tissue_name} ({eid}): downloading...")
    urllib.request.urlretrieve(url, gz_path)
    with gzip.open(gz_path, "rb") as fi, open(bed_path, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(gz_path)
    return bed_path


def intersect(cpg_bed, chromhmm_bed):
    """Intersect CpG BED with ChromHMM BED. Tries bedtools, falls back to pandas."""
    try:
        r = subprocess.run(
            ["bedtools", "intersect", "-a", cpg_bed, "-b", chromhmm_bed, "-wa", "-wb"],
            capture_output=True, text=True, check=True)
        out = []
        for line in r.stdout.strip().split("\n"):
            if not line: continue
            f = line.split("\t")
            out.append((f[0], int(f[1]), f[6]))
        return out
    except (FileNotFoundError, subprocess.CalledProcessError):
        return intersect_pandas(cpg_bed, chromhmm_bed)


def intersect_pandas(cpg_bed, chromhmm_bed):
    print("    (pandas fallback — slower)")
    cpgs = pd.read_csv(cpg_bed, sep="\t", header=None, names=["chrom","start","end"])
    hmm = pd.read_csv(chromhmm_bed, sep="\t", header=None,
                       names=["chrom","hmm_start","hmm_end","state"])
    out = []
    for _, row in cpgs.iterrows():
        m = hmm[(hmm["chrom"]==row["chrom"]) &
                (hmm["hmm_start"]<=row["start"]) &
                (hmm["hmm_end"]>row["start"])]
        state = m.iloc[0]["state"] if len(m) > 0 else "NA"
        out.append((row["chrom"], row["start"], state))
    return out


def annotate_set(conc_bed, disc_bed, tissue_beds, set_label):
    """Annotate a concordant/discordant BED pair across all tissues."""
    rows = []
    for label, bed_file in [("concordant", conc_bed), ("discordant", disc_bed)]:
        if not os.path.exists(bed_file):
            print(f"  SKIPPED: {bed_file} not found")
            continue
        for tissue, bed_path in tissue_beds.items():
            states = intersect(bed_file, bed_path)
            for chrom, pos, state in states:
                rows.append({
                    "chrom": chrom, "position": pos,
                    "cpg_set": label, "tissue": tissue,
                    "state_raw": state,
                    "state_group": STATE_GROUPS.get(state, "Other"),
                    "input_set": set_label,
                })
    return pd.DataFrame(rows)


def make_summary(df):
    """Compute % per state_group per tissue per cpg_set."""
    rows = []
    for tissue in df["tissue"].unique():
        for cpg_set in ["concordant", "discordant"]:
            sub = df[(df["tissue"]==tissue) & (df["cpg_set"]==cpg_set)]
            total = len(sub)
            for grp in GROUP_ORDER:
                n = (sub["state_group"]==grp).sum()
                pct = round(100*n/total, 1) if total > 0 else 0
                rows.append({"tissue": tissue, "cpg_set": cpg_set,
                             "state_group": grp, "n": n, "pct": pct})
    return pd.DataFrame(rows)


def make_figure(summary, tissues_to_plot, title_suffix, outfile):
    """Grouped bar chart of state distributions."""
    n_t = len(tissues_to_plot)
    fig, axes = plt.subplots(1, n_t, figsize=(3.5*n_t, 5), sharey=True)
    if n_t == 1: axes = [axes]
    bw = 0.35
    for ax, tissue in zip(axes, tissues_to_plot):
        cp, dp = [], []
        for grp in GROUP_ORDER:
            c = summary[(summary["tissue"]==tissue)&(summary["cpg_set"]=="concordant")&(summary["state_group"]==grp)]
            d = summary[(summary["tissue"]==tissue)&(summary["cpg_set"]=="discordant")&(summary["state_group"]==grp)]
            cp.append(c["pct"].values[0] if len(c) > 0 else 0)
            dp.append(d["pct"].values[0] if len(d) > 0 else 0)
        x = np.arange(len(GROUP_ORDER))
        ax.bar(x-bw/2, cp, bw, color=OI["concordant"], label="Concordant", alpha=0.85)
        ax.bar(x+bw/2, dp, bw, color=OI["discordant"], label="Discordant", alpha=0.85)
        ax.set_title(tissue.replace("_"," "), fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([g.replace(" ","\n") for g in GROUP_ORDER], fontsize=5.5, ha="center")
        ax.set_ylim(0, 70)
        if ax == axes[0]:
            ax.set_ylabel("% of CpGs", fontsize=11)
            ax.legend(fontsize=8)
    fig.suptitle(f"Chromatin state context — {title_suffix}\n(Roadmap 15-state ChromHMM)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile}")


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ── Download ─────────────────────────────────────────────
    print("Downloading ChromHMM files...")
    tissue_beds = {}
    for name, eid in TISSUES.items():
        tissue_beds[name] = download_chromhmm(eid, name)

    # ── Define input sets ────────────────────────────────────
    sets = {
        "CGI": ("cgi_concordant.bed", "cgi_discordant.bed"),
        "ALL": ("all_concordant.bed", "all_discordant.bed"),
    }

    all_dfs = []

    for set_label, (conc_bed, disc_bed) in sets.items():
        print(f"\n{'═'*60}")
        print(f"Annotating {set_label} set...")
        print(f"{'═'*60}")
        df = annotate_set(conc_bed, disc_bed, tissue_beds, set_label)
        if len(df) == 0:
            print(f"  No results for {set_label} — skipping")
            continue
        all_dfs.append(df)

        # Per-CpG annotations
        outf = f"chromhmm_{set_label.lower()}_state_counts.tsv"
        df.to_csv(outf, sep="\t", index=False)
        print(f"  Saved: {outf}")

        # Summary
        summary = make_summary(df)
        outf = f"chromhmm_{set_label.lower()}_summary.tsv"
        summary.to_csv(outf, sep="\t", index=False)
        print(f"  Saved: {outf}")

        # Main figure (blood + brain + liver + colon)
        main_tissues = [t for t in ["Blood_PBMC", "Brain_AG", "Liver", "Colon"]
                        if t in df["tissue"].unique()]
        make_figure(summary, main_tissues,
                    f"{set_label} CpGs", f"chromhmm_{set_label.lower()}_fig_main.png")

        # Blood-only figure (PBMC + neutrophil + monocyte + T cell)
        blood_tissues = [t for t in ["Blood_PBMC", "Neutrophil_CD15", "Monocyte_CD14", "T_cell"]
                         if t in df["tissue"].unique()]
        if len(blood_tissues) > 1:
            make_figure(summary, blood_tissues,
                        f"{set_label} CpGs — blood subtypes",
                        f"chromhmm_{set_label.lower()}_fig_blood.png")

        # Print Polycomb comparison
        print(f"\n  Polycomb repressed (%) — {set_label}:")
        for tissue in df["tissue"].unique():
            for cs in ["concordant", "discordant"]:
                row = summary[(summary["tissue"]==tissue)&
                              (summary["cpg_set"]==cs)&
                              (summary["state_group"]=="Polycomb repressed")]
                pct = row["pct"].values[0] if len(row) > 0 else "?"
                print(f"    {tissue:20s} {cs:12s} {pct}%")

    # ── PCDHGA-specific ──────────────────────────────────────
    print(f"\n{'═'*60}")
    print("PCDHGA CpG (chr5:140821425)")
    print(f"{'═'*60}")
    if len(all_dfs) > 0:
        combined = pd.concat(all_dfs, ignore_index=True)
        pcdhga = combined[combined["position"] == 140821425]
        if len(pcdhga) > 0:
            for _, row in pcdhga.drop_duplicates(["tissue","state_raw"]).iterrows():
                print(f"  {row['tissue']:20s}  {row['state_raw']:20s}  ({row['state_group']})")
            pcdhga.to_csv("chromhmm_pcdhga_states.tsv", sep="\t", index=False)
            print("  Saved: chromhmm_pcdhga_states.tsv")
        else:
            print("  Not found in annotated CpGs")

    print(f"\n{'═'*60}")
    print("Done.")
    print(f"{'═'*60}")
