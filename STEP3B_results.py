#!/usr/bin/env python
"""
STEP3B_results.py — gene-level analysis, epigenetic clock overlap, BED export.

INPUT:
    cgi_annotated_full.tsv                  (from STEP3A_results.py)
    step2_significant_overlap.tsv           (from STEP2_results.py)
    gb-2013-14-10-r115-S3.csv               (Horvath clock)
    NIHMS418935-supplement-03.csv           (Hannum clock)

OUTPUT:
    step3b_gene_table.tsv                   (supplementary; underlies Figure 11)
    step3b_clock_overlap.tsv                (supplementary; underlies clock paragraph in §3.4)
    step3b_literature_table.tsv             (curated manually — displayed inline in §3.4)
    step3b_fig1_top_genes.png               (Figure 11)
    step3b_fig2_clock_overlap.png           (referenced inline in §3.4)
    great_all_cgi.bed                       (for GREAT web submission — §2.6.5)
    great_concordant_cgi.bed                (for GREAT web submission — §2.6.5)
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({"font.family":"Arial","font.size":11,"axes.titlesize":13,
    "axes.titleweight":"bold","axes.labelsize":12,"axes.spines.top":False,
    "axes.spines.right":False,"figure.dpi":200,"savefig.dpi":300,"savefig.bbox":"tight"})

CB_BLUE="#0072B2";CB_ORANGE="#E69F00";CB_TEAL="#009E73";CB_GREY="#999999"

df = pd.read_csv("cgi_annotated_full.tsv", sep="\t")
conc=df[df["concordant"]==True]; disc=df[df["concordant"]==False]
n_t=len(df);n_c=len(conc);n_d=len(disc)

# ── 1. GENE-LEVEL CONSISTENCY ─────────────────────────────────────────────────
rows=[]
for _,r in df.iterrows():
    genes=str(r["gene_clean"]).split(";") if r["gene_clean"] else [""]
    for g in genes:
        g=g.strip()
        if g: rows.append({"gene":g,"concordant":r["concordant"],"direction":r["direction"],"cpg":r["cpg"]})
gene_df=pd.DataFrame(rows)

gene_summary=gene_df.groupby("gene").agg(
    n_cpgs=("cpg","nunique"),n_concordant=("concordant","sum"),
    directions=("direction",lambda x:";".join(sorted(set(x)))),
    cpg_ids=("cpg",lambda x:";".join(sorted(set(x))))
).reset_index()
gene_summary["n_discordant"]=gene_summary["n_cpgs"]-gene_summary["n_concordant"]
gene_summary["pct_concordant"]=(100*gene_summary["n_concordant"]/gene_summary["n_cpgs"]).round(1)
gene_summary["consistency"]=gene_summary.apply(
    lambda r:"All concordant" if r["n_concordant"]==r["n_cpgs"] else "All discordant" if r["n_concordant"]==0 else "Mixed",axis=1)
gene_summary=gene_summary.sort_values("n_cpgs",ascending=False).reset_index(drop=True)
gene_summary.to_csv("step3b_gene_table.tsv",sep="\t",index=False)

multi=gene_summary[gene_summary["n_cpgs"]>=2]
print(f"Genes >=2 CpGs: {len(multi)} | All conc: {(multi['consistency']=='All concordant').sum()} | Mixed: {(multi['consistency']=='Mixed').sum()} | All disc: {(multi['consistency']=='All discordant').sum()}")

# ── FIG 1: top genes ──────────────────────────────────────────────────────────
def top_genes(subset,n=15):
    return subset["gene_clean"].dropna().str.split(";").explode().str.strip().replace("",pd.NA).dropna().value_counts().head(n)

top_c=top_genes(conc);top_d=top_genes(disc)
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,6))
ax1.barh(top_c.index[::-1],top_c.values[::-1],color=CB_TEAL,edgecolor="white")
for i,(g,v) in enumerate(zip(top_c.index[::-1],top_c.values[::-1])): ax1.text(v+0.05,i,str(v),va="center",fontsize=9,fontweight="bold")
ax1.set_xlabel("Number of CpGs");ax1.set_title(f"Concordant genes ({n_c} CpGs)",color=CB_TEAL);ax1.set_xlim(0,top_c.values.max()*1.3)

ax2.barh(top_d.index[::-1],top_d.values[::-1],color=CB_ORANGE,edgecolor="white")
for i,(g,v) in enumerate(zip(top_d.index[::-1],top_d.values[::-1])): ax2.text(v+0.05,i,str(v),va="center",fontsize=9,fontweight="bold")
ax2.set_xlabel("Number of CpGs");ax2.set_title(f"Discordant genes ({n_d} CpGs)",color=CB_ORANGE);ax2.set_xlim(0,top_d.values.max()*1.3)
for label in ax2.get_yticklabels():
    if "PCDHGA" in label.get_text(): label.set_fontweight("bold");label.set_color("#7b241c")
fig.suptitle("Gene-level analysis of CGI CpGs",fontweight="bold")
plt.tight_layout();plt.savefig("step3b_fig1_top_genes.png");plt.close()

# ── 2. CLOCK CROSS-REFERENCE ─────────────────────────────────────────────────
all_sig=pd.read_csv("step2_significant_overlap.tsv",sep="\t")
clock_results=[]

for cname,cfile,cpg_col,skip in [
    ("Horvath","gb-2013-14-10-r115-S3.csv","CpGmarker",2),
    ("Hannum","NIHMS418935-supplement-03.csv","Marker",0)]:
    if not Path(cfile).exists():
        print(f"  {cfile} not found — skipping {cname}")
        continue
    try:
        cdf=pd.read_csv(cfile,skiprows=skip)
        cpgs=set(cdf[cpg_col].dropna().astype(str).str.strip())
        cpgs.discard("(Intercept)")
    except Exception as e:
        print(f"  Error loading {cfile}: {e}"); continue

    print(f"{cname}: {len(cpgs)} clock CpGs")
    ov=set(all_sig["cpg"].dropna())&cpgs
    print(f"  vs all significant ({len(all_sig)}): {len(ov)} overlap")
    if ov:
        ov_df=all_sig[all_sig["cpg"].isin(ov)]
        nc=ov_df["concordant"].sum();nd=len(ov_df)-nc
        print(f"  Concordant: {nc}/{len(ov_df)}")
        for _,r in ov_df.iterrows():
            clock_results.append({"clock":cname,"cpg":r["cpg"],"direction":r["direction"],
                "concordant":r["concordant"],"in_cgi":r.get("in_cgi",False)})

if clock_results:
    crdf=pd.DataFrame(clock_results)
    crdf.to_csv("step3b_clock_overlap.tsv",sep="\t",index=False)

    fig,ax=plt.subplots(figsize=(8,4))
    for i,ck in enumerate(crdf["clock"].unique()):
        sub=crdf[crdf["clock"]==ck];nc=sub["concordant"].sum();nd=len(sub)-nc
        ax.barh(i-0.15,nc,0.3,color=CB_TEAL,edgecolor="white",label="Concordant" if i==0 else "")
        ax.barh(i+0.15,nd,0.3,color=CB_ORANGE,edgecolor="white",label="Discordant" if i==0 else "")
        ax.text(max(nc,0.1)+0.1,i-0.15,str(nc),va="center",fontsize=10,fontweight="bold")
        ax.text(max(nd,0.1)+0.1,i+0.15,str(nd),va="center",fontsize=10,fontweight="bold")
    ax.set_yticks(range(len(crdf["clock"].unique())))
    ax.set_yticklabels([f"{c} clock" for c in crdf["clock"].unique()],fontsize=11)
    ax.set_xlabel("Overlapping CpGs");ax.set_title("Overlap with epigenetic clocks\n(all significant CpGs, p<0.05)")
    ax.legend(loc="lower right");plt.savefig("step3b_fig2_clock_overlap.png");plt.close()

# ── 3. BED FILES ──────────────────────────────────────────────────────────────
for label,subset,fname in [("All CGI",df,"great_all_cgi.bed"),("Concordant",conc,"great_concordant_cgi.bed")]:
    bed=subset[["chrom","position"]].copy();bed["end"]=bed["position"]+1;bed["name"]="."
    bed[["chrom","position","end","name"]].to_csv(fname,sep="\t",header=False,index=False)
    print(f"BED: {label} -> {fname} ({len(bed)} regions)")

# ── 4. LITERATURE TABLE ──────────────────────────────────────────────────────
pd.DataFrame([
    {"gene":"RUNX1","status":"Concordant","function":"Haematopoietic TF; age-related hypermethylation reflects myeloid skewing","ref":"Chen 2025"},
    {"gene":"FAM19A4","status":"Concordant","function":"Tumour suppressor silenced by methylation","ref":"Kremer 2019"},
    {"gene":"LINC00857","status":"Concordant","function":"lncRNA; aging role not established","ref":"Chen 2022"},
    {"gene":"SOX11","status":"Concordant","function":"Developmental TF silenced in adult tissue","ref":"Gustavsson 2010"},
    {"gene":"PCDHGA1-12","status":"Discordant","function":"Most conserved cross-tissue aging locus (17 tissues)","ref":"Jacques 2025"},
    {"gene":"FOXR1","status":"Discordant","function":"Neural expression; neuroblastoma context","ref":"Waxman 2025"},
    {"gene":"SLC6A16","status":"Discordant","function":"Neurotransmitter transporter","ref":"Jomura 2022"},
]).to_csv("step3b_literature_table.tsv",sep="\t",index=False)

print(f"\nDone. gene_table + top_genes fig + clock overlap + BED files + literature table")
