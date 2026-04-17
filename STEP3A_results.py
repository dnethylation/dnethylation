#!/usr/bin/env python
"""
STEP3A_results.py — CpG island + gene feature analysis.

INPUT:
    step2_significant_overlap.tsv           (from STEP2_results.py)
    cpgIslandExt.txt                        (UCSC hg19 track)
    infinium-methylationepic-v-1-0-b5-manifest-file.csv  (Illumina)

OUTPUT:
    cgi_annotated_full.tsv                  (used by STEP3B_results.py, STEP4_prep.py)
    step3a_enrichment_table.tsv             (supplementary)
    step3a_fig1_cgi_concordance.png         (Figure 7)
    step3a_fig2_feature_enrichment.png      (Figure 8A / 8B)
    step3a_fig3_direction_by_feature.png    (Figure 9)
    
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt
from scipy.stats import binomtest

plt.rcParams.update({"font.family":"Arial","font.size":11,"axes.titlesize":13,
    "axes.titleweight":"bold","axes.labelsize":12,"axes.spines.top":False,
    "axes.spines.right":False,"axes.grid":False,"figure.dpi":200,
    "savefig.dpi":300,"savefig.bbox":"tight"})

CB_BLUE="#0072B2";CB_ORANGE="#E69F00";CB_TEAL="#009E73";CB_GREY="#999999";CB_BLACK="#000000"
FEAT_ORDER=["TSS200","TSS1500","5'UTR","1stExon","Body","3'UTR","Intergenic"]
FEAT_COLORS=["#c0392b","#e67e22","#f39c12","#d4ac0d","#27ae60","#2980b9","#7f8c8d"]

def wilson_ci(k,n,z=1.96):
    p=k/n;d=1+z**2/n;c=(p+z**2/(2*n))/d;m=z*np.sqrt(p*(1-p)/n+z**2/(4*n**2))/d
    return 100*(c-m),100*(c+m)

def primary_feature(val):
    if pd.isna(val) or str(val).strip()=="": return "Intergenic"
    parts=[p.strip() for p in str(val).split(";")]
    for f in ["TSS200","TSS1500","5'UTR","1stExon","Body","3'UTR"]:
        if f in parts: return f
    return parts[0]

def clean_gene(val):
    if pd.isna(val) or str(val).strip()=="": return ""
    return ";".join(sorted(set(p.strip() for p in str(val).split(";") if p.strip())))

# ── LOAD ──────────────────────────────────────────────────────────────────────
sig = pd.read_csv("step2_significant_overlap.tsv", sep="\t")

# ── CGI TAGGING (was in step2, now here) ──────────────────────────────────────
cgi = pd.read_csv("cpgIslandExt.txt", sep="\t", header=None, usecols=[1,2,3], names=["chrom","start","end"])
sig = sig.sort_values(["chrom","position"]).reset_index(drop=True)
in_cgi = np.zeros(len(sig), dtype=bool)
for chrom, grp in sig.groupby("chrom"):
    islands = cgi[cgi["chrom"]==chrom].sort_values("start")
    if len(islands)==0: continue
    pos = grp["position"].values
    idx = np.searchsorted(islands["start"].values, pos, side="right") - 1
    valid = idx >= 0
    in_cgi[grp.index] = valid & (pos <= islands["end"].values[np.where(valid, idx, 0)])
sig["in_cgi"] = in_cgi

cgi_sig = sig[sig["in_cgi"]==True].copy().reset_index(drop=True)
noncgi_sig = sig[sig["in_cgi"]==False].copy().reset_index(drop=True)
print(f"CGI: {len(cgi_sig):,}  |  Non-CGI: {len(noncgi_sig):,}")

# ── MANIFEST ANNOTATION ──────────────────────────────────────────────────────
manifest = pd.read_csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv",
    skiprows=7, low_memory=False,
    usecols=["Name","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","450k_Enhancer"])
manifest = manifest[manifest["Name"].str.startswith("cg",na=False)].copy()
manifest = manifest.rename(columns={"Name":"cpg","UCSC_RefGene_Name":"manifest_gene",
    "UCSC_RefGene_Group":"gene_feature_raw","Relation_to_UCSC_CpG_Island":"island_relation",
    "450k_Enhancer":"is_enhancer"})
manifest["island_relation"] = manifest["island_relation"].fillna("Island").replace("","Island")

cgi_ann = cgi_sig.merge(manifest, on="cpg", how="left")
cgi_ann["primary_feature"] = cgi_ann["gene_feature_raw"].apply(primary_feature).replace({"5UTR":"5'UTR","3UTR":"3'UTR"})
cgi_ann["gene_clean"] = cgi_ann["manifest_gene"].apply(clean_gene)
cgi_ann["is_enhancer_bool"] = cgi_ann["is_enhancer"].notna() & (cgi_ann["is_enhancer"]!="")
cgi_ann.to_csv("cgi_annotated_full.tsv", sep="\t", index=False)

n_total=len(cgi_ann); n_conc=cgi_ann["concordant"].sum(); n_disc=n_total-n_conc
print(f"Annotated: {n_total} CGI CpGs ({n_conc} concordant, {n_disc} discordant)")

# ── FIG 1: CGI vs Non-CGI concordance ─────────────────────────────────────────
fig,ax=plt.subplots(figsize=(9,5.5))
groups=[]
for label,subset in [("CGI\nAll",cgi_sig),("CGI\nHypo",cgi_sig[cgi_sig["direction"]=="hypo"]),
    ("CGI\nHyper",cgi_sig[cgi_sig["direction"]=="hyper"]),
    ("Non-CGI\nAll",noncgi_sig),("Non-CGI\nHypo",noncgi_sig[noncgi_sig["direction"]=="hypo"]),
    ("Non-CGI\nHyper",noncgi_sig[noncgi_sig["direction"]=="hyper"])]:
    k=subset["concordant"].sum();n=len(subset);pct=100*k/n if n>0 else 0
    lo,hi=wilson_ci(k,n) if n>0 else (0,0)
    groups.append({"label":label,"pct":pct,"n":n,"lo":lo,"hi":hi})
g=pd.DataFrame(groups);colors=[CB_TEAL,CB_BLUE,CB_ORANGE]*2
bars=ax.bar(range(len(g)),g["pct"],color=colors,width=0.65,edgecolor="white",linewidth=1.5,zorder=3)
for i,(bar,row) in enumerate(zip(bars,g.itertuples())):
    ax.errorbar(i,row.pct,yerr=[[row.pct-row.lo],[row.hi-row.pct]],fmt="none",color=CB_BLACK,capsize=5,linewidth=1.5,zorder=4)
    ax.text(i,row.hi+1,f"{row.pct:.1f}%\n(n={row.n:,})",ha="center",va="bottom",fontsize=9,fontweight="bold")
ax.set_xticks(range(len(g)));ax.set_xticklabels(g["label"],fontsize=9)
ax.axhline(50,color=CB_GREY,linestyle="--",linewidth=1.5,label="50% (chance)")
ax.axvline(2.5,color=CB_GREY,linestyle=":",linewidth=1,alpha=0.5)
ax.set_ylabel("Concordance with blood EWAS (%)");ax.set_ylim(25,80);ax.legend(loc="upper right")
ax.set_title("Concordance by CpG island status and direction")
plt.savefig("step3a_fig1_cgi_concordance.png");plt.close()

# ── FIG 2: enrichment + concordance per feature ──────────────────────────────
bg_cgi=manifest[manifest["island_relation"]=="Island"].copy()
bg_cgi["feature"]=bg_cgi["gene_feature_raw"].apply(primary_feature).replace({"5UTR":"5'UTR","3UTR":"3'UTR"})
bg_pct=bg_cgi["feature"].value_counts(normalize=True).reindex(FEAT_ORDER,fill_value=0)*100
ag_pct=cgi_ann["primary_feature"].value_counts(normalize=True).reindex(FEAT_ORDER,fill_value=0)*100

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5.5))
x=np.arange(len(FEAT_ORDER));w=0.35
ax1.bar(x-w/2,bg_pct.values,w,label="EPIC CGI (background)",color=CB_GREY,edgecolor="white")
ax1.bar(x+w/2,ag_pct.values,w,label="Aging CGI CpGs",color=CB_TEAL,edgecolor="white")
for i,(b,a) in enumerate(zip(bg_pct.values,ag_pct.values)):
    d=a-b
    if abs(d)>2: ax1.text(i+w/2,a+1,f"{'+' if d>0 else ''}{d:.1f}%",ha="center",fontsize=8,fontweight="bold",color=CB_TEAL if d>0 else CB_ORANGE)
ax1.set_xticks(x);ax1.set_xticklabels(FEAT_ORDER,fontsize=9,rotation=15,ha="right")
ax1.set_ylabel("% of CpGs");ax1.set_title("A. Gene feature distribution\nvs EPIC array background");ax1.legend(fontsize=9)

crs,clo,chi,ns=[],[],[],[]
for f in FEAT_ORDER:
    sub=cgi_ann[cgi_ann["primary_feature"]==f];n=len(sub);ns.append(n)
    if n==0: crs.append(np.nan);clo.append(0);chi.append(0);continue
    k=sub["concordant"].sum();p=100*k/n;lo,hi=wilson_ci(k,n);crs.append(p);clo.append(lo);chi.append(hi)
bars2=ax2.bar(FEAT_ORDER,crs,color=FEAT_COLORS,width=0.6,edgecolor="white",zorder=3)
for i,(p,lo,hi,n) in enumerate(zip(crs,clo,chi,ns)):
    if np.isnan(p):continue
    ax2.errorbar(i,p,yerr=[[p-lo],[hi-p]],fmt="none",color=CB_BLACK,capsize=5,linewidth=1.5,zorder=4)
    ax2.text(i,hi+1.5,f"n={n}",ha="center",fontsize=8,color="#555")
ax2.axhline(50,color=CB_GREY,linestyle="--",linewidth=1.5)
ax2.set_ylabel("Concordance (%)");ax2.set_title("B. Concordance per gene feature\n(CGI CpGs, 95% Wilson CI)");ax2.set_ylim(20,95)
ax2.set_xticklabels(FEAT_ORDER,fontsize=9,rotation=15,ha="right")
plt.tight_layout(w_pad=3);plt.savefig("step3a_fig2_feature_enrichment.png");plt.close()

pd.DataFrame({"feature":FEAT_ORDER,"EPIC_bg_pct":bg_pct.values.round(1),"aging_CGI_pct":ag_pct.values.round(1),
    "difference":(ag_pct.values-bg_pct.values).round(1),"concordance_pct":[round(c,1) if not np.isnan(c) else "-" for c in crs],"N":ns
}).to_csv("step3a_enrichment_table.tsv",sep="\t",index=False)

# ── FIG 3: direction split ────────────────────────────────────────────────────
conc_only=cgi_ann[cgi_ann["concordant"]==True]
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5))
for ax,subset,title in [(ax1,cgi_ann,f"All CGI (n={n_total})"),(ax2,conc_only,f"Concordant (n={n_conc})")]:
    hp,yp,ts=[],[],[]
    for f in FEAT_ORDER:
        sub=subset[subset["primary_feature"]==f];t=len(sub);h=(sub["direction"]=="hypo").sum()
        ts.append(t);hp.append(100*h/t if t>0 else 0);yp.append(100*(t-h)/t if t>0 else 0)
    x=np.arange(len(FEAT_ORDER))
    ax.bar(x,hp,color=CB_BLUE,edgecolor="white",label="Hypo")
    ax.bar(x,yp,bottom=hp,color=CB_ORANGE,edgecolor="white",label="Hyper")
    for i,t in enumerate(ts):
        if t>0:ax.text(i,50,f"n={t}",ha="center",va="center",fontsize=8,color="white",fontweight="bold")
    ax.axhline(50,color="white",linestyle="--",linewidth=1,alpha=0.7)
    ax.set_xticks(x);ax.set_xticklabels(FEAT_ORDER,fontsize=9,rotation=15,ha="right")
    ax.set_ylabel("% of CpGs");ax.set_ylim(0,100);ax.set_title(title,fontweight="bold");ax.legend(fontsize=9,loc="lower right")
plt.suptitle("Direction of age-association by gene feature",fontweight="bold",y=1.01)
plt.tight_layout();plt.savefig("step3a_fig3_direction_by_feature.png");plt.close()

print(f"Done. cgi_annotated_full.tsv ({n_total} rows) + 3 figures + enrichment table")
