#!/usr/bin/env python
"""
STEP2_results.py — concordance figures and tables.

INPUT:
    overlap_chinese+scottish.tsv            (from STEP2.py)

OUTPUT:
    step2_significant_overlap.tsv           (used by STEP3A_results.py, STEP4_prep.py)
    step2_concordance_table.tsv             (supplementary table)
    step2_pipeline_summary.tsv              (internal log; not cited in dissertation)
    step2_fig1_concordance_overview.png     (Figure 4)
    step2_fig2_stringency_effectsize.png    (Figure 5A / 5B)
    step2_fig3_chromosome_concordance.png   (Figure 6)
    
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binomtest
from scipy.ndimage import uniform_filter1d

plt.rcParams.update({
    "font.family": "Arial", "font.size": 11, "axes.titlesize": 13,
    "axes.titleweight": "bold", "axes.labelsize": 12,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
})

CB_BLUE="#0072B2"; CB_ORANGE="#E69F00"; CB_TEAL="#009E73"
CB_GREY="#999999"; CB_PURPLE="#CC79A7"; CB_BLACK="#000000"

def wilson_ci(k, n, z=1.96):
    p=k/n; d=1+z**2/n; c=(p+z**2/(2*n))/d; m=z*np.sqrt(p*(1-p)/n+z**2/(4*n**2))/d
    return 100*(c-m), 100*(c+m)

def fmt_p(p):
    if p<1e-10: return "p < 1e-10"
    elif p<0.001: return f"p = {p:.1e}"
    elif p<0.01: return f"p = {p:.4f}"
    else: return f"p = {p:.3f}"

df = pd.read_csv("overlap_chinese+scottish.tsv", sep="\t")
sig = df[df["pval"] < 0.05]
hypo_sig = sig[sig["direction"]=="hypo"]; hyper_sig = sig[sig["direction"]=="hyper"]
n_sig=len(sig); n_conc=sig["concordant"].sum()
hypo_conc=hypo_sig["concordant"].sum(); hyper_conc=hyper_sig["concordant"].sum()
binom_all=binomtest(n_conc,n_sig,0.5,alternative="greater")
binom_hypo=binomtest(hypo_conc,len(hypo_sig),0.5,alternative="greater")
binom_hyper=binomtest(hyper_conc,len(hyper_sig),0.5,alternative="greater")

pd.DataFrame([
    {"Stage":"Overlapping CpGs","N":len(df)},
    {"Stage":"Significant in cfDNA (p<0.05)","N":n_sig},
    {"Stage":"  Hypomethylating","N":len(hypo_sig)},
    {"Stage":"  Hypermethylating","N":len(hyper_sig)},
    {"Stage":"Concordant","N":n_conc},
    {"Stage":"  Hypo concordant","N":hypo_conc},
    {"Stage":"  Hyper concordant","N":hyper_conc},
    {"Stage":"Discordant","N":n_sig-n_conc},
]).to_csv("step2_pipeline_summary.tsv",sep="\t",index=False)

# FIG 1: concordance overview
fig,ax=plt.subplots(figsize=(7,5.5))
cats=["All CpGs","Hypomethylating","Hypermethylating"]
pcts=[100*n_conc/n_sig,100*hypo_conc/len(hypo_sig),100*hyper_conc/len(hyper_sig)]
ks=[n_conc,hypo_conc,hyper_conc]; ns=[n_sig,len(hypo_sig),len(hyper_sig)]
pvs=[binom_all.pvalue,binom_hypo.pvalue,binom_hyper.pvalue]; cols=[CB_TEAL,CB_BLUE,CB_ORANGE]
bars=ax.bar(cats,pcts,color=cols,width=0.55,edgecolor="white",linewidth=1.5,zorder=3)
for bar,pct,n,k,pv in zip(bars,pcts,ns,ks,pvs):
    lo,hi=wilson_ci(k,n); cx=bar.get_x()+bar.get_width()/2
    ax.errorbar(cx,pct,yerr=[[pct-lo],[hi-pct]],fmt="none",color=CB_BLACK,capsize=6,linewidth=1.5,zorder=4)
    ax.text(cx,hi+1.2,f"{pct:.1f}%\n(n={n:,})",ha="center",va="bottom",fontsize=10,fontweight="bold")
    ax.text(cx,31.5,fmt_p(pv),ha="center",va="bottom",fontsize=9,color="#222",style="italic")
ax.axhline(50,color=CB_GREY,linestyle="--",linewidth=1.5,label="50% (chance)")
ax.set_ylabel("Concordance with blood EWAS (%)"); ax.set_ylim(30,72)
ax.legend(loc="upper right")
ax.set_title("Directional concordance between cfDNA and blood\nage-associated CpGs (cfDNA p < 0.05)")
plt.savefig("step2_fig1_concordance_overview.png"); plt.close()

# FIG 2: stringency + effect size
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.5))
thresholds=[0.05,0.01,0.005,0.001,0.0005,0.0001]
sa,sh,sy,sn=[],[],[],[]
for t in thresholds:
    sub=df[df["pval"]<t]
    if len(sub)<5: sa.append(np.nan);sh.append(np.nan);sy.append(np.nan);sn.append(0);continue
    sa.append(100*sub["concordant"].mean()); sn.append(len(sub))
    s1=sub[sub["direction"]=="hypo"]; s2=sub[sub["direction"]=="hyper"]
    sh.append(100*s1["concordant"].mean() if len(s1)>=5 else np.nan)
    sy.append(100*s2["concordant"].mean() if len(s2)>=5 else np.nan)
ax1.plot(thresholds,sa,"o-",color=CB_TEAL,linewidth=2.2,markersize=8,label="All",zorder=4)
ax1.plot(thresholds,sh,"s--",color=CB_BLUE,linewidth=1.5,markersize=6,label="Hypo",zorder=3)
ax1.plot(thresholds,sy,"^--",color=CB_ORANGE,linewidth=1.5,markersize=6,label="Hyper",zorder=3)
ax1.axhline(50,color=CB_GREY,linestyle="--",linewidth=1.2)
for xp,n,pct in zip(thresholds,sn,sa):
    if n>0 and not np.isnan(pct): ax1.annotate(f"n={n:,}",(xp,pct),textcoords="offset points",xytext=(0,-18),ha="center",fontsize=8,color="#555")
ax1.set_xscale("log");ax1.invert_xaxis();ax1.set_xticks(thresholds)
ax1.set_xticklabels(["5e-2","1e-2","5e-3","1e-3","5e-4","1e-4"],fontsize=9)
ax1.set_xlabel("cfDNA p-value threshold");ax1.set_ylabel("Concordance (%)")
ax1.set_title("A. Concordance vs cfDNA stringency");ax1.set_ylim(40,80);ax1.legend(loc="upper left")

sc=sig.copy();sc["abs_scot_beta"]=sc["cpg_beta"].abs()
sc=sc.sort_values("abs_scot_beta").reset_index(drop=True)
np.random.seed(42);jit=np.random.uniform(-0.08,0.08,len(sc))
ax2.scatter(sc["abs_scot_beta"],sc["concordant"].astype(int)+jit,s=3,alpha=0.08,color=CB_GREY,zorder=2,rasterized=True)
w=500;run=uniform_filter1d(sc["concordant"].values.astype(float),size=w,mode="nearest")*100
ax2.plot(sc["abs_scot_beta"].values,run,color=CB_PURPLE,linewidth=2.5,zorder=4,label=f"Running avg (n={w})")
pr=run/100;ci=1.96*np.sqrt(pr*(1-pr)/w)*100
ax2.fill_between(sc["abs_scot_beta"].values,run-ci,run+ci,alpha=0.2,color=CB_PURPLE,zorder=3)
ax2.axhline(50,color=CB_GREY,linestyle="--",linewidth=1.2)
ax2.set_xlabel("Blood EWAS effect size (|B|)");ax2.set_ylabel("Concordance (%)")
ax2.set_title("B. Concordance by blood effect size")
ax2.set_ylim(40,70);ax2.set_xlim(0,sc["abs_scot_beta"].quantile(0.99));ax2.legend(loc="upper left",fontsize=9)
plt.tight_layout(w_pad=3);plt.savefig("step2_fig2_stringency_effectsize.png");plt.close()

# FIG 3: chromosome
chroms=[f"chr{i}" for i in range(1,23)]; cd=[]
for c in chroms:
    sub=sig[sig["chrom"]==c]
    if len(sub)==0:continue
    k=sub["concordant"].sum();n=len(sub)
    cd.append({"chr":c.replace("chr",""),"n":n,"pct":100*k/n,"mean_beta":sub["beta"].mean(),"pct_hypo":100*(sub["direction"]=="hypo").mean()})
cd=pd.DataFrame(cd)
fig,axes=plt.subplots(3,1,figsize=(12,10),sharex=True,gridspec_kw={"hspace":0.12})
axes[0].bar(cd["chr"],cd["pct"],color=[CB_TEAL if p>50 else CB_ORANGE for p in cd["pct"]],width=0.7,edgecolor="white",zorder=3)
axes[0].axhline(50,color=CB_GREY,linestyle="--",linewidth=1.2);axes[0].set_ylabel("Concordance (%)");axes[0].set_title("A. Concordance by chromosome (autosomes, p<0.05)");axes[0].set_ylim(30,70)
axes[1].bar(cd["chr"],cd["mean_beta"],color=[CB_BLUE if b<0 else CB_ORANGE for b in cd["mean_beta"]],width=0.7,edgecolor="white",zorder=3)
axes[1].axhline(0,color=CB_BLACK,linewidth=0.8);axes[1].set_ylabel("Mean cfDNA B\n(% meth/year)");axes[1].set_title("B. Mean effect direction")
axes[2].bar(cd["chr"],cd["pct_hypo"],color=CB_BLUE,width=0.7,edgecolor="white",label="Hypo")
axes[2].bar(cd["chr"],100-cd["pct_hypo"],bottom=cd["pct_hypo"],color=CB_ORANGE,width=0.7,edgecolor="white",label="Hyper")
axes[2].axhline(50,color="white",linestyle="--",linewidth=1,alpha=0.7)
axes[2].set_ylabel("% of significant CpGs");axes[2].set_title("C. Hypo/hyper split");axes[2].set_ylim(0,100);axes[2].legend(loc="lower right");axes[2].set_xlabel("Chromosome")
plt.savefig("step2_fig3_chromosome_concordance.png");plt.close()

# Concordance table
tbl=[]
for t in [0.05,0.01,0.001]:
    sub=df[df["pval"]<t]
    if len(sub)==0:continue
    k=sub["concordant"].sum();n=len(sub);bp=binomtest(k,n,0.5,alternative="greater").pvalue
    s1=sub[sub["direction"]=="hypo"];s2=sub[sub["direction"]=="hyper"]
    hk=s1["concordant"].sum();yk=s2["concordant"].sum()
    bph=binomtest(hk,len(s1),0.5,alternative="greater").pvalue if len(s1)>0 else np.nan
    bpy=binomtest(yk,len(s2),0.5,alternative="greater").pvalue if len(s2)>0 else np.nan
    tbl.append({"p_threshold":t,"N":n,"N_concordant":k,"pct_concordant":round(100*k/n,1),
        "binomial_p":f"{bp:.2e}","N_hypo":len(s1),"pct_hypo_conc":round(100*hk/len(s1),1) if len(s1)>0 else "-",
        "binomial_p_hypo":f"{bph:.2e}" if not np.isnan(bph) else "-",
        "N_hyper":len(s2),"pct_hyper_conc":round(100*yk/len(s2),1) if len(s2)>0 else "-",
        "binomial_p_hyper":f"{bpy:.2e}" if not np.isnan(bpy) else "-"})
pd.DataFrame(tbl).to_csv("step2_concordance_table.tsv",sep="\t",index=False)

sig.to_csv("step2_significant_overlap.tsv",sep="\t",index=False)
print(f"Done. Figures: 3, Tables: 2, Step3 input: step2_significant_overlap.tsv ({len(sig):,} rows)")
