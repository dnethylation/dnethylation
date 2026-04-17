#!/usr/bin/env python
"""
QC_STEP2.py — overlap cfDNA age QC with the Generation Scotland blood EWAS.

INPUT:
    ewas_age_linear.tsv                     (Bernabeu et al., n = 18,413)
    QC_chinese_age_association.tsv          (from STEP1.py)

OUTPUT:
    QC_ overlap_chinese+scottish.tsv           (not used, leftover)
    QC_concordance_summary.tsv                 (minimal log of the QC overlap state)

                                             
Applies the EPIC epigenome-wide threshold p < 3.6e-8 to blood and performs an
inner join on chr/hg19-position. Direction of effect recorded for each dataset
and a boolean concordance column added.
"""

import pandas as pd
import numpy as np
from scipy.stats import binomtest

SCOT_FILE  = "ewas_age_linear.tsv"
CHINA_FILE = "QC_chinese_age_association.tsv"
OUT_FILE   = "QC_overlap_chinese+scottish.tsv"
OUT_SUMM   = "QC_concordance_summary.tsv"
PTHRESH    = 3.6e-8

# ── Scottish EWAS ─────────────────────────────────────────────────────────────
scot = pd.read_csv(SCOT_FILE, sep="\t")
scot = scot[scot["cpg_p"] < PTHRESH].reset_index(drop=True)
scot["chr"] = "chr" + scot["chr"].astype(str)
scot["direction_scot"] = np.where(scot["cpg_beta"] > 0, "hyper", "hypo")

# ── Chinese age association ───────────────────────────────────────────────────
chin = pd.read_csv(CHINA_FILE, sep="\t")

# ── Inner join on chr + position ──────────────────────────────────────────────
merged = chin.merge(
    scot[["cpg", "chr", "pos", "cpg_beta", "cpg_p", "direction_scot"]],
    left_on=["chrom", "position"],
    right_on=["chr", "pos"],
    how="inner"
).drop(columns=["chr", "pos"])

merged["concordant"] = merged["direction"] == merged["direction_scot"]

merged.to_csv(OUT_FILE, sep="\t", index=False)

# ── Concordance summary ───────────────────────────────────────────────────────
rows = []
for t in [0.05, 0.01, 0.001]:
    sub = merged[merged["pval"] < t]
    if len(sub) == 0:
        continue
    k = sub["concordant"].sum()
    n = len(sub)
    bp = binomtest(k, n, 0.5, alternative="greater").pvalue
    hypo = sub[sub["direction"] == "hypo"]
    hyper = sub[sub["direction"] == "hyper"]
    rows.append({
        "p_thresh": t, "N": n, "concordant": k,
        "pct_concordant": round(100*k/n, 2),
        "binom_p": bp,
        "pct_hypo_conc":  round(100*hypo["concordant"].mean(), 2) if len(hypo) > 0 else np.nan,
        "pct_hyper_conc": round(100*hyper["concordant"].mean(), 2) if len(hyper) > 0 else np.nan,
    })

summary = pd.DataFrame(rows)
summary.to_csv(OUT_SUMM, sep="\t", index=False)

print(f"Overlap: {len(merged):,} CpGs")
print(f"Saved: {OUT_FILE}")
print(f"Saved: {OUT_SUMM}")
print(summary[summary["p_thresh"]==0.05].to_string(index=False))
