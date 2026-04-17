#!/usr/bin/env python
"""
STEP1_QC_reverse.py 

INPUT:
    chinese_cpg_forward_strand.tsv          (from STEP0.py)

OUTPUT:
    QC_chinese_age_association.tsv          (used by QC_STEP2.py)

Pipeline:
  1. Three QC filters: non-zero, SD >= 5, bidirectionality.
  2. Vectorised OLS of methylation (%) on age (midpoints: 26.5 / 60.5 / 81.0).
  3. Two-tailed t-test (df = 33), writes beta, p, direction, group means.
"""
import pandas as pd
import numpy as np
from scipy.stats import t as t_dist

INPUT_FILE  = "chinese_cpg_forward_strand.tsv"
OUTPUT_FILE = "QC_chinese_age_association.tsv"
CHUNK_SIZE  = 500_000  

MIN_SD       = 5.0
MIN_ACTIVE   = 3
MIN_INACTIVE = 3

AGE_OLD   = 81.0
AGE_MID   = 60.5
AGE_YOUNG = 26.5

header      = pd.read_csv(INPUT_FILE, sep="\t", nrows=0)
meta_cols   = ["chrom", "position"]
sample_cols = [c for c in header.columns if c not in meta_cols]

def get_age(col):
    c = col.lower()
    if "old"   in c: return AGE_OLD
    if "mid"   in c: return AGE_MID
    if "young" in c: return AGE_YOUNG
    raise ValueError(f"Cannot assign age: {col}")

ages = np.array([get_age(c) for c in sample_cols], dtype=float)
print(f"Young:{(ages==AGE_YOUNG).sum()}  Mid:{(ages==AGE_MID).sum()}  Old:{(ages==AGE_OLD).sum()}")

# ── Precompute regression constants ONCE (x never changes between CpGs) ──────
x_dev = ages - ages.mean()
ss_xx = (x_dev ** 2).sum()
n     = len(ages)

def regress_chunk(meth):
    """
    Vectorized OLS for all rows at once.
    meth: [N_cpgs x N_samples], no NaNs allowed (filter before calling)
    Returns: beta [N], pval [N]
    """
    y_mean    = meth.mean(axis=1)                              # [N]
    beta      = (meth - y_mean[:, None]) @ x_dev / ss_xx       # [N]
    intercept = y_mean - beta * ages.mean()                    # [N]
    residuals = meth - (intercept[:, None] + beta[:, None] * ages)  # [N x S]
    ss_res    = (residuals ** 2).sum(axis=1)                   # [N]
    se_beta   = np.sqrt(ss_res / (n - 2) / ss_xx)              # [N]
    t_stat    = beta / np.where(se_beta == 0, np.nan, se_beta) # [N]
    pval      = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)        # [N]
    return beta, pval

tot_in = tot_f1 = tot_f2 = tot_f3 = 0
tot_sig = tot_hyper = tot_hypo = 0
first_write = True

for chunk_i, df in enumerate(pd.read_csv(INPUT_FILE, sep="\t",
                                          chunksize=CHUNK_SIZE)):

    meth = df[sample_cols].values.astype(float)
    tot_in += len(df)

    # F1: not all-zero
    keep = np.nanmax(meth, axis=1) > 0
    tot_f1 += keep.sum()

    # F2: SD >= 5
    keep[keep] &= np.nanstd(meth[keep], axis=1) >= MIN_SD
    tot_f2 += keep.sum()

    # F3: bidirectional
    keep[keep] &= ((meth[keep] > 0).sum(axis=1) >= MIN_ACTIVE) & \
                  ((meth[keep] < 100).sum(axis=1) >= MIN_INACTIVE)
    tot_f3 += keep.sum()

    meth_f = meth[keep]
    df_f   = df[keep].reset_index(drop=True)

    if len(df_f) == 0:
        continue

    # Drop rows with any NaN before vectorized regression
    # (NaN rows = samples with zero coverage at that position)
    has_nan = np.isnan(meth_f).any(axis=1)
    meth_clean  = meth_f[~has_nan]
    df_clean    = df_f[~has_nan].reset_index(drop=True)

    if len(df_clean) == 0:
        continue

    beta, pval = regress_chunk(meth_clean)

    out = pd.DataFrame({
        "chrom":      df_clean["chrom"].values,
        "position":   df_clean["position"].values,
        "beta":       beta,
        "pval":       pval,
        "mean_young": meth_clean[:, ages == AGE_YOUNG].mean(axis=1),
        "mean_mid":   meth_clean[:, ages == AGE_MID].mean(axis=1),
        "mean_old":   meth_clean[:, ages == AGE_OLD].mean(axis=1),
        "direction":  np.where(beta > 0, "hyper", "hypo")
    })

    sig = out[out["pval"] < 0.05]
    tot_sig   += len(sig)
    tot_hyper += (sig["direction"] == "hyper").sum()
    tot_hypo  += (sig["direction"] == "hypo").sum()

    out.to_csv(OUTPUT_FILE, sep="\t", index=False,
               mode="w" if first_write else "a", header=first_write)
    first_write = False

    if chunk_i % 10 == 0:
        print(f"Chunk {chunk_i:>4} | in:{len(df):>8,} | "
              f"F1:{(np.nanmax(meth,axis=1)>0).sum():>8,} | "
              f"kept:{len(df_clean):>8,}", flush=True)

print(f"\n── FILTER SUMMARY ───────────────────────────")
print(f"  Input:                   {tot_in:>12,}")
print(f"  After F1 (not all-zero): {tot_f1:>12,}  (-{tot_in-tot_f1:,})")
print(f"  After F2 (SD >= 5):      {tot_f2:>12,}  (-{tot_f1-tot_f2:,})")
print(f"  After F3 (bidirectional):{tot_f3:>12,}  (-{tot_f2-tot_f3:,})")
print(f"\n── RESULTS ──────────────────────────────────")
print(f"  Significant p<0.05:      {tot_sig:>12,}")
print(f"  Hypermethylating:        {tot_hyper:>12,}")
print(f"  Hypomethylating:         {tot_hypo:>12,}")
print(f"\nSaved: {OUTPUT_FILE}")