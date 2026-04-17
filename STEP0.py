#!/usr/bin/env python
"""
STEP0_strand_separation.py — split Bismark CpG output into forward and reverse strands.

INPUT:
    GSE259312_Bismark_Methylated_rate.tsv   ( GEO accession GSE259312)

OUTPUT:
    chinese_cpg_forward_strand.tsv          (used by STEP1.py and STEP3C_deconvolution.py)
    chinese_cpg_reverse_strand.tsv          (used by STEP1_QC_reverse.py only)

The forward strand is used for all primary analyses; the reverse strand is
processed independently through STEP1_QC_reverse and STEP2_QC_reverse as a
validation check .
"""

import pandas as pd
import numpy as np

INPUT = "GSE259312_Bismark_Methylated_rate.tsv"

df = pd.read_csv(INPUT, sep="\t")
print(f"Total rows: {len(df):,}")

df = df.sort_values(["chrom", "position"]).reset_index(drop=True)

pos   = df["position"].values
chrom = df["chrom"].values

# Mark reverse strand rows (G-strand)
is_g = np.zeros(len(df), dtype=bool)
is_g[1:] = (pos[1:] == pos[:-1] + 1) & (chrom[1:] == chrom[:-1])

# Forward strand (C-strand) — keep where NOT G
fwd = df[~is_g].reset_index(drop=True)

# Reverse strand (G-strand) — keep where IS G
rev = df[is_g].copy().reset_index(drop=True)
# Shift reverse strand position back by 1 so both strands share same coordinate key
rev["position"] = rev["position"] - 1

print(f"Forward strand rows: {len(fwd):,}")
print(f"Reverse strand rows: {len(rev):,}")

fwd.to_csv("chinese_cpg_forward_strand.tsv", sep="\t", index=False)
rev.to_csv("chinese_cpg_reverse_strand.tsv", sep="\t", index=False)
print("Both saved.")
