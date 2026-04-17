#!/usr/bin/env python
"""
STEP3C_deconvolution.py — tissue-of-origin NNLS deconvolution .

INPUT:
    atlas_markers.tsv                       (Loyfer et al. Supplementary Table 4a)
    chinese_cpg_forward_strand.tsv          (from STEP0.py)

OUTPUT:
    deconvolution_results.tsv               (used by STEP3D_deconvolution_figs.py)
    deconvolution_age_summary.tsv           (internal log; not cited in dissertation)
"""

import pandas as pd
import numpy as np
from scipy.optimize import nnls

# ── INPUTS ────────────────────────────────────────────────────────────────────
MARKERS_FILE = "atlas_markers.tsv"        
CHINESE_FILE = "chinese_cpg_forward_strand.tsv"  # your 20.5M row methylation file
OUTPUT_FILE  = "deconvolution_results.tsv"
CHUNKSIZE    = 500000

# Age group labels (based on column name containing young/mid/old)
AGE_MAP = {"young": 26.5, "mid": 60.5, "old": 81.0}

# ════════════════════════════════════════════════════════════════════════════
# 1. LOAD MARKERS AND BUILD REFERENCE ATLAS MATRIX
# ════════════════════════════════════════════════════════════════════════════
print("Loading atlas markers...")
markers = pd.read_csv(MARKERS_FILE, sep="\t")

# Standardise column names in case of slight differences
markers.columns = markers.columns.str.strip()

# Rename columns if needed
if "Target meth." in markers.columns:
    markers = markers.rename(columns={"Target meth.": "target_meth",
                                       "Background meth.": "bg_meth"})

# Ensure chr format matches your Chinese data (needs "chr" prefix)
markers["chr"] = markers["chr"].apply(lambda x: x if str(x).startswith("chr") else "chr" + str(x))

# Get sorted list of all cell types
cell_types = sorted(markers["Type"].unique())
print(f"  Cell types in atlas: {len(cell_types)}")
print(f"  Marker regions total: {len(markers)}")

# Build reference matrix: rows = markers, cols = cell types
# For each marker: target cell type gets target_meth, all others get bg_meth
ref_matrix = pd.DataFrame(index=range(len(markers)), columns=cell_types, dtype=float)

for i, row in markers.iterrows():
    for ct in cell_types:
        if ct == row["Type"]:
            ref_matrix.loc[i, ct] = row["target_meth"]
        else:
            ref_matrix.loc[i, ct] = row["bg_meth"]

ref_array = ref_matrix.values  # shape: (n_markers, n_cell_types)
print(f"  Reference matrix shape: {ref_array.shape}")

# ════════════════════════════════════════════════════════════════════════════
# 2. LOAD CHINESE DATA HEADER — GET SAMPLE NAMES AND AGE GROUPS
# ════════════════════════════════════════════════════════════════════════════
print("\nReading Chinese data header...")
header = pd.read_csv(CHINESE_FILE, sep="\t", nrows=0)
meta_cols = ["chrom", "position"]
sample_cols = [c for c in header.columns if c not in meta_cols]

def get_age(col):
    c = col.lower()
    for key, val in AGE_MAP.items():
        if key in c:
            return val
    return None

sample_ages = {s: get_age(s) for s in sample_cols}
print(f"  Samples found: {len(sample_cols)}")
for age_group, age_val in AGE_MAP.items():
    n = sum(1 for v in sample_ages.values() if v == age_val)
    print(f"    {age_group}: {n} samples")

# ════════════════════════════════════════════════════════════════════════════
# 3. EXTRACT METHYLATION AT EACH MARKER REGION PER SAMPLE
# ════════════════════════════════════════════════════════════════════════════
# For each marker region, find all CpGs within its coordinates
# and compute the average methylation per sample

print("\nExtracting methylation at marker regions...")
print("  (This scans 20M rows in chunks — takes ~5-10 minutes)")

# Initialise accumulator: for each marker, sum of methylation and count of CpGs
n_markers = len(markers)
meth_sum   = np.zeros((n_markers, len(sample_cols)))  # total methylation
meth_count = np.zeros(n_markers)                       # number of CpGs found

# Build a fast lookup: for each chromosome, list of (start, end, marker_index)
from collections import defaultdict
chr_markers = defaultdict(list)
for idx, row in markers.iterrows():
    chr_markers[row["chr"]].append((row["start"], row["end"], idx))

# Sort each chromosome's markers by start position for fast search
for chrom in chr_markers:
    chr_markers[chrom].sort()

chunk_n = 0
for chunk in pd.read_csv(CHINESE_FILE, sep="\t", chunksize=CHUNKSIZE):
    chunk_n += 1
    if chunk_n % 10 == 0:
        print(f"  chunk {chunk_n}...", flush=True)

    for chrom, grp in chunk.groupby("chrom"):
        if chrom not in chr_markers:
            continue
        c_markers = chr_markers[chrom]
        positions  = grp["position"].values
        meth_vals  = grp[sample_cols].values  # shape: (n_positions, n_samples)

        for (start, end, midx) in c_markers:
            # Find positions within [start, end]
            in_region = (positions >= start) & (positions <= end)
            if not in_region.any():
                continue
            region_meth = meth_vals[in_region]        # shape: (n_cpgs, n_samples)
            # Average across CpGs (ignoring NaN)
            meth_sum[midx]   += np.nanmean(region_meth, axis=0)
            meth_count[midx] += 1

print(f"  Done. Markers with ≥1 CpG covered: {(meth_count > 0).sum()} / {n_markers}")

# ════════════════════════════════════════════════════════════════════════════
# 4. BUILD OBSERVED METHYLATION MATRIX
# ════════════════════════════════════════════════════════════════════════════
# Average across chunks where marker was found
obs_matrix = np.full((n_markers, len(sample_cols)), np.nan)
covered = meth_count > 0
obs_matrix[covered] = meth_sum[covered] / meth_count[covered, np.newaxis]

# Convert to 0-1 scale (your data is 0-100)
obs_matrix = obs_matrix / 100.0

# Drop markers with no coverage
covered_idx   = np.where(covered)[0]
obs_covered   = obs_matrix[covered_idx]   # shape: (n_covered_markers, n_samples)
ref_covered   = ref_array[covered_idx]    # shape: (n_covered_markers, n_cell_types)

print(f"\nMarkers covered: {len(covered_idx)}")
print(f"Proceeding with {len(covered_idx)} markers for deconvolution")

# ════════════════════════════════════════════════════════════════════════════
# 5. NNLS DECONVOLUTION PER SAMPLE
# ════════════════════════════════════════════════════════════════════════════
print("\nRunning NNLS deconvolution per sample...")

results = []
for si, sample in enumerate(sample_cols):
    y = obs_covered[:, si]                # observed methylation at covered markers
    valid = ~np.isnan(y)
    if valid.sum() < 10:
        print(f"  WARNING: {sample} has <10 valid markers, skipping")
        continue

    # NNLS: find proportions x such that ref_covered[valid] @ x ≈ y[valid]
    proportions, residual = nnls(ref_covered[valid], y[valid])

    # Normalise to sum to 1
    total = proportions.sum()
    if total > 0:
        proportions = proportions / total

    row = {"sample": sample, "age": sample_ages[sample]}
    for ct, prop in zip(cell_types, proportions):
        row[ct] = round(prop, 4)
    results.append(row)
    print(f"  {sample} (age {sample_ages[sample]}): top tissue = "
          f"{cell_types[np.argmax(proportions)]} ({max(proportions):.1%})")

results_df = pd.DataFrame(results)

# ════════════════════════════════════════════════════════════════════════════
# 6. SUMMARISE BY AGE GROUP
# ════════════════════════════════════════════════════════════════════════════
print("\n══ Mean tissue proportions by age group ══")
summary = results_df.groupby("age")[cell_types].mean().round(3)
print(summary.T.to_string())  # transpose so tissues are rows

# ════════════════════════════════════════════════════════════════════════════
# 7. SAVE
# ════════════════════════════════════════════════════════════════════════════
results_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"\nSaved: {OUTPUT_FILE}")

# Also save the age-group summary
summary.T.to_csv("deconvolution_age_summary.tsv", sep="\t")
print("Saved: deconvolution_age_summary.tsv")
