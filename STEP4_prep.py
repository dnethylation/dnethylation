#!/usr/bin/env python3
"""
STEP4_prep.py — prepare probe lists and BED files for missMethyl and ChromHMM.


Generates probe lists and BED files for:
  (A) CGI CpGs only (n=239: 143 concordant, 96 discordant)
  (B) ALL significant overlap CpGs (n=4,291: concordant vs discordant)

These files feed into:
  - missMethyl (R): probe-bias-corrected GO enrichment
  - ChromHMM annotation (Python): chromatin state context

INPUT FILES (must be in working directory):
  - cgi_annotated_full.tsv          (239 rows, from STEP3A)
  - step2_significant_overlap.tsv   (4,291 rows, from STEP2_results)

OUTPUT FILES:
  CGI set:
    cgi_concordant_probes.txt       probe IDs for missMethyl
    cgi_discordant_probes.txt       probe IDs for missMethyl
    cgi_all_probes.txt              background for missMethyl
    cgi_concordant.bed              for ChromHMM
    cgi_discordant.bed              for ChromHMM

  ALL significant set:
    all_concordant_probes.txt       probe IDs for missMethyl
    all_discordant_probes.txt       probe IDs for missMethyl
    all_significant_probes.txt      background for missMethyl
    all_concordant.bed              for ChromHMM
    all_discordant.bed              for ChromHMM
"""

import pandas as pd
import sys
import os

# ══════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════

print("═" * 60)
print("STEP4_prep.py — Generating input files")
print("═" * 60)

# --- CGI file ---
cgi_file = "cgi_annotated_full.tsv"
if not os.path.exists(cgi_file):
    sys.exit(f"ERROR: {cgi_file} not found. Run STEP3A_results.py first.")

cgi = pd.read_csv(cgi_file, sep="\t")
print(f"\nLoaded {cgi_file}: {len(cgi)} rows")
print(f"  Concordant: {cgi['concordant'].sum()}")
print(f"  Discordant: {(~cgi['concordant']).sum()}")

# --- Full significant overlap file ---
sig_file = "step2_significant_overlap.tsv"
if not os.path.exists(sig_file):
    sys.exit(f"ERROR: {sig_file} not found. Run STEP2_results.py first.")

sig = pd.read_csv(sig_file, sep="\t")
print(f"\nLoaded {sig_file}: {len(sig)} rows")

# Check required columns exist
# The file should have: chrom, position, cpg (probe ID), concordant
required_cols = ["chrom", "position", "concordant"]
for col in required_cols:
    if col not in sig.columns:
        sys.exit(f"ERROR: column '{col}' not found in {sig_file}")

# Check for probe ID column (might be 'cpg' or 'probe_id')
probe_col = None
for candidate in ["cpg", "probe_id", "IlmnID", "Name"]:
    if candidate in sig.columns:
        probe_col = candidate
        break

if probe_col is None:
    print("\n  WARNING: No probe ID column found in significant overlap file.")
    print("  Columns available:", list(sig.columns))
    print("  missMethyl probe lists will NOT be generated for the ALL set.")
    print("  ChromHMM BED files will still be generated.\n")
    has_probes_all = False
else:
    print(f"  Probe ID column: '{probe_col}'")
    has_probes_all = True

conc_all = sig[sig["concordant"] == True]
disc_all = sig[sig["concordant"] == False]
print(f"  Concordant: {len(conc_all)}")
print(f"  Discordant: {len(disc_all)}")


# ══════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════

def save_probes(df, probe_col, outfile):
    """Save probe IDs to text file (one per line)."""
    df[probe_col].dropna().to_csv(outfile, index=False, header=False)
    print(f"  Saved: {outfile} ({len(df)} probes)")

def save_bed(df, outfile):
    """Save BED file (chrom, start, end) sorted."""
    bed = df[["chrom", "position"]].copy()
    bed["end"] = bed["position"] + 1
    bed = bed.rename(columns={"position": "start"})
    bed = bed[["chrom", "start", "end"]].sort_values(["chrom", "start"])
    bed.to_csv(outfile, sep="\t", index=False, header=False)
    print(f"  Saved: {outfile} ({len(bed)} regions)")


# ══════════════════════════════════════════════════════════════
# 3. CGI SET (n=239)
# ══════════════════════════════════════════════════════════════

print("\n── CGI set (n=239) ──────────────────────────────────")

cgi_conc = cgi[cgi["concordant"] == True]
cgi_disc = cgi[cgi["concordant"] == False]

# Probe lists for missMethyl
save_probes(cgi_conc, "cpg", "cgi_concordant_probes.txt")
save_probes(cgi_disc, "cpg", "cgi_discordant_probes.txt")
save_probes(cgi, "cpg", "cgi_all_probes.txt")

# BED files for ChromHMM
save_bed(cgi_conc, "cgi_concordant.bed")
save_bed(cgi_disc, "cgi_discordant.bed")


# ══════════════════════════════════════════════════════════════
# 4. ALL SIGNIFICANT SET (n=4,291)
# ══════════════════════════════════════════════════════════════

print("\n── ALL significant set (n=4,291) ────────────────────")

# Probe lists for missMethyl (if probe IDs available)
if has_probes_all:
    save_probes(conc_all, probe_col, "all_concordant_probes.txt")
    save_probes(disc_all, probe_col, "all_discordant_probes.txt")
    save_probes(sig, probe_col, "all_significant_probes.txt")
else:
    print("  SKIPPED: probe lists (no probe ID column)")

# BED files for ChromHMM
save_bed(conc_all, "all_concordant.bed")
save_bed(disc_all, "all_discordant.bed")


# ══════════════════════════════════════════════════════════════
# 5. SUMMARY
# ══════════════════════════════════════════════════════════════

print("\n" + "═" * 60)
print("OUTPUT FILES")
print("═" * 60)
print("\nCGI set:")
print("  cgi_concordant_probes.txt    →  missMethyl sig list")
print("  cgi_discordant_probes.txt    →  missMethyl sig list")
print("  cgi_all_probes.txt           →  missMethyl background")
print("  cgi_concordant.bed           →  ChromHMM input")
print("  cgi_discordant.bed           →  ChromHMM input")
print("\nALL significant set:")
if has_probes_all:
    print("  all_concordant_probes.txt    →  missMethyl sig list")
    print("  all_discordant_probes.txt    →  missMethyl sig list")
    print("  all_significant_probes.txt   →  missMethyl background")
print("  all_concordant.bed           →  ChromHMM input")
print("  all_discordant.bed           →  ChromHMM input")
print("\nNext steps:")
print("  1. Run chromhmm_annotation.py  (uses .bed files)")
print("  2. Run missmethyl_enrichment.R (uses .txt probe lists)")
print("═" * 60)
