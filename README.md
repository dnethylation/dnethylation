# cfDNA × blood EWAS concordance analysis

BSc Medical Genetics dissertation codebase. Queen Mary University of London, 2025/26.

Tests directional concordance between age-associated CpG methylation in cfDNA
(Li et al. 2024, GSE259312, n = 35) and the Generation Scotland blood EWAS
(Bernabeu et al. 2023, n = 18,413), characterises the genomic / chromatin
context of concordant and discordant positions, and estimates age-associated
changes in cfDNA tissue composition (Loyfer et al. 2023 atlas).

## Pipeline

Run scripts in order. Each documents its inputs, outputs, and downstream
consumers at the top of the file.

1. `STEP0.py`                      — separate Bismark output by strand
2. `STEP1.py`                      — per-CpG OLS regression on forward strand
3. `QC_STEP1.py`                   — same as STEP1 on reverse strand (validation only)
4. `STEP2.py`                      — overlap with blood EWAS + concordance
5. `QC_STEP2.py`                   — reverse-strand concordance (validation only)
6. `STEP2_results.py`              — concordance figures and tables
7. `STEP3A_results.py`             — CpG island + gene feature analysis
8. `STEP3B_results.py`             — gene-level + epigenetic clock overlap
9. `STEP3C_deconvolution.py`       — tissue-of-origin NNLS deconvolution
10. `STEP3C_results.py`            — deconvolution figures and Spearman
11. `STEP4_prep.py`                — prepare probe lists + BED files
12. `STEP4A_missmethyl.R`          — probe-bias-corrected GO / KEGG enrichment
13. `STEP4B_chromhmm.py`           — Roadmap ChromHMM chromatin state annotation

| Step / Script        | Output files                               | Used by / Appears in                  |
| -------------------- | ------------------------------------------ | ------------------------------------- |
| STEP0                | chinese_cpg_forward_strand.tsv             | STEP1, STEP3C                         |
|                      | chinese_cpg_reverse_strand.tsv             | QC_STEP1                      |
| STEP1                | chinese_age_association.tsv                | STEP2                                 |
| QC_STEP1             | QC_chinese_age_association.tsv             | QC_STEP2                      |
| STEP2                | overlap_chinese+scottish.tsv               | STEP2_results                         |
|                      | concordance_summary.tsv                    | Superseded (kept as log)              |
| QC_STEP2             | QC_overlap_chinese+scottish.tsv            | Internal log                          |
|                      | QC_concordance_summary.tsv                 | Dissertation §3.1 / §4.1              |
| STEP2_results        | step2_significant_overlap.tsv              | STEP3A, STEP4_prep                    |
|                      | step2_concordance_table.tsv                | Dissertation (Supplementary Table)    |
|                      | step2_fig1/fig2/fig3                       | Figures 4, 5, 6                       |
| STEP3A_results       | cgi_annotated_full.tsv                     | STEP3B, STEP4_prep                    |
|                      | step3a_enrichment_table.tsv                | Dissertation (Supplementary Table)    |
|                      | step3a_fig1/fig2/fig3                      | Figures 7, 8, 9                       |
| STEP3B_results       | step3b_gene_table.tsv                      | Dissertation (Supplementary Table)    |
|                      | step3b_clock_overlap.tsv                   | Dissertation §3.4                     |
|                      | step3b_literature_table.tsv                | Curated reference material            |
|                      | step3b_fig1/fig2                           | Figure 11 + inline reference          |
|                      | great_*.bed                                | GREAT web submission (Methods §2.6.5) |
| STEP3C_deconvolution | deconvolution_results.tsv                  | STEP3C                                |
|                      | deconvolution_age_summary.tsv              | Internal log                          |
| STEP3C_results       | step3c_spearman_table.tsv                  | Dissertation (Supplementary Table)    |
|                      | step3c_fig1/fig2/fig3                      | Figures 12, 13, 14                    |
| STEP4_prep           | cgi_*_probes.txt, all_*_probes.txt         | STEP4A                                |
|                      | cgi_*.bed, all_*.bed                       | STEP4B                                |
| STEP4A_missmethyl    | missmethyl_*_GO.tsv, missmethyl_*_KEGG.tsv | Dissertation §3.5.2                   |
|                      | missmethyl_bias_plot.pdf                   | QC visual                             |
| STEP4B_chromhmm      | chromhmm_*_summary.tsv                     | Dissertation §3.3.4                   |
|                      | chromhmm_*_state_counts.tsv                | Supplementary               |
|                      | chromhmm_*_fig_main/blood.png              | Figure 10, Figures S1/S2/S3           |
|                      | chromhmm_pcdhga_states.tsv                 | Dissertation §3.4 / §4.3              |
## External data sources

The following files are required but NOT included in this repository:

| File | Source |
|------|--------|
| `GSE259312_Bismark_Methylated_rate.tsv` | GEO accession GSE259312 (Li et al. 2024) |
| `ewas_age_linear.tsv` | Bernabeu et al. 2023 (Generation Scotland EWAS summary stats) |
| `atlas_markers.tsv` | Loyfer et al. 2023 Supplementary Table 4a |
| `cpgIslandExt.txt` | UCSC Genome Browser (hg19) |
| `infinium-methylationepic-v-1-0-b5-manifest-file.csv` | Illumina  |
| `gb-2013-14-10-r115-S3.csv` | Horvath 2013 supplement |
| `NIHMS418935-supplement-03.csv` | Hannum et al. 2013 supplement |

Place these in the project root before running.

Roadmap Epigenomics ChromHMM segmentations are auto-downloaded by `STEP4B_chromhmm.py`.


## Notes on methodology

Full methods and results are described in the accompanying BSc dissertation.
Key design choices:

- Age is encoded as group-midpoint (26.5 / 60.5 / 81.0 years) because
  individual-level ages were not available in GSE259312.
- The EPIC epigenome-wide threshold p < 3.6e-8 is used for the blood EWAS
  filter (Saffari et al. 2018).
- Concordance is tested with a one-sided binomial test against 50%;
  limitations of this assumption are discussed in Section 4.1 of the
  dissertation.

