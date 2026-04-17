#!/usr/bin/env Rscript
# =============================================================
# STEP4A_missmethyl.R — Probe-bias-corrected GO enrichment
#
# Runs gometh() on TWO input sets:
#   (A) CGI CpGs: 143 concordant vs 96 discordant (bg = 239)
#   (B) ALL CpGs: ~2344 concordant vs ~1947 discordant (bg = 4291)
#
# For each set, tests:
#   - GO (CGI/all-specific background)
#   - GO (full EPIC background)
#   - KEGG
#
# INPUT:
#   cgi_concordant_probes.txt, cgi_discordant_probes.txt, cgi_all_probes.txt
#   all_concordant_probes.txt, all_discordant_probes.txt, all_significant_probes.txt
#   - All from STEP4_prep.py
#
# OUTPUT:
#   missmethyl_{cgi,all}_{concordant,discordant}_GO.tsv
#   missmethyl_{cgi,all}_{concordant,discordant}_GO_fullbg.tsv
#   missmethyl_{cgi,all}_{concordant,discordant}_KEGG.tsv
#   missmethyl_bias_plot.pdf
# =============================================================
library(missMethyl)

run_enrichment <- function(sig_file, bg_file, label, set_name) {
    # Load probes
    sig <- readLines(sig_file)
    bg  <- readLines(bg_file)
    cat(sprintf("\n══ %s — %s ══\n", set_name, label))
    cat(sprintf("  Sig: %d probes | Background: %d probes\n", length(sig), length(bg)))

    prefix <- paste0(set_name, "_", label)

    # --- GO with specific background ---
    cat("  Running GO (specific background)...\n")
    go <- gometh(sig.cpg = sig, all.cpg = bg,
                 collection = "GO", array.type = "EPIC",
                 plot.bias = FALSE, prior.prob = TRUE)
    go <- go[order(go$P.DE), ]
    go_top <- head(go, 30)
    outf <- paste0("missmethyl_", prefix, "_GO.tsv")
    write.table(go_top, outf, sep="\t", quote=FALSE, row.names=TRUE)
    cat(sprintf("  Saved: %s  (top P.DE = %.4f, FDR = %.4f)\n",
                outf, go_top$P.DE[1], go_top$FDR[1]))

    # Print top 5
    cat("  Top 5 GO terms:\n")
    for (i in 1:min(5, nrow(go_top))) {
        cat(sprintf("    %s (P=%.4f, FDR=%.4f)\n",
                    go_top$TERM[i], go_top$P.DE[i], go_top$FDR[i]))
    }

    # --- GO with full EPIC background ---
    cat("  Running GO (full EPIC background)...\n")
    go_full <- gometh(sig.cpg = sig,
                      collection = "GO", array.type = "EPIC",
                      plot.bias = FALSE, prior.prob = TRUE)
    go_full <- go_full[order(go_full$P.DE), ]
    go_full_top <- head(go_full, 30)
    outf <- paste0("missmethyl_", prefix, "_GO_fullbg.tsv")
    write.table(go_full_top, outf, sep="\t", quote=FALSE, row.names=TRUE)
    cat(sprintf("  Saved: %s  (top P.DE = %.4f, FDR = %.4f)\n",
                outf, go_full_top$P.DE[1], go_full_top$FDR[1]))

    # Check for granulocyte terms
    gran_terms <- go_full[grepl("granulocyte|myeloid|neutrophil|leukocyte|hematopoie|haematopoie",
                                 go_full$TERM, ignore.case=TRUE), ]
    if (nrow(gran_terms) > 0) {
        cat("  Immune/haematopoietic terms found:\n")
        for (i in 1:min(5, nrow(gran_terms))) {
            cat(sprintf("    %s (P=%.4f, FDR=%.4f, rank=%d/%d)\n",
                        gran_terms$TERM[i], gran_terms$P.DE[i],
                        gran_terms$FDR[i],
                        which(rownames(go_full) == rownames(gran_terms)[i]),
                        nrow(go_full)))
        }
    }

    # --- KEGG ---
    cat("  Running KEGG (specific background)...\n")
    kegg <- gometh(sig.cpg = sig, all.cpg = bg,
                   collection = "KEGG", array.type = "EPIC",
                   plot.bias = FALSE, prior.prob = TRUE)
    kegg <- kegg[order(kegg$P.DE), ]
    kegg_top <- head(kegg, 20)
    outf <- paste0("missmethyl_", prefix, "_KEGG.tsv")
    write.table(kegg_top, outf, sep="\t", quote=FALSE, row.names=TRUE)
    cat(sprintf("  Saved: %s  (top P.DE = %.4f, FDR = %.4f)\n",
                outf, kegg_top$P.DE[1], kegg_top$FDR[1]))
}

# ── Save bias plot for CGI concordant ─────────────────────────
cat("Generating bias plot...\n")
pdf("missmethyl_bias_plot.pdf", width=6, height=5)
sig_cgi <- readLines("cgi_concordant_probes.txt")
bg_cgi  <- readLines("cgi_all_probes.txt")
gometh(sig.cpg = sig_cgi, all.cpg = bg_cgi,
       collection = "GO", array.type = "EPIC",
       plot.bias = TRUE, prior.prob = TRUE)
dev.off()
cat("Saved: missmethyl_bias_plot.pdf\n")

# ── CGI SET ───────────────────────────────────────────────────
run_enrichment("cgi_concordant_probes.txt", "cgi_all_probes.txt",
               "concordant", "cgi")
run_enrichment("cgi_discordant_probes.txt", "cgi_all_probes.txt",
               "discordant", "cgi")

# ── ALL SIGNIFICANT SET ───────────────────────────────────────
# Only run if probe files exist (STEP4_prep.py may skip these)
if (file.exists("all_concordant_probes.txt")) {
    run_enrichment("all_concordant_probes.txt", "all_significant_probes.txt",
                   "concordant", "all")
    run_enrichment("all_discordant_probes.txt", "all_significant_probes.txt",
                   "discordant", "all")
} else {
    cat("\nSKIPPED: ALL set probe files not found.\n")
}

cat("\n══════════════════════════════════════\n")
cat("Done. Check output TSV files.\n")
cat("══════════════════════════════════════\n")
