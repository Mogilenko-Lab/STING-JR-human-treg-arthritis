#!/usr/bin/env Rscript
# 03_pseudobulk_volcano_viz.R: VIZ. One figure — the sorted-Treg synovial-fluid-versus-paired-blood
# donor-pseudobulk contrast, drawn with the canonical volcano grammar.
# =============================================================================
#   treg_volcano   the committed limma-voom Treg contrast on the standard volcano
#
# WHY THIS IS AN R SCRIPT, SPLIT OUT OF 03_pseudobulk_de_viz.py:
# The volcano grammar is the RNAseq-toolkit plotter
# (01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R), which is R. Drawing this
# contrast by hand in matplotlib alongside the stage's other two panels put a second, divergent
# volcano dialect in the compartment: the hand-rolled version plotted -log10(padj) on the y axis,
# which collapses many raw p onto one adjusted value and renders the stair-step/plateau artefact
# the toolkit plotter exists to avoid, and it coloured two categories where the toolkit names
# four. 18_tf_selective_volcano_viz.R already draws THIS SAME contrast through the toolkit, so the
# compartment was showing one DE table under two different volcano conventions. It now has one.
#
# What changes on the canvas, relative to the retired matplotlib panel:
#   - y is raw p (-log10), decision is FDR — resolution on the axis, multiple-testing control on
#     the colour. The dashed horizontal rule sits at the raw p that realises the FDR cut.
#   - four significance categories (NS / Log2FC / p-value / p-value & Log2FC), not two.
#   - labels are the top genes BY SIGNIFICANCE on each side, not the largest |log2FC|.
# The decision rule, its thresholds and the gene set they select are unchanged: the committed
# thresholds.de_fdr / thresholds.de_logfc, applied to the committed DE table.
#
# project_theme() is applied afterwards, since it is a complete theme and has to be added last.
#
# Reads:
#   03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv   the committed limma-voom table
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/03_pseudobulk_volcano_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")  # create_standard_volcano

STAGE  <- "03_pseudobulk"
SCRIPT <- "02_analysis/scripts/03_pseudobulk_volcano_viz.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

CFG <- FIG_CFG
FIG <- CFG$figures
THR <- CFG$thresholds

# The compartment's committed DE decision rule, applied unchanged.
FDR_CUT   <- as.numeric(THR$de_fdr %||% 0.05)
LOGFC_CUT <- as.numeric(THR$de_logfc %||% 1.0)

# `volcano_label_top` is a TOTAL label budget; the plotter's `top_n` is per side, so it halves.
LABEL_TOP <- as.integer(FIG$volcano_label_top %||% 10L)
TOP_N     <- LABEL_TOP %/% 2L
# Rows carried in the neighbour table. The full contrast is the committed DE table next door.
TABLE_N   <- 500L

RES <- CFG$paths$results %||% "03_results/"
DE  <- readr::read_csv(file.path(RES, STAGE, "tables", "de_SFvsPB_treg.csv"),
                       show_col_types = FALSE)

# create_standard_volcano() reads limma's column names off the frame and takes gene ids from
# ROWNAMES, so the committed table is renamed to fit the plotter. One row per HGNC symbol is what
# makes rownames a safe key; assert it rather than trust it.
DE <- DE %>% filter(!is.na(padj), !is.na(log2FoldChange))
stopifnot(!any(duplicated(DE$gene_symbol)), !any(is.na(DE$gene_symbol)))

de_v <- DE %>%
  transmute(gene_symbol, logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
  as.data.frame()
rownames(de_v) <- de_v$gene_symbol
de_v$gene_symbol <- NULL

# --- counts for the caption, read off the same frame the figure draws ----------------------
sig_fc  <- abs(de_v$logFC) >= LOGFC_CUT
sig_dec <- de_v$adj.P.Val <= FDR_CUT
sig     <- sig_fc & sig_dec
n_sig   <- sum(sig)
n_up    <- sum(sig & de_v$logFC > 0)
n_down  <- sum(sig & de_v$logFC < 0)
n_fdr   <- sum(sig_dec)
# The raw p the FDR cut resolves to — the height of the dashed horizontal rule.
p_bound <- if (n_fdr > 0) max(de_v$P.Value[sig_dec], na.rm = TRUE) else NA_real_

# --- the figure -----------------------------------------------------------------------------
message("[fig] Treg synovial-fluid-versus-blood contrast on the standard volcano ...")

# x_breaks = 2: the contrast spans roughly -7..+9 log2FC, and unit ticks over that range crowd
# the axis past the legibility floor.
p <- create_standard_volcano(
  de_v,
  decision_by     = "fdr",
  p_cutoff        = FDR_CUT,
  fc_cutoff       = LOGFC_CUT,
  top_n           = TOP_N,
  x_breaks        = 2,
  annotate_counts = TRUE,
  # max.overlaps = Inf: a requested label that ggrepel silently drops makes the figure
  # under-report its own label budget.
  max.overlaps    = Inf,
  title    = "Sorted JIA Tregs, synovial fluid versus paired blood",
  # Two lines on purpose: as one line this runs past the right edge of the canvas and the
  # threshold clause is the half that gets cut.
  subtitle = sprintf(paste0("Donor pseudobulk, limma-voom, %d donors paired across both arms\n",
                            "%s genes clear FDR <= %.2g and |log2FC| >= %.1f (%s up, %s down)"),
                     DE$n_paired_donors[1], format(n_sig, big.mark = ","), FDR_CUT, LOGFC_CUT,
                     format(n_up, big.mark = ","), format(n_down, big.mark = ","))
) +
  # project_theme is a COMPLETE theme, so it must come after the plotter's own theming or it
  # would be reset by it; anything element-level goes after this line.
  project_theme(config = CFG) +
  theme(legend.position = "bottom")

# --- the neighbour table --------------------------------------------------------------------
# Carries the significance category the figure colours by and whether the gene is one of the
# named ones, so the table can be read against the canvas rather than merely beside it.
labelled <- de_v %>%
  tibble::rownames_to_column("gene_symbol") %>%
  mutate(side = sign(logFC)) %>%
  filter(side != 0) %>%
  group_by(side) %>%
  arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>%
  ungroup() %>%
  pull(gene_symbol)

fig_tbl <- DE %>%
  transmute(gene_symbol, log2FoldChange, stat, pvalue, padj,
            category = case_when(
              abs(log2FoldChange) >= LOGFC_CUT & padj <= FDR_CUT ~ "p-value & Log2FC",
              abs(log2FoldChange) >= LOGFC_CUT                   ~ "Log2FC",
              padj <= FDR_CUT                                    ~ "p-value",
              TRUE                                               ~ "NS"),
            labelled_on_figure = gene_symbol %in% labelled) %>%
  arrange(pvalue) %>%
  slice_head(n = TABLE_N)

save_overview(
  p, STAGE, "treg_volcano", table = fig_tbl,
  finding = sprintf(paste("Synovial-fluid Tregs carry a reproducible transcriptional program",
                          "against the same donors' blood: %s genes clear FDR <= %.2g and",
                          "|log2FC| >= %.1f, %s up in the joint and %s down, and %s genes clear",
                          "the FDR cut on its own. This contrast is the ranking every enrichment",
                          "in the Treg gate is computed on."),
                    format(n_sig, big.mark = ","), FDR_CUT, LOGFC_CUT,
                    format(n_up, big.mark = ","), format(n_down, big.mark = ","),
                    format(n_fdr, big.mark = ",")),
  script = SCRIPT, fn = "create_standard_volcano",
  config_kv = paste0("thresholds.de_fdr=", FDR_CUT, "; thresholds.de_logfc=", LOGFC_CUT,
                     "; figures.volcano_label_top=", LABEL_TOP),
  input = "03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv",
  how_to_read = paste(
    "x is log2 fold change, synovial fluid over paired blood; y is raw p on a -log10 scale, while",
    "significance is decided on FDR — the axis keeps the per-gene resolution that -log10(FDR)",
    "would collapse. Colour gives four categories: neither cut, fold change only, FDR only, both.",
    sprintf("The dashed horizontal rule is the raw p realising FDR <= %.2g (p <= %.3g), the",
            FDR_CUT, p_bound),
    sprintf("vertical rules |log2FC| >= %.1f, the legend arrows the up/down split of genes",
            LOGFC_CUT),
    sprintf("clearing both. The %d named genes are the %d most significant per side, not the",
            length(labelled), TOP_N),
    sprintf("largest fold changes. The neighbour table holds the top %d by p with their category.",
            TABLE_N),
    "Correlative donor-pseudobulk DE."),
  config = CFG, height = 7.5
)

message("\n[DONE] 03_pseudobulk_volcano_viz complete.")
