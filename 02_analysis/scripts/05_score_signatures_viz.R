#!/usr/bin/env Rscript
# 05_score_signatures_viz.R — VIZ (no statistics). Per-population running-sum.
# =====================================================================
# Drives the CANONICAL RNAseq-toolkit running-sum plotter
# (01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R)
# off the clusterProfiler `gseaResult` S4 objects written by the compute step
# (helpers/fgsea_prerank.R -> gsea_pseudobulk_{tag}.rds). The plotter overlays
# N gene SETS from ONE ranked list, so we emit ONE figure PER POPULATION, each
# overlaying WT_heat_up AND WT_heat_down. (The retired hand-rolled 3-population
# single-axis ggplot is gone: it cannot mix a set-overlay with a pop-overlay.)
#
# Input  (03_results/05_scoring/tables/):
#   gsea_pseudobulk_{treg,tcon,cd8}.rds   ← clusterProfiler gseaResult S4
#
# Output (03_results/05_scoring/):
#   figures/_overview/wt_heat_running_sum_{treg,tcon,cd8}.{pdf,png}
#   tables/_overview/wt_heat_running_sum_{treg,tcon,cd8}.csv   ← same-stem source table
#   README.md captions (via save_overview)
#
# Run from the compartment root, AFTER 05_score_signatures.py:
#   Rscript 02_analysis/scripts/05_score_signatures_viz.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("02_analysis/helpers/figure_style.R")
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")

STAGE  <- "05_scoring"
SCRIPT <- "02_analysis/scripts/05_score_signatures_viz.R"
TDIR   <- "03_results/05_scoring/tables"

POP_TAG    <- c(Treg = "treg", Tcon = "tcon", CD8 = "cd8")
GENE_SETS  <- c("WT_heat_up", "WT_heat_down")
SET_LABELS <- c(WT_heat_up = "WT_heat up", WT_heat_down = "WT_heat down")
# Diverging cue: heat-up = warm = brown, heat-down = cool = blue (semantic,
# theme-consistent). UNNAMED on purpose — the toolkit plotter requires an
# unnamed color vector; order follows GENE_SETS (up first, then down).
SET_PAL <- c("#A6611A", "#2166AC")

ylim    <- as.numeric(unlist(FIG_CFG$figures$running_sum_ylim     %||% c(-1, 1)))
heights <- as.numeric(unlist(FIG_CFG$figures$running_sum_heights  %||% c(2.4, 0.7, 0.9)))
base_thm <- if (exists("project_theme")) project_theme(config = FIG_CFG) else theme_classic(base_size = 14)

# Own the running-sum figure namespace: drop the retired single 3-population
# figure (wt_heat_running_sum.*) and any stale per-pop files before rewriting.
if (exists("purge_figures"))
  purge_figures(STAGE, "wt_heat_running_sum", overview = TRUE, config = FIG_CFG)

FINDING <- paste(
  "Per-population leading-edge view: where the mouse 39 °C WT_heat up- and",
  "down-programs concentrate along each population's SF-vs-PB pseudobulk ranking.",
  "The Treg up-curve is the go/no-go money panel; Tcon/CD8 are the",
  "Treg-specificity controls.")
HOW_TO_READ <- paste(
  "Top panel = weighted running enrichment score (ES) walking the ranked list",
  "from SF-enriched (left) to PB-enriched (right); a positive, left-shifted peak",
  "= SF enrichment. Middle rug = gene-set member positions; bottom = the signed",
  "Wald ranking metric. Two curves per panel: WT_heat up (warm) and down (cool).",
  "ES y clamped to [-1, 1] for cross-population comparability. Display of compute",
  "output (clusterProfiler gseaResult); correlative, not causal.")
CONFIG_KV <- paste0("gsea_min_size=", FIG_CFG$thresholds$gsea_min_size %||% 5,
                    "; gsea_max_size=", FIG_CFG$thresholds$gsea_max_size %||% 500,
                    "; running_sum_ylim=[", ylim[1], ",", ylim[2], "]",
                    "; engine=clusterProfiler::GSEA(by=fgsea)")

n_written <- 0L
for (pop in names(POP_TAG)) {
  tag <- POP_TAG[[pop]]
  rds <- file.path(TDIR, sprintf("gsea_pseudobulk_%s.rds", tag))
  if (!file.exists(rds)) {
    warning(sprintf("[05_viz_R] gseaResult RDS not found, skipping %s: %s", pop, rds))
    next
  }
  g   <- readRDS(rds)
  ids <- intersect(GENE_SETS, g@result$ID)
  if (length(ids) == 0) {
    warning(sprintf("[05_viz_R] no WT_heat sets in %s result, skipping", pop))
    next
  }

  p <- gsea_running_sum_plot(
    g,
    gene_set_ids    = ids,
    palette         = SET_PAL[seq_along(ids)],
    labels          = SET_LABELS[ids],
    es_ylim         = ylim,
    panel_heights   = heights,
    legend_position = "right",
    base_theme      = base_thm,
    title           = sprintf("%s — WT_heat running enrichment (SF vs PB)", pop))

  # Source table = the gseaResult summary rows behind THIS figure.
  res <- g@result
  tbl <- res[res$ID %in% ids,
             c("ID", "Description", "setSize", "enrichmentScore", "NES",
               "pvalue", "p.adjust", "core_enrichment")]
  tbl$population <- pop

  save_overview(
    p, STAGE, sprintf("wt_heat_running_sum_%s", tag), table = tbl,
    finding = FINDING, script = SCRIPT, fn = "main",
    config_kv = CONFIG_KV,
    input = sprintf("03_results/05_scoring/tables/gsea_pseudobulk_%s.rds", tag),
    how_to_read = HOW_TO_READ,
    config = FIG_CFG, height = 7)
  n_written <- n_written + 1L
}

if (n_written == 0L)
  stop("[05_viz_R] No gseaResult RDS found in ", TDIR, ". Run 05_score_signatures.py first.")
message(sprintf("[05_viz_R] wrote %d per-population running-sum overview(s)", n_written))
