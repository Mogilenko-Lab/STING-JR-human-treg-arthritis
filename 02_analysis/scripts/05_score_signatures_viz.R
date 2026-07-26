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
# theme-consistent). Keyed BY LEGEND LABEL, which is the value the colour
# aesthetic actually carries — see `apply_set_palette()` for why a positional
# vector cannot express this mapping safely.
SET_PAL <- c("WT_heat up" = "#A6611A", "WT_heat down" = "#2166AC")

#' Re-assert the up/down colour mapping BY NAME on every panel.
#'
#' Defence in depth, kept deliberately. The original defect was upstream: the
#' toolkit plotter forwarded `palette` UNNAMED into `enrichplot::gseaplot2()`,
#' which maps `colour = Description` and calls `scale_color_manual(values =)`.
#' Unnamed values match by POSITION against levels ggplot has sorted
#' ALPHABETICALLY — "WT_heat down" precedes "WT_heat up" — so a positional
#' palette in GENE_SETS order put brown on the down curve. The `labels` lookup
#' is named, so the legend TEXT stayed correct while the colours swapped, and it
#' read as a deliberate bad choice rather than a bug.
#'
#' The toolkit now re-keys the palette to the plotted label itself and owns the
#' scale (RNAseq-toolkit `fix(GSEA): key running-sum palette by plotted label`),
#' so this is no longer load-bearing — verified idempotent: toolkit-only and
#' toolkit-plus-this both bind brown to the up curve. It stays because the
#' parent pins the toolkit by gitlink and that fix carries no release tag yet,
#' so an older pin would silently reintroduce the swap. Retire it once the
#' compartment pins a toolkit tag containing the fix. `limits` fixes legend
#' order to up-then-down.
apply_set_palette <- function(p, ids) {
  labs_in_order <- unname(SET_LABELS[ids])
  suppressMessages(
    p & ggplot2::scale_color_manual(values = SET_PAL, limits = labs_in_order,
                                    name = NULL))
}

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
  "The Treg up-curve carries the claim; Tcon and CD8 test whether it is",
  "Treg-selective.")
HOW_TO_READ <- paste(
  "Top panel = weighted running enrichment score (ES) walking the ranked list",
  "from SF-enriched (left) to PB-enriched (right); a positive, left-shifted peak",
  "= SF enrichment. Middle rug = gene-set member positions; bottom = the signed",
  "Wald ranking metric. Two curves per panel, same colour in curve and rug:",
  "WT_heat up = warm brown, WT_heat down = cool blue.",
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
    palette         = unname(SET_PAL[SET_LABELS[ids]]),
    labels          = SET_LABELS[ids],
    es_ylim         = ylim,
    panel_heights   = heights,
    legend_position = "right",
    base_theme      = base_thm,
    title           = sprintf("%s — WT_heat running enrichment (SF vs PB)", pop))
  p <- apply_set_palette(p, ids)

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
