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
source("02_analysis/helpers/source_hash_manifest.R")
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")

STAGE  <- "05_scoring"
SCRIPT <- "02_analysis/scripts/05_score_signatures_viz.R"
TDIR   <- "03_results/05_scoring/tables"

POP_TAG    <- c(Treg = "treg", Tcon = "tcon", CD8 = "cd8")
GENE_SETS  <- c("WT_heat_up", "WT_heat_down")
SET_LABELS <- c(WT_heat_up = "WT_heat up", WT_heat_down = "WT_heat down")
SET_ARM    <- c(WT_heat_up = "up", WT_heat_down = "down")
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
apply_set_palette <- function(p, ids, labs) {
  labs_in_order <- unname(labs[ids])
  pal <- stats::setNames(unname(SET_PAL[SET_LABELS[ids]]), labs_in_order)
  suppressMessages(
    p & ggplot2::scale_color_manual(values = pal, limits = labs_in_order,
                                    name = NULL))
}

#' Nominal size of each frozen mouse arm — a line count, not a statistic.
#'
#' An effective set size only means something against the nominal one, so the
#' legend carries both. Read from the frozen mouse->human projection contract that
#' the compute step scored, so the denominator cannot drift from the numerator.
nominal_set_sizes <- function() {
  contract <- FIG_CFG$paths$signature_contract %||%
    "../mouse_anchor/03_results/human_projection/"
  sig_dir <- file.path(contract, "signatures", "WT_heat")
  vapply(GENE_SETS, function(id) {
    path <- file.path(sig_dir, paste0(id, ".txt"))
    if (!file.exists(path))
      stop("[05_viz_R] frozen signature not found: ", path)
    verify_source_hash(
      path, id, file.path(TDIR, "source_hash_manifest.csv"),
      root = normalizePath("..", mustWork = FALSE))
    genes <- trimws(readLines(path, warn = FALSE))
    length(unique(genes[nzchar(genes)]))
  }, FUN.VALUE = integer(1))
}

#' Render an FDR for an in-figure label: fixed below 3 decimals, else scientific.
fmt_fdr <- function(p) if (is.na(p)) "FDR n/a" else if (p >= 0.001)
  sprintf("FDR %.3f", p) else sprintf("FDR %.0e", p)

ylim    <- as.numeric(unlist(FIG_CFG$figures$running_sum_ylim     %||% c(-1, 1)))
heights <- as.numeric(unlist(FIG_CFG$figures$running_sum_heights  %||% c(2.4, 0.7, 0.9)))
base_thm <- if (exists("project_theme")) project_theme(config = FIG_CFG) else theme_classic(base_size = 14)

# Own the running-sum figure namespace: drop the retired single 3-population
# figure (wt_heat_running_sum.*) and any stale per-pop files before rewriting.
if (exists("purge_figures"))
  purge_figures(STAGE, "wt_heat_running_sum", overview = TRUE, config = FIG_CFG)

#' Build this population's finding from THIS population's own summary rows.
#'
#' A caption may quote only its same-stem table, and this figure's same-stem table
#' holds one population. The retired caption quoted all three populations' NES from
#' the dot plot's table instead, which made a three-population claim from a
#' one-population panel; the cross-population comparison belongs on the dot plot
#' and this caption points there rather than restating it.
finding_for <- function(pop, tbl, nominal) {
  parts <- vapply(seq_len(nrow(tbl)), function(i) {
    r <- tbl[i, ]
    sprintf("the %s arm reaches NES %+.4f at %s with %d of its %d genes in the ranked list",
            SET_ARM[[r$ID]], r$NES, fmt_fdr(r$p.adjust), r$setSize, nominal[[r$ID]])
  }, FUN.VALUE = character(1))
  paste0(
    "In ", pop, " ", paste(parts, collapse = ", and "), ". ",
    "The curve gives the place along this population's synovial-fluid-versus-blood ",
    "ranking where each arm concentrates. The cross-population comparison — whether ",
    "one sorted population separates more than another — is read off the ordered NES ",
    "dot plot.")
}

HOW_TO_READ <- paste(
  "One population per panel, showing the donor-pseudobulk fgsea result behind the",
  "confirmatory answer. The top trace walks from SF-enriched to PB-enriched genes;",
  "a positive left peak indicates SF enrichment. The middle rug marks set members",
  "and the bottom shows the signed moderated-t ranking. Warm brown is the up arm and",
  "cool blue the down arm. Legends report effective and nominal size, NES, and FDR.",
  "The shared [-1, 1] enrichment-score range supports shape comparison. Read the",
  "cross-population result from the ordered NES dot plot, which establishes the",
  "pan-T pattern. Display of compute output; correlative.")
CONFIG_KV <- paste0("gsea_min_size=", FIG_CFG$thresholds$gsea_min_size %||% 5,
                    "; gsea_max_size=", FIG_CFG$thresholds$gsea_max_size %||% 500,
                    "; running_sum_ylim=[", ylim[1], ",", ylim[2], "]",
                    "; engine=clusterProfiler::GSEA(by=fgsea)")

NOMINAL <- nominal_set_sizes()

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

  # Source table = the gseaResult summary rows behind THIS figure.
  res <- g@result
  tbl <- res[res$ID %in% ids,
             c("ID", "Description", "setSize", "enrichmentScore", "NES",
               "pvalue", "p.adjust", "core_enrichment")]
  tbl$population  <- pop
  tbl$n_nominal   <- unname(NOMINAL[tbl$ID])
  tbl <- tbl[match(ids, tbl$ID), ]

  # Effective set size travels with every NES, on the FACE. The legend label is the
  # only string on this plot wide enough to hold it, and the toolkit wraps labels
  # only past `max_name_length`, so the limit is raised past the longest label
  # rather than left to truncate one.
  labs <- stats::setNames(
    vapply(ids, function(id) {
      r <- tbl[tbl$ID == id, ][1, ]
      sprintf("%s   %d of %d genes ranked,  NES %+.2f,  %s",
              SET_LABELS[[id]], r$setSize, r$n_nominal, r$NES, fmt_fdr(r$p.adjust))
    }, FUN.VALUE = character(1)),
    ids)

  p <- gsea_running_sum_plot(
    g,
    gene_set_ids    = ids,
    palette         = stats::setNames(unname(SET_PAL[SET_LABELS[ids]]), ids),
    labels          = labs,
    max_name_length = max(nchar(labs)) + 1L,
    es_ylim         = ylim,
    panel_heights   = heights,
    legend_position = "right",
    base_theme      = base_thm,
    title           = sprintf(
      "%s — where the mouse 39 °C-derived arms sit\nalong this population's SF-vs-PB ranking",
      pop))
  p <- apply_set_palette(p, ids, labs)

  save_overview(
    p, STAGE, sprintf("wt_heat_running_sum_%s", tag), table = tbl,
    finding = finding_for(pop, tbl, NOMINAL), script = SCRIPT, fn = "main",
    config_kv = CONFIG_KV,
    input = sprintf("03_results/05_scoring/tables/gsea_pseudobulk_%s.rds", tag),
    how_to_read = HOW_TO_READ,
    config = FIG_CFG, width = 11, height = 7)
  n_written <- n_written + 1L
}

if (n_written == 0L)
  stop("[05_viz_R] No gseaResult RDS found in ", TDIR, ". Run 05_score_signatures.py first.")
message(sprintf("[05_viz_R] wrote %d per-population running-sum overview(s)", n_written))
