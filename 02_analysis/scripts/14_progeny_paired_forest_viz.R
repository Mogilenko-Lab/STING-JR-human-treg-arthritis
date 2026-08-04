#!/usr/bin/env Rscript
# 14_progeny_paired_forest_viz.R: VIZ (no statistics)
# =============================================================================
# The donor-paired reading of the PROGENy footprints in the JIA synovial-fluid-
# versus-paired-blood contrast, drawn as a forest with a 95% interval per point.
#
# This is a DIFFERENT figure from progeny_activity_panel in
# 14_unbiased_enrichment_viz.R, and both are kept. That panel plots the
# single-contrast multivariate-linear-model statistic from progeny_activity.csv,
# one number per population x pathway with no interval, and it puts all three
# populations on the identical y so the three markers overlay. This one plots the
# donor-paired test from progeny_sf_vs_pb.csv: six donors each voting once,
# synovial fluid minus that same donor's own blood, with the interval the paired
# t-test already reported. Every number here is read from the CSV; nothing is
# recomputed.
#
# READABILITY CHOICES, each one a reaction to how the overlaid version failed:
#   * each population gets its own VERTICAL OFFSET inside the pathway row
#     (+0.24 / 0 / -0.24 row units, the offset the rest of this project uses), so
#     three points and three intervals stay separately readable;
#   * significance is a FILLED versus OPEN marker as well as a position, so the
#     distinction survives greyscale printing and a colour-blind reader;
#   * faint separators run between pathway rows so a reader can track one pathway
#     across a wide panel without losing the row;
#   * the exact numbers stay in the same-stem source CSV rather than being printed
#     forty-two times on the canvas.
#
# Input  (03_results/14_unbiased_enrichment/tables/):
#   progeny_sf_vs_pb.csv     42 rows = 3 sorted populations x 14 PROGENy pathways
#
# Output (03_results/14_unbiased_enrichment/):
#   figures/_overview/progeny_paired_forest.{pdf,png}
#   tables/_overview/progeny_paired_forest.csv
#   README.md caption                              (via save_overview)
#
# Run from the compartment root, AFTER 14_unbiased_enrichment.R:
#   Rscript 02_analysis/scripts/14_progeny_paired_forest_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview

STAGE  <- "14_unbiased_enrichment"
SCRIPT <- "02_analysis/scripts/14_progeny_paired_forest_viz.R"
STEM   <- "progeny_paired_forest"
TDIR   <- file.path(FIG_CFG$paths$results %||% "03_results/", STAGE,
                    FIG_CFG$paths$stage_tables_subdir %||% "tables")

FDR   <- as.numeric(FIG_CFG$thresholds$gsea_fdr %||% 0.05)
OI    <- FIG_CFG$colors$okabe_ito
PT    <- as.numeric(FIG_CFG$figures$point_size %||% 2.4)
LW    <- as.numeric(FIG_CFG$figures$line_width %||% 1.0)
THEME <- project_theme(config = FIG_CFG)

POP_LEVELS <- c("Treg", "Tcon", "CD8")
## Vertical offset per sorted population inside one pathway row. Three points on a
## shared y overlay each other and their intervals merge into one bar; 0.24 row units
## is the separation the rest of this project's per-pathway panels use.
POP_DY <- c(Treg = 0.24, Tcon = 0.0, CD8 = -0.24)
POP_COL <- c(Treg = OI$vermillion, Tcon = OI$blue, CD8 = OI$bluish_green)

#' Render an FDR for prose: fixed above 0.01, else one-digit scientific. The cut sits at
#' 0.01 rather than 0.001 because a fixed-decimal 0.00115 renders as "0.001" and loses
#' the leading digit a reader compares against a neighbouring pathway.
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p >= 0.01, sprintf("%.3f", p), sprintf("%.1e", p)))

#' Hard-wrap a title / subtitle / caption so it cannot be clipped at the canvas edge.
wrap_at <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

p_path <- file.path(TDIR, "progeny_sf_vs_pb.csv")
if (!file.exists(p_path))
  stop("[14_forest] ", p_path, " not found. Run 02_analysis/scripts/14_unbiased_enrichment.R first.")

pg <- readr::read_csv(p_path, show_col_types = FALSE, progress = FALSE)
stopifnot("progeny_sf_vs_pb.csv must carry its 95% interval columns" =
            all(c("ci_low", "ci_high") %in% names(pg)),
          "every row needs a populated interval" =
            all(!is.na(pg$ci_low) & !is.na(pg$ci_high)))

## Rows are ordered by the Treg paired difference, which the caption states, so a
## reader knows the ordering is one population's and not a consensus.
ord <- pg |> dplyr::filter(population == "Treg") |>
  dplyr::arrange(mean_difference) |> dplyr::pull(pathway_name)
stopifnot("the Treg arm must cover every pathway drawn" =
            setequal(ord, unique(pg$pathway_name)))

pg <- pg |>
  dplyr::mutate(population  = factor(population, levels = POP_LEVELS),
                pathway_name = factor(pathway_name, levels = ord),
                significant  = !is.na(padj) & padj < FDR,
                y = as.integer(pathway_name) + POP_DY[as.character(population)])

N_ROW  <- length(ord)
N_PAIR <- sort(unique(pg$n_paired_donors))
xr     <- range(c(pg$ci_low, pg$ci_high))
## The widest interval reaches +15.6, so a thin pad leaves that whisker flush against
## the panel edge and it reads as clipped. Pad both sides and give the axis breaks every
## five units so a reader can place a point without counting from zero.
pad    <- diff(xr) * 0.07
xbrk   <- seq(-15, 15, by = 5)

p <- ggplot(pg, aes(x = mean_difference, y = y, colour = population)) +
  geom_hline(yintercept = seq_len(N_ROW - 1L) + 0.5, linewidth = 0.3, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = LW * 0.8, alpha = 0.8) +
  geom_point(aes(shape = significant), size = PT * 1.7, fill = "white", stroke = LW * 1.1) +
  scale_colour_manual(values = POP_COL, name = "sorted population") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     breaks = c(TRUE, FALSE),
                     labels = c(`TRUE` = sprintf("FDR < %.2g", FDR),
                                `FALSE` = sprintf("FDR \u2265 %.2g", FDR)),
                     name = "donor-paired test", drop = FALSE) +
  scale_y_continuous(breaks = seq_len(N_ROW), labels = ord,
                     limits = c(0.55, N_ROW + 0.45), expand = c(0, 0)) +
  scale_x_continuous(limits = c(xr[1] - pad, xr[2] + pad), breaks = xbrk) +
  guides(colour = guide_legend(override.aes = list(size = PT * 1.5, shape = 16)),
         shape  = guide_legend(override.aes = list(size = PT * 1.5, colour = "grey20"))) +
  labs(
    title = "PROGENy footprints, paired synovial fluid minus blood",
    subtitle = wrap_at(paste("Mean within-donor difference in multivariate-linear-model activity,",
                             "synovial fluid minus that same donor's blood, with the 95% interval",
                             "of the paired t-test. Positive is higher in synovial fluid."), 118),
    x = "mean paired difference in PROGENy activity score",
    y = NULL,
    caption = wrap_at(paste("Six paired donors per population. A footprint is inferred from",
                            "target-gene expression, so it reads downstream transcription and",
                            "the reading stays correlative."),
                      118)) +
  THEME +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

tbl <- pg |>
  dplyr::select(population, contrast, pathway_name, n_paired_donors, mean_sf, mean_pb,
                mean_difference, ci_low, ci_high, t_statistic, pvalue, padj, direction,
                significant) |>
  dplyr::arrange(population, dplyr::desc(mean_difference))

## Checks named in the task, read back off the table rather than restated by hand, so a
## regenerated table that moved would move this sentence with it.
gv <- function(pop, pw, col) pg[[col]][pg$population == pop & pg$pathway_name == pw][1]
egfr  <- gv("Treg", "EGFR", "mean_difference");    egfr_p <- gv("Treg", "EGFR", "padj")
nfkb  <- gv("Treg", "NFkB", "mean_difference");    nfkb_p <- gv("Treg", "NFkB", "padj")
hyp   <- gv("Treg", "Hypoxia", "mean_difference"); hyp_p  <- gv("Treg", "Hypoxia", "padj")
tnfa  <- gv("Treg", "TNFa", "mean_difference");    tnfa_p <- gv("Treg", "TNFa", "padj")
n_sig <- pg |> dplyr::filter(population == "Treg", significant) |> nrow()

save_overview(
  p, STAGE, STEM, table = tbl,
  finding = sprintf(paste0("Tested one donor at a time, %d of the fourteen PROGENy footprints ",
                           "separate the JIA synovial Treg pool from paired blood at FDR < %.2g. ",
                           "EGFR is the largest at %+.1f (FDR %s), followed by NFkB at %+.1f (FDR ",
                           "%s) and Hypoxia at %+.1f (FDR %s) across all six paired donors, while ",
                           "TNFa sits lower on the synovial-fluid side at %+.2f (FDR %s), so the ",
                           "footprints that separate split across both directions."),
                    n_sig, FDR, egfr, fmt_p(egfr_p), nfkb, fmt_p(nfkb_p),
                    hyp, fmt_p(hyp_p), tnfa, fmt_p(tnfa_p)),
  script = SCRIPT, fn = "main",
  config_kv = sprintf("progeny.organism=%s; progeny.top=%s; progeny.minsize=%s; gsea_fdr=%s; population_offset=%.2f row units",
                      FIG_CFG$unbiased_enrichment$progeny$organism,
                      FIG_CFG$unbiased_enrichment$progeny$top,
                      FIG_CFG$unbiased_enrichment$progeny$minsize, FDR, POP_DY[["Treg"]]),
  input = "03_results/14_unbiased_enrichment/tables/progeny_sf_vs_pb.csv",
  how_to_read = paste(
    "One row per PROGENy pathway, rows ordered by the Treg paired difference so the",
    "ordering belongs to one population. Inside a row the three sorted populations are",
    "offset vertically, Treg above, Tcon centre, CD8 below, so their markers and",
    "intervals stay separately readable; colour repeats the same key. Horizontal",
    "position is the mean within-donor difference in activity, synovial fluid minus that",
    "donor's own blood, and the bar through it is the 95% interval of the paired t-test.",
    "A filled marker reaches FDR < 0.05 in that population and an open marker sits above",
    "that threshold, so the distinction survives a greyscale print. Right of zero the",
    "footprint is higher in synovial fluid, left of it higher in blood. Two limits on the",
    "reading. The",
    "activity score is computed on expression centred within a population, so a",
    "difference is comparable between pathways of the same population and carries no",
    "meaning compared across populations; a Treg point sitting further right than a CD8",
    "point on the same row is not a between-population effect size. And six pairs is a",
    "small n, so an open marker here leaves the question open rather than settling it.",
    "The six donors are the same count in every population but not the same six people:",
    "Treg pairs JIA_patient_3 and skips 5, Tcon and CD8 do the reverse. Correlative; a",
    "footprint is inferred from target-gene expression."),
  config = FIG_CFG, width = 12, height = 9)

message(sprintf("[14_forest] wrote %s with its source table and README caption (%d pairs per test)",
                STEM, N_PAIR[1]))
