#!/usr/bin/env Rscript
# 14_sweep_named_sets_viz.R: VIZ (no statistics)
# =============================================================================
# The closing panel of this compartment's narrative. Every earlier figure tested a
# gene set that was brought to the data. This one scores no favourite: every set in
# eleven collections was run on the same donor-level ranked lists, and the sets the
# narrative names are placed inside that full distribution, per sorted population, so
# a reader sees both the effect and how ordinary or extraordinary it is.
#
# WHICH SETS ARE DRAWN IS A COMMITTED DECISION, NOT A CHOICE MADE HERE. The selection
# and the reason for every row live in tables/sweep_named_sets.csv, written by
# 14_sweep_named_sets.R, which also audits the substring trap that a naive search for
# "STING" falls into. This script reads that table and plots it.
#
# FOUR THINGS THE PANEL AND ITS CAPTION CARRY, because leaving any of them out turns
# the figure into an overclaim by omission:
#   1. all three sorted populations, since sting_specific_up reaches pooled FDR 0.018
#      in Tcon while missing in Treg and CD8, and a Treg-only panel would read as a
#      settled negative for a set that is significant in one of the three;
#   2. KO_heat_up beside WT_heat_up, since the cGAS-knockout comparator beats the
#      primary arm on pooled FDR in all three populations and shares 182 of its genes;
#   3. set size, on the marker area and in the caption's baseline rates, because how
#      often a set of a given size reaches significance at all differs by nearly
#      sevenfold between the band HALLMARK_HYPOXIA sits in and the band the cGAS-STING
#      sets sit in;
#   4. how many sets reach significance at all in each population, so a headline is
#      read against the size of the family it came out of.
#
# The grey cloud on the bottom row of each panel is every set tested in that
# population. It is the calibration: a marker far to the right inside a dense cloud is
# strong and ordinary at the same time.
#
# Input  (03_results/14_unbiased_enrichment/tables/):
#   gsea_all.csv                 the full sweep, drawn as the background cloud
#   sweep_named_sets.csv         the committed selection and its per-row reasons
#   sweep_named_sets_stats.csv   per population x named set statistics and ranks
#   sweep_setsize_baseline.csv   pooled-significance rate by set-size band
#
# Output (03_results/14_unbiased_enrichment/):
#   figures/_overview/named_sets_in_sweep.{pdf,png}
#   tables/_overview/named_sets_in_sweep.csv
#   README.md caption                                  (via save_overview)
#
# Run from the compartment root, AFTER 14_sweep_named_sets.R:
#   Rscript 02_analysis/scripts/14_sweep_named_sets_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview

STAGE  <- "14_unbiased_enrichment"
SCRIPT <- "02_analysis/scripts/14_sweep_named_sets_viz.R"
STEM   <- "named_sets_in_sweep"
TDIR   <- file.path(FIG_CFG$paths$results %||% "03_results/", STAGE,
                    FIG_CFG$paths$stage_tables_subdir %||% "tables")

FDR    <- as.numeric(FIG_CFG$thresholds$gsea_fdr %||% 0.05)
NESCAP <- as.numeric(FIG_CFG$figures$nes_cap %||% 3.5)
OI     <- FIG_CFG$colors$okabe_ito
PT     <- as.numeric(FIG_CFG$figures$point_size %||% 2.4)
LBL    <- as.numeric(FIG_CFG$figures$label_size %||% 4)
LW     <- as.numeric(FIG_CFG$figures$line_width %||% 1.0)
RDPI   <- as.numeric(FIG_CFG$figures$rasterized_dpi %||% FIG_CFG$figures$dpi %||% 600)
THEME  <- project_theme(config = FIG_CFG)

POP_LEVELS <- c("Treg", "Tcon", "CD8")
set.seed(as.integer(FIG_CFG$thresholds$gsea_seed %||% 123))   # the background jitter

read_tbl <- function(f) {
  p <- file.path(TDIR, f)
  if (!file.exists(p))
    stop("[14_named] ", p, " not found. Run 02_analysis/scripts/14_sweep_named_sets.R first.")
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)
}

sweep    <- read_tbl("gsea_all.csv")
named    <- read_tbl("sweep_named_sets.csv")
stats    <- read_tbl("sweep_named_sets_stats.csv")
baseline <- read_tbl("sweep_setsize_baseline.csv")

#' Render an FDR for an on-face label: fixed above 0.01, else one-digit scientific.
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p >= 0.01, sprintf("%.3f", p), sprintf("%.1e", p)))
#' Thousands separator, for counts a reader is meant to compare at a glance.
cma <- function(x) format(x, big.mark = ",", trim = TRUE)
#' Hard-wrap a title / subtitle / caption so it cannot be clipped at the canvas edge.
wrap_at <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

# =============================================================================
# 1. Row geometry: eleven named rows, plus one row holding the whole sweep
# =============================================================================
## Threads run top to bottom in the order the selection table fixed, and within a
## thread the rows are ordered by the Treg NES, so the ordering belongs to one
## population and the caption says which. The background cloud gets row 0, below a
## heavier separator, because it is the reference distribution rather than a member of
## the comparison.
sel_order <- stats |>
  dplyr::filter(population == "Treg") |>
  dplyr::arrange(thread_order, dplyr::desc(nes)) |>
  dplyr::pull(display_label)
N_SET  <- length(sel_order)
YPOS   <- setNames(seq(N_SET, 1L), sel_order)

## MSigDB identifiers carry no spaces, so str_wrap on the raw id leaves a single
## 58-character line that pushes the panel off the canvas; the selection table already
## turned underscores into spaces for the all-caps ids, so wrapping works and nothing
## is truncated.
wrap_lbl <- function(x, width = 32) stringr::str_wrap(x, width = width)

pts <- stats |>
  dplyr::mutate(population = factor(population, levels = POP_LEVELS),
                y = YPOS[display_label],
                nes_plot = pmax(pmin(nes, NESCAP), -NESCAP))
drawn <- pts |> dplyr::filter(tested)

## Thread boundaries, so a reader can see the four comparisons as blocks.
thread_rows <- pts |> dplyr::filter(population == "Treg") |>
  dplyr::group_by(thread_order, thread) |>
  dplyr::summarise(y_min = min(y), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(y_min))
sep_y <- thread_rows$y_min[-nrow(thread_rows)] - 0.5

bg <- sweep |>
  dplyr::mutate(population = factor(population, levels = POP_LEVELS),
                nes_plot = pmax(pmin(nes, NESCAP), -NESCAP))

# =============================================================================
# 2. The per-panel annotation gutter
# =============================================================================
## Two lines per row, not one: a single "FDR 3.7e-12, rank 31 of 11,236" ran past the
## gutter and lost its own tail to the canvas edge in the first render. The denominator
## is dropped from every row and stated once per panel on the background row instead,
## since it is the same number for all eleven. The rank is on pooled FDR and only on
## pooled FDR, which the caption names, because the same set ranks very differently by
## NES and mixing the two in one sentence misreads both.
ann_sets <- pts |>
  dplyr::transmute(population, y,
                   label = ifelse(tested,
                                  sprintf("FDR %s\nrank %s", fmt_p(padj_pooled),
                                          cma(rank_padj_pooled)),
                                  sprintf("not tested\nin %s", population)))
## The count that matters on the background row differs per population, so it is
## annotated inside each panel; a shared axis label could only carry a sum that would
## be wrong for every panel beside it.
ann_bg <- stats |>
  dplyr::group_by(population) |>
  dplyr::summarise(n_tests_pooled = dplyr::first(n_tests_pooled),
                   n_sig_pooled = dplyr::first(n_sig_pooled), .groups = "drop") |>
  dplyr::transmute(population = factor(population, levels = POP_LEVELS), y = 0,
                   label = sprintf("%s of %s\nat FDR < %.2g",
                                   cma(n_sig_pooled), cma(n_tests_pooled), FDR))
ann <- dplyr::bind_rows(ann_sets, ann_bg) |> dplyr::mutate(x = NESCAP * 1.12)

THREAD_COL <- c(OI$vermillion, OI$blue, OI$reddish_purple, OI$bluish_green)
names(THREAD_COL) <- thread_rows$thread

## Alpha and point size are deliberately low. At alpha 0.18 the cloud saturated into a
## solid grey slab running out to the clamp, which hid the one thing the row exists to
## show: the bulk of the distribution sits between NES 0.6 and 1.9 in Treg and only
## about one set in a hundred passes 2.25, so a marker out at 2.6 is genuinely in the
## tail. A translucent cloud makes that thinning visible instead of asserting it.
bg_layer <- geom_jitter(data = bg, aes(x = nes_plot, y = 0), inherit.aes = FALSE,
                        height = 0.32, width = 0, size = PT * 0.13,
                        alpha = 0.08, colour = "grey30")
## Only the CLOUD is rasterised, not the named markers: at nearly 34,000 points the
## cloud would embed one vector glyph per set and make the PDF unopenable, while the
## eleven markers per panel are what a collaborator will want to select and move in a
## vector editor.
if (requireNamespace("ggrastr", quietly = TRUE))
  bg_layer <- ggrastr::rasterise(bg_layer, dpi = RDPI)

p <- ggplot(drawn, aes(x = nes_plot, y = y)) +
  bg_layer +
  geom_hline(yintercept = 0.5, linewidth = 0.5, linetype = "22", colour = "grey45") +
  geom_hline(yintercept = sep_y, linewidth = 0.3, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_point(aes(colour = thread, shape = significant_pooled, size = set_size),
             fill = "white", stroke = LW * 1.2) +
  geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, size = LBL * 0.78, lineheight = 0.9, colour = "grey25") +
  scale_colour_manual(values = THREAD_COL, breaks = thread_rows$thread,
                      name = "comparison thread") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     breaks = c(TRUE, FALSE),
                     labels = c(`TRUE` = sprintf("pooled FDR < %.2g", FDR),
                                `FALSE` = sprintf("pooled FDR \u2265 %.2g", FDR)),
                     name = "pooled significance", drop = FALSE) +
  scale_size_continuous(range = c(PT * 0.8, PT * 2.6), breaks = c(5, 20, 60, 140),
                        name = "genes in the ranked list") +
  scale_x_continuous(limits = c(-NESCAP * 1.03, NESCAP * 1.66),
                     breaks = seq(-3, 3, by = 1.5)) +
  scale_y_continuous(breaks = c(0, unname(YPOS)),
                     labels = c("all sets tested", wrap_lbl(names(YPOS))),
                     limits = c(-0.55, N_SET + 0.55), expand = c(0, 0)) +
  facet_wrap(~ population, nrow = 1) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = PT * 1.5, shape = 16)),
         shape  = guide_legend(order = 2, override.aes = list(size = PT * 1.5, colour = "grey20")),
         size   = guide_legend(order = 3, override.aes = list(colour = "grey20", shape = 16))) +
  labs(
    title = "Named sets against every set tested, by population",
    subtitle = wrap_at(paste("Positive is enriched toward synovial fluid, negative toward",
                             "paired blood. Marker area is the number of the set's genes",
                             "present in that population's ranked list."), 132),
    x = sprintf("normalized enrichment score, clamped to ±%.1f", NESCAP),
    y = NULL,
    caption = wrap_at(paste("The bottom row of each panel is every set tested in that",
                            "population, one grey point per set, so a marker can be read",
                            "against the distribution it came out of. Correlative: a set",
                            "enriching says its gene content moves with the synovial-fluid side",
                            "of this contrast, and carries no claim that the program it is named",
                            "for is present."), 124)) +
  THEME +
  # Legends STACKED: laid out horizontally the three guides span the whole canvas and
  # the last label ends flush against the right edge.
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey92"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(lineheight = 0.9))

# =============================================================================
# 3. Source table and caption
# =============================================================================
tbl <- stats |>
  dplyr::left_join(named |> dplyr::select(pathway_id, why_included, source_collection,
                                          n_genes_shared_with_WT_heat_up),
                   by = "pathway_id") |>
  dplyr::select(population, thread, pathway_id, display_label, source_collection,
                tested, direction, nes, pvalue, padj, padj_pooled, significant_pooled,
                set_size, leading_edge_size, rank_padj_pooled, rank_nes_signed,
                rank_nes_abs, n_tests_pooled, n_sig_pooled,
                n_genes_shared_with_WT_heat_up, why_included) |>
  dplyr::arrange(population, thread, padj_pooled)

gv <- function(pop, pid, col) stats[[col]][stats$population == pop & stats$pathway_id == pid][1]
bv <- function(pop, band, col) baseline[[col]][baseline$population == pop &
                                                baseline$band_label == band][1]

hyp_nes  <- gv("Treg", "HALLMARK_HYPOXIA", "nes")
hyp_p    <- gv("Treg", "HALLMARK_HYPOXIA", "padj_pooled")
hyp_rank <- gv("Treg", "HALLMARK_HYPOXIA", "rank_padj_pooled")
wt_p     <- gv("Treg", "WT_heat_up", "padj_pooled")
wt_rank  <- gv("Treg", "WT_heat_up", "rank_padj_pooled")
ko_p     <- gv("Treg", "KO_heat_up", "padj_pooled")
ko_rank  <- gv("Treg", "KO_heat_up", "rank_padj_pooled")
st_tcon  <- gv("Tcon", "sting_specific_up", "padj_pooled")
st_size  <- gv("Tcon", "sting_specific_up", "set_size")
n_treg   <- gv("Treg", "WT_heat_up", "n_tests_pooled")
sig_treg <- gv("Treg", "WT_heat_up", "n_sig_pooled")
n_shared <- named$n_genes_shared_with_WT_heat_up[named$pathway_id == "KO_heat_up"][1]

sting_treg <- stats |> dplyr::filter(population == "Treg", thread == "cGAS-STING")
sting_best <- sting_treg |> dplyr::arrange(padj_pooled) |> dplyr::slice_head(n = 1)
big_band   <- "130 to 150 genes"; small_band <- "10 to 22 genes"

save_overview(
  p, STAGE, STEM, table = tbl,
  finding = sprintf(paste0("Scored with no favourite against all %s sets tested in the JIA Treg ",
                           "contrast, HALLMARK_HYPOXIA reaches NES %+.2f at pooled FDR %s, rank ",
                           "%s of %s by pooled FDR, while the best-placed of the six cGAS-STING ",
                           "sets reaches pooled FDR %s at rank %s. sting_specific_up ",
                           "does reach pooled FDR %s in Tcon on %d genes, so the cGAS-STING ",
                           "reading turns on which population is read. Set size tracks the ",
                           "outcome closely: %.1f%% of Treg sets of %s reach pooled significance ",
                           "against %.1f%% of sets of %s, the band five of the six cGAS-STING ",
                           "sets fall in, and %s of %s Treg tests are significant at all."),
                    cma(n_treg), hyp_nes, fmt_p(hyp_p), cma(hyp_rank), cma(n_treg),
                    fmt_p(sting_best$padj_pooled[1]), cma(min(sting_treg$rank_padj_pooled)),
                    fmt_p(st_tcon), st_size,
                    100 * bv("Treg", big_band, "frac_pooled_significant"), big_band,
                    100 * bv("Treg", small_band, "frac_pooled_significant"), small_band,
                    cma(sig_treg), cma(n_treg)),
  script = SCRIPT, fn = "main",
  config_kv = sprintf("gsea_min_size=%s; gsea_max_size=%s; gsea_fdr=%s; nes_cap=%s; padj_pooled_method=%s; set_selection=tables/sweep_named_sets.csv",
                      FIG_CFG$thresholds$gsea_min_size, FIG_CFG$thresholds$gsea_max_size,
                      FDR, NESCAP,
                      FIG_CFG$unbiased_enrichment$padj_pooled_method %||% "BH"),
  input = "03_results/14_unbiased_enrichment/tables/{gsea_all,sweep_named_sets,sweep_named_sets_stats,sweep_setsize_baseline}.csv",
  ## Every one of the four bounds below is load-bearing, so the block is long by the
  ## README's usual standard. It is tightened as far as it goes without dropping a number.
  how_to_read = sprintf(paste(
    "Columns are the three sorted populations on one shared row axis. Each of the eleven",
    "upper rows is one named set, coloured by comparison thread and ordered inside a",
    "thread by its Treg NES. Below the dashed separator the bottom row is every set",
    "tested in that population, one grey point each, which is the distribution a marker",
    "is read against. Horizontal position is NES clamped to plus or minus %.1f; right of",
    "zero the set's genes concentrate on the synovial-fluid side of the ranking. A filled",
    "marker reaches pooled FDR < %.2g, an open marker sits above it. Marker area is genes",
    "reaching the ranked list, so a large NES on a small marker rests on few genes. Grey",
    "text gives pooled FDR and rank within that population's whole sweep, ranked on",
    "pooled FDR alone; by NES the same sets order differently, and the two orderings",
    "answer different questions. A cell reading 'not tested' had fewer than the minimum",
    "five of its genes in that ranked list, so it carries no result. Four bounds on the",
    "reading. This contrast moves many programs at once: %s of %s tests reach pooled FDR",
    "< %.2g in Treg, %s of %s in Tcon and %s of %s in CD8. Set size drives that rate: in",
    "Treg %s of %s sets of %s are pooled-significant (%.1f%%) against %s of %s sets of %s",
    "(%.1f%%); five of the six cGAS-STING sets sit in that smaller band and the sixth at",
    "five genes, while HALLMARK_HYPOXIA carries %d testable genes. KO_heat_up is drawn",
    "beside WT_heat_up because that comparator reaches pooled FDR %s against %s, ranks %s",
    "against %s, and the two lists share %d genes. And the cGAS-STING family's own signs",
    "disagree: two regulation-of terms carry opposite sign and the positive-regulation",
    "term runs negative. The selection, its reason per row, and the two excluded",
    "substring matches are committed in tables/sweep_named_sets.csv. Correlative",
    "throughout."),
    NESCAP, FDR,
    cma(sig_treg), cma(n_treg), FDR,
    cma(gv("Tcon", "WT_heat_up", "n_sig_pooled")), cma(gv("Tcon", "WT_heat_up", "n_tests_pooled")),
    cma(gv("CD8", "WT_heat_up", "n_sig_pooled")), cma(gv("CD8", "WT_heat_up", "n_tests_pooled")),
    cma(bv("Treg", big_band, "n_pooled_significant")), cma(bv("Treg", big_band, "n_sets_tested")),
    big_band, 100 * bv("Treg", big_band, "frac_pooled_significant"),
    cma(bv("Treg", small_band, "n_pooled_significant")), cma(bv("Treg", small_band, "n_sets_tested")),
    small_band, 100 * bv("Treg", small_band, "frac_pooled_significant"),
    gv("Treg", "HALLMARK_HYPOXIA", "set_size"),
    fmt_p(ko_p), fmt_p(wt_p), cma(ko_rank), cma(wt_rank), n_shared),
  config = FIG_CFG, width = 16.5, height = 11)

message(sprintf("[14_named] wrote %s: %d named rows x %d populations over %s background sets",
                STEM, N_SET, length(POP_LEVELS), cma(nrow(sweep))))
