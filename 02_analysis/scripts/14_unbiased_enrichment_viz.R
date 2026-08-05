#!/usr/bin/env Rscript
# 14_unbiased_enrichment_viz.R — VIZ (no statistics)
# =============================================================================
# Reads the tables 14_unbiased_enrichment.R wrote and renders three figures. It
# computes nothing: every number on every panel is read from a CSV, and each
# figure's source table is written as its same-stem neighbour so a reader can
# check the panel against the numbers behind it.
#
# The three panels answer three different questions, which is why there are three
# and not one:
#   pooled_overview_by_population  How much of the curated universe moves with the
#       synovial-fluid side of this contrast, collection by collection, and where do
#       the mouse-derived arms sit inside that distribution? This is the calibration
#       panel — the one that says whether the mouse-derived enrichment is
#       distinctive or ordinary.
#   treg_top_sets                  In the Treg compartment specifically, WHICH sets
#       are strongest in each direction, named in full.
#   progeny_activity_panel         The same contrast read by a method that uses no
#       gene-set list at all, so the answer does not inherit set-size and curation
#       choices.
#
# READABILITY IS A REQUIREMENT HERE, NOT A NICETY. The equivalent mouse panel
# stacked fourteen sub-panels into one column and clipped its own title, which made
# it unreadable at half width. Every choice below is a reaction to a real failure
# seen in a rendered draft of this stage:
#   * populations go SIDE BY SIDE on one shared row axis, never stacked;
#   * MSigDB identifiers carry no spaces, so str_wrap alone cannot break them —
#     underscores are turned into spaces BEFORE wrapping. Nothing is truncated;
#   * a count that differs per population is annotated INSIDE each panel, because a
#     shared y-axis label can only carry a number summed over panels, and that
#     number would be wrong for every panel it sits beside;
#   * the three mouse-derived arms get one row EACH rather than three labels repelled
#     off one row, which is what collided and overran the panel edge in the draft;
#   * titles, subtitles and captions are hard-wrapped, because a long single-line
#     subtitle is silently clipped at the canvas edge.
#
# Input  (03_results/14_unbiased_enrichment/tables/):
#   gsea_all.csv, gsea_pooled_summary_by_db.csv, progeny_activity.csv,
#   progeny_sf_vs_pb.csv
#
# Output (03_results/14_unbiased_enrichment/):
#   figures/_overview/{pooled_overview_by_population,treg_top_sets,
#                      progeny_activity_panel}.{pdf,png}
#   tables/_overview/<same stem>.csv        source table for each figure
#   README.md captions                     (via save_overview)
#
# Run from the compartment root, AFTER 14_unbiased_enrichment.R:
#   Rscript 02_analysis/scripts/14_unbiased_enrichment_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview

STAGE  <- "14_unbiased_enrichment"
SCRIPT <- "02_analysis/scripts/14_unbiased_enrichment_viz.R"
TDIR   <- file.path(FIG_CFG$paths$results %||% "03_results/", STAGE,
                    FIG_CFG$paths$stage_tables_subdir %||% "tables")

FDR    <- as.numeric(FIG_CFG$thresholds$gsea_fdr %||% 0.05)
TOP_N  <- as.integer(FIG_CFG$figures$top_n %||% 20L)
NESCAP <- as.numeric(FIG_CFG$figures$nes_cap %||% 3.5)
DIV    <- FIG_CFG$colors$diverging
OI     <- FIG_CFG$colors$okabe_ito
THEME  <- project_theme(config = FIG_CFG)

POP_LEVELS <- c("Treg", "Tcon", "CD8")
## The one population palette, read from analysis_config.yaml::colors.populations, so a
## population keeps one colour across every figure this compartment ships.
POP_COL    <- population_colors(FIG_CFG)
GATE_DB    <- FIG_CFG$unbiased_enrichment$mouse_projection$name %||% "mouse_projection"

read_tbl <- function(f) {
  p <- file.path(TDIR, f)
  if (!file.exists(p))
    stop("[14_viz] ", p, " not found. Run 02_analysis/scripts/14_unbiased_enrichment.R first.")
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)
}

sweep   <- read_tbl("gsea_all.csv") |>
  dplyr::mutate(population = factor(population, levels = POP_LEVELS))
by_db   <- read_tbl("gsea_pooled_summary_by_db.csv")
progeny <- read_tbl("progeny_activity.csv")
pg_pair <- read_tbl("progeny_sf_vs_pb.csv")

#' Render a p-value for an on-face label: fixed above 0.001, else one-digit scientific.
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.0e", p)))

#' Hard-wrap a title / subtitle / caption so it cannot be clipped at the canvas edge.
wrap_at <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

#' Make an MSigDB-style identifier legible on a categorical axis WITHOUT truncating it.
#'
#' str_wrap breaks on whitespace, and these identifiers contain none — a raw
#' str_wrap leaves a 110-character single line that pushes the panel off the canvas.
#' Underscores become spaces first, which is display prettifying, so every token
#' survives and the full identifier stays in the source table and in the tooltip-grade
#' `pathway_id` column.
pretty_set <- function(id, width = 38) stringr::str_wrap(gsub("_", " ", id), width = width)

# Own this stage's figure namespace so a rename never leaves an orphan behind.
if (exists("purge_figures"))
  for (stem in c("pooled_overview_by_population", "treg_top_sets", "progeny_activity_panel"))
    purge_figures(STAGE, stem, overview = TRUE, config = FIG_CFG)

# =============================================================================
# FIGURE 1 — pooled_overview_by_population
# =============================================================================
## One row per collection, one COLUMN per sorted population, so the eye compares
## populations horizontally along a shared row axis instead of scrolling a
## fourteen-panel stack. Each collection row is the NES distribution of its
## pooled-significant sets. The three mouse-derived arms are lifted out onto rows of
## their OWN, one arm each, so a reader can see at a glance whether an arm sits at
## the edge of the curated distribution or inside it — and so their labels cannot
## collide, which is what wrecked the first draft.

sig <- sweep |> dplyr::filter(padj_pooled < FDR)
ARMS <- sweep |> dplyr::filter(database == GATE_DB) |> dplyr::pull(pathway_id) |> unique()

db_order <- sig |> dplyr::filter(database != GATE_DB) |>
  dplyr::count(database, name = "n") |> dplyr::arrange(n) |> dplyr::pull(database)
db_order <- c(setdiff(unique(sweep$database[sweep$database != GATE_DB]), db_order), db_order)
ROW_LEVELS <- c(db_order, rev(sort(ARMS)))   # arms pinned to the top of the panel

strip <- sig |>
  dplyr::filter(database != GATE_DB) |>
  dplyr::mutate(row = factor(database, levels = ROW_LEVELS),
                nes_plot = pmax(pmin(nes, NESCAP), -NESCAP))
arms <- sweep |>
  dplyr::filter(database == GATE_DB) |>
  dplyr::mutate(row = factor(pathway_id, levels = ROW_LEVELS),
                nes_plot = pmax(pmin(nes, NESCAP), -NESCAP),
                survives = padj_pooled < FDR)

## The count that matters differs per population, so it is annotated INSIDE each
## panel rather than on the shared axis label. Collection rows carry "significant of
## tested"; the arm rows carry their own NES and pooled FDR, which is the number a
## reader wants for a single set.
ann_db <- by_db |>
  dplyr::filter(database != GATE_DB) |>
  dplyr::transmute(population = factor(population, levels = POP_LEVELS),
                   row = factor(database, levels = ROW_LEVELS),
                   label = sprintf("%d of %d", sig_pooled, n_tests_in_db))
## Two lines, not one: a single "NES +1.47, FDR 0.235" ran past the annotation gutter
## and lost its own FDR value to the canvas edge in the first render.
ann_arm <- arms |>
  dplyr::transmute(population, row,
                   label = sprintf("NES %+.2f\nFDR %s", nes, fmt_p(padj_pooled)))
ann <- dplyr::bind_rows(ann_db, ann_arm) |> dplyr::mutate(x = NESCAP * 1.12)

p1 <- ggplot(strip, aes(x = nes_plot, y = row)) +
  geom_hline(yintercept = length(db_order) + 0.5, linewidth = 0.5,
             linetype = "22", colour = "grey45") +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_jitter(aes(colour = direction), height = 0.26, width = 0, alpha = 0.32, size = 1.0) +
  geom_point(data = arms, aes(shape = survives), size = 4.0, stroke = 1.1,
             fill = OI$yellow, colour = "black") +
  geom_text(data = ann, aes(x = x, y = row, label = label), inherit.aes = FALSE,
            hjust = 0, size = 3.1, lineheight = 0.9, colour = "grey25") +
  scale_colour_manual(values = c(up = DIV$up, down = DIV$down),
                      labels = c(up = "enriched toward synovial fluid",
                                 down = "enriched toward paired blood"),
                      name = NULL, drop = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 5),
                     labels = c(`TRUE` = sprintf("mouse-derived up arm, FDR < %.2g", FDR),
                                `FALSE` = "mouse-derived up arm, not significant"),
                     name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-NESCAP * 1.02, NESCAP * 1.62),
                     breaks = seq(-3, 3, by = 1.5)) +
  facet_wrap(~ population, nrow = 1) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.6)),
         shape  = guide_legend(override.aes = list(size = 3.4))) +
  labs(
    title = "What the synovial-fluid-versus-blood contrast contains, one collection at a time",
    subtitle = wrap_at(sprintf(paste0("Every set in eleven collections scored on the same frozen ",
                                      "ranked lists, with Benjamini-Hochberg pooled across all ",
                                      "%s to %s tests asked of one population. Small points are ",
                                      "the sets reaching FDR < %.2g; the text on each row is how ",
                                      "many of that collection's sets did, out of how many were ",
                                      "tested in that population."),
                              format(min(sweep$n_tests_pooled), big.mark = ","),
                              format(max(sweep$n_tests_pooled), big.mark = ","), FDR), 128),
    x = sprintf("normalized enrichment score, clamped to ±%.1f", NESCAP),
    y = NULL,
    caption = wrap_at(paste("Above the dashed line each row is a single mouse-derived up arm.",
                            "Correlative: a set enriching says its gene content moves with the",
                            "synovial-fluid side of this contrast, not that the program it is",
                            "named for is present."), 120)) +
  THEME +
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10.5))

tbl1 <- sweep |>
  dplyr::filter(padj_pooled < FDR | database == GATE_DB) |>
  dplyr::left_join(by_db |> dplyr::select(population, database,
                                          n_sig_pooled_in_database = sig_pooled,
                                          n_sig_per_database_only = sig_per_database),
                   by = c("population", "database")) |>
  dplyr::select(population, database, pathway_id, direction, nes, padj, padj_pooled,
                set_size, n_tests_in_db, n_tests_pooled,
                n_sig_pooled_in_database, n_sig_per_database_only) |>
  dplyr::arrange(population, database, padj_pooled)

n_sig_treg  <- sum(sig$population == "Treg")
wt_treg     <- sweep |> dplyr::filter(population == "Treg", pathway_id == "WT_heat_up")
bigger_treg <- sig |> dplyr::filter(population == "Treg",
                                    abs(nes) > abs(wt_treg$nes[1])) |> nrow()
## The count of larger-|NES| sets needs its DIRECTION SPLIT to be readable. Ranking on
## absolute NES is the fair comparison — a signed ranking would flatter the arm by
## ignoring everything moving the other way — but almost all the sets that beat the arm
## on magnitude are on the blood side, so the bare count reads as "unremarkable" when
## conditioned on direction the arm is near the top. Both numbers go on the caption.
bigger_dn_treg <- sig |> dplyr::filter(population == "Treg", direction == "down",
                                       abs(nes) > abs(wt_treg$nes[1])) |> nrow()
sf_side_treg   <- sig |> dplyr::filter(population == "Treg", direction == "up")
wt_rank_sf     <- which(sf_side_treg$pathway_id[order(-abs(sf_side_treg$nes))] == "WT_heat_up")[1]

save_overview(
  p1, STAGE, "pooled_overview_by_population", table = tbl1,
  finding = sprintf(paste0("The synovial-fluid-versus-paired-blood contrast moves a great many ",
                           "curated programs at once — %s of %s tests reach FDR < %.2g in Treg ",
                           "after pooling — and while %d pooled-significant sets carry a larger ",
                           "absolute NES than the mouse-derived WT_heat_up arm (NES %+.2f, pooled ",
                           "FDR %s), %d of those %d are on the blood side, so among the %s sets ",
                           "moving toward synovial fluid the arm ranks %d."),
                    format(n_sig_treg, big.mark = ","),
                    format(sweep$n_tests_pooled[sweep$population == "Treg"][1], big.mark = ","),
                    FDR, bigger_treg, wt_treg$nes[1], fmt_p(wt_treg$padj_pooled[1]),
                    bigger_dn_treg, bigger_treg,
                    format(nrow(sf_side_treg), big.mark = ","), wt_rank_sf),
  script = SCRIPT, fn = "main",
  config_kv = sprintf("gsea_min_size=%s; gsea_max_size=%s; gsea_fdr=%s; nes_cap=%s; padj_pooled_method=%s",
                      FIG_CFG$thresholds$gsea_min_size, FIG_CFG$thresholds$gsea_max_size,
                      FDR, NESCAP,
                      FIG_CFG$unbiased_enrichment$padj_pooled_method %||% "BH"),
  input = "03_results/14_unbiased_enrichment/tables/{gsea_all,gsea_pooled_summary_by_db}.csv",
  how_to_read = paste(
    "Columns are the three sorted populations side by side. Below the dashed line each",
    "row is a collection, ordered by how many of its sets reach significance; above it",
    "each row is one mouse-derived up arm. A small point is a set at FDR < 0.05 after",
    "pooling across every test asked of that ranking. Brown concentrates on the",
    "synovial-fluid side, blue on paired blood. Horizontal position is the exact NES,",
    "clamped to plus or minus 3.5. Yellow diamonds are the arms, filled when",
    "significant. Grey text is that population's own count, or an arm's NES and FDR.",
    "Read for calibration: far right in a dense row is strong and ordinary at once.",
    "Correlative."),
  config = FIG_CFG, width = 14, height = 9)

# =============================================================================
# FIGURE 2 — treg_top_sets
# =============================================================================
## Treg is the compartment this project exists for, so it gets its sets NAMED. The
## cap is split evenly between the two directions: taking the top 20 by absolute NES
## returned twenty translation and ribosome sets and hid the synovial-fluid side
## entirely, which is a presentation artefact of one arm carrying larger magnitudes,
## not a finding. A lollipop rather than a bar because the readable quantity is a
## position on the NES axis, not an area.

treg_sig <- sweep |> dplyr::filter(population == "Treg", padj_pooled < FDR)
half <- max(as.integer(TOP_N / 2), 1L)
top2 <- treg_sig |>
  dplyr::group_by(direction) |>
  dplyr::arrange(dplyr::desc(abs(nes)), padj_pooled, .by_group = TRUE) |>
  dplyr::slice_head(n = half) |>
  dplyr::ungroup() |>
  dplyr::arrange(nes) |>
  dplyr::mutate(is_mouse = database == GATE_DB,
                label = pretty_set(sprintf("%s [%s]", pathway_id, database), width = 44),
                label = factor(label, levels = label))

xr <- range(top2$nes)
x_txt <- xr[2] + diff(xr) * 0.06

p2 <- ggplot(top2, aes(x = nes, y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_segment(aes(x = 0, xend = nes, yend = label, colour = direction),
               linewidth = 1.5, show.legend = FALSE) +
  geom_point(aes(colour = direction, size = set_size)) +
  geom_point(data = dplyr::filter(top2, is_mouse), shape = 21, size = 6.6,
             stroke = 1.4, colour = "black", fill = NA) +
  geom_text(aes(x = x_txt, label = fmt_p(padj_pooled)), hjust = 0, size = 3.3,
            colour = "grey25") +
  scale_colour_manual(values = c(up = DIV$up, down = DIV$down),
                      labels = c(up = "toward synovial fluid", down = "toward paired blood"),
                      name = NULL) +
  scale_size_continuous(range = c(2.2, 5.6), name = "genes in the ranked list") +
  scale_x_continuous(limits = c(xr[1] - diff(xr) * 0.04, x_txt + diff(xr) * 0.16)) +
  labs(
    title = "The strongest programs in the sorted JIA Treg niche contrast",
    subtitle = wrap_at(sprintf(paste0("Top %d in each direction, by absolute NES, among the %s ",
                                      "sets reaching FDR < %.2g after pooling across the whole ",
                                      "sweep. %s further significant sets are not shown, so the ",
                                      "cap must not be read as completeness. The two directions ",
                                      "are capped separately because the paired-blood arm carries ",
                                      "larger magnitudes and would otherwise fill every row."),
                              half, format(nrow(treg_sig), big.mark = ","), FDR,
                              format(max(nrow(treg_sig) - nrow(top2), 0), big.mark = ",")), 112),
    x = "normalized enrichment score, synovial fluid versus paired blood",
    y = NULL,
    caption = wrap_at(paste("The grey number beside each point is the pooled FDR.",
                            "A black ring marks a mouse-derived up arm. Underscores in",
                            "identifiers are shown as spaces; nothing is truncated.",
                            "Correlative."), 108)) +
  THEME +
  # Legends STACKED, not side by side: laid out horizontally the two guides span the
  # full canvas and the last legend label ends flush against the right edge.
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.y = element_text(size = 9, lineheight = 0.92),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

tbl2 <- top2 |>
  dplyr::select(population, database, pathway_id, pathway_name, direction, nes,
                pvalue, padj, padj_pooled, set_size, leading_edge_size,
                n_tests_in_db, n_tests_pooled) |>
  dplyr::arrange(direction, dplyr::desc(abs(nes)))

best_up <- top2 |> dplyr::filter(direction == "up") |>
  dplyr::arrange(dplyr::desc(nes)) |> dplyr::slice_head(n = 1)
best_dn <- top2 |> dplyr::filter(direction == "down") |>
  dplyr::arrange(nes) |> dplyr::slice_head(n = 1)

save_overview(
  p2, STAGE, "treg_top_sets", table = tbl2,
  finding = sprintf(paste0("In the JIA Treg contrast the largest shifts are downward: %s reaches ",
                           "NES %+.2f toward paired blood, against %s at %+.2f toward synovial ",
                           "fluid, so the niche difference this compartment reads as an ",
                           "inflammatory gain is accompanied by an at least equally large loss of ",
                           "translation and ribosomal programs."),
                    best_dn$pathway_id[1], best_dn$nes[1],
                    best_up$pathway_id[1], best_up$nes[1]),
  script = SCRIPT, fn = "main",
  config_kv = sprintf("figures.top_n=%d (split %d per direction); gsea_fdr=%s; gsea_min_size=%s; gsea_max_size=%s",
                      TOP_N, half, FDR, FIG_CFG$thresholds$gsea_min_size,
                      FIG_CFG$thresholds$gsea_max_size),
  input = "03_results/14_unbiased_enrichment/tables/gsea_all.csv",
  how_to_read = paste(
    "One row per gene set, capped at the top ten in each direction by absolute NES among",
    "sets at FDR < 0.05 after pooling; the subtitle states how many the cap leaves out, so",
    "the panel is a top-N view. Identifiers are shown with underscores as spaces and",
    "wrapped in full, each with its collection in brackets. Right of zero the",
    "set's genes concentrate on the synovial-fluid side of this ranking, left on paired",
    "blood. Point size is how many of the set's genes reach the ranked list, so a large",
    "NES on a small point rests on few genes. The grey number is the pooled FDR; a black",
    "ring marks a mouse-derived arm. Correlative."),
  config = FIG_CFG, width = 12.5, height = 11)

# =============================================================================
# FIGURE 3 — progeny_activity_panel
# =============================================================================
## The gene-set-free read of the same contrast. Fourteen pathways is few enough to
## name every row, so this is a dot panel and not a heatmap: a position on a shared
## score axis is read more accurately than a fill. Two tests are shown on one panel
## because they check each other — the model score on the contrast statistics and a
## donor-paired test of the per-donor activities — and a pathway significant by only
## one of them must not look the same as one significant by both.

pg <- progeny |>
  dplyr::mutate(population = factor(population, levels = POP_LEVELS)) |>
  dplyr::left_join(
    pg_pair |> dplyr::transmute(population, pathway_name,
                                donor_padj = padj, donor_difference = mean_difference,
                                n_paired_donors),
    by = c("population", "pathway_name")) |>
  dplyr::mutate(sig_contrast = padj < FDR,
                sig_donor    = !is.na(donor_padj) & donor_padj < FDR)

ord <- pg |> dplyr::filter(population == "Treg") |> dplyr::arrange(nes) |>
  dplyr::pull(pathway_name)
pg <- pg |> dplyr::mutate(pathway_name = factor(pathway_name, levels = ord))

p3 <- ggplot(pg, aes(x = nes, y = pathway_name)) +
  geom_vline(xintercept = 0, linewidth = 0.6, colour = "grey55") +
  geom_line(aes(group = pathway_name), colour = "grey80", linewidth = 0.7) +
  geom_point(aes(colour = population, alpha = sig_contrast), size = 3.8) +
  geom_point(data = dplyr::filter(pg, sig_donor), shape = 21, size = 6.2,
             stroke = 1.2, colour = "black", fill = NA) +
  scale_colour_manual(values = POP_COL, name = "sorted population") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.32),
                     labels = c(`TRUE` = sprintf("FDR < %.2g", FDR),
                                `FALSE` = "not significant"),
                     name = "contrast-statistic test", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3.4)),
         alpha  = guide_legend(override.aes = list(size = 3.4, colour = "grey20"))) +
  labs(
    title = "Signalling footprints in the same contrast, read without a gene-set list",
    subtitle = wrap_at(paste("decoupleR multivariate linear model on the human PROGENy model -",
                             "fourteen pathways, continuous weights, top 500 target genes each -",
                             "run on the donor-pseudobulk moderated-t contrast statistics.",
                             "Rows are ordered by the Treg score."), 112),
    x = "PROGENy activity score, synovial fluid versus paired blood",
    y = NULL,
    caption = wrap_at(paste("A black ring marks a pathway that is also significant in the",
                            "independent donor-paired test of per-donor activities.",
                            "Correlative: an activity score is a footprint inferred from target-gene",
                            "expression, not a measurement of pathway activity."), 108)) +
  THEME +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey88"),
        panel.grid.minor = element_blank())

tbl3 <- pg |>
  dplyr::select(population, contrast, pathway_name, database,
                contrast_score = nes, contrast_pvalue = pvalue, contrast_padj = padj,
                direction, set_size, donor_difference, donor_padj, n_paired_donors,
                sig_contrast, sig_donor) |>
  dplyr::arrange(population, contrast_padj)

hyp  <- pg |> dplyr::filter(population == "Treg", pathway_name == "Hypoxia")
topA <- pg |> dplyr::filter(population == "Treg") |> dplyr::arrange(dplyr::desc(nes)) |>
  dplyr::slice_head(n = 1)
n_both <- pg |> dplyr::filter(population == "Treg", sig_contrast, sig_donor) |> nrow()

save_overview(
  p3, STAGE, "progeny_activity_panel", table = tbl3,
  finding = sprintf(paste0("Read without any gene-set list, the JIA Treg synovial-fluid-versus-blood ",
                           "contrast carries its largest PROGENy footprint in %s (score %+.2f, FDR ",
                           "%s) while the Hypoxia footprint scores %+.2f at FDR %s, so the ",
                           "inflammatory and low-oxygen readouts rise together in the same niche ",
                           "contrast; which of the two drives the other is untested here. %d of ",
                           "the fourteen pathways are significant on both tests in Treg."),
                    as.character(topA$pathway_name[1]), topA$nes[1], fmt_p(topA$padj[1]),
                    hyp$nes[1], fmt_p(hyp$padj[1]), n_both),
  script = SCRIPT, fn = "main",
  config_kv = sprintf("progeny.organism=%s; progeny.top=%s; progeny.minsize=%s; gsea_fdr=%s",
                      FIG_CFG$unbiased_enrichment$progeny$organism,
                      FIG_CFG$unbiased_enrichment$progeny$top,
                      FIG_CFG$unbiased_enrichment$progeny$minsize, FDR),
  input = "03_results/14_unbiased_enrichment/tables/{progeny_activity,progeny_sf_vs_pb}.csv",
  how_to_read = paste(
    "One row per PROGENy pathway, ordered by its Treg score, with a grey line joining the",
    "three populations so the spread within a row is the between-population difference.",
    "Horizontal position is the model activity score on the donor-pseudobulk moderated-t",
    "contrast statistics: right of zero the pathway's footprint genes move with synovial",
    "fluid, left with paired blood. Colour is the population. A solid point reaches",
    "FDR < 0.05 on that test; a black ring marks one that also reaches it in the",
    "independent donor-paired test, so ringed and solid is corroborated twice while faded",
    "and unringed is corroborated by neither. A footprint is inferred from target-gene",
    "expression; pathway activity itself is untested. Correlative; no causal reading."),
  config = FIG_CFG, width = 12.5, height = 8.5)

message("[14_viz] wrote 3 overview figures with same-stem source tables and README captions")
