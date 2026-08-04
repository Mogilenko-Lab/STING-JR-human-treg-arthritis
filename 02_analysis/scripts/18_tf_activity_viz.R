#!/usr/bin/env Rscript
# 18_tf_activity_viz.R: VIZ. Three figures off the tables 18_tf_activity.R wrote.
# =============================================================================
#   tf_rank_cascade              rank of each focus TF across every network variant
#                                crossed with every estimator, plus the committed
#                                unsigned-regulon fgsea rank
#   tf_target_promiscuity        per-target signed contribution against how many other
#                                CollecTRI regulons claim the same target
#   tf_activity_vs_regulon_size  activity against regulon size across every TF tested,
#                                with the size-conditional expectation and the
#                                size-and-expression-matched random-regulon null
#
# One claim each, and no statistics here. Every figure goes through save_overview(), which
# emits the vector PDF, the raster PNG, the same-stem source CSV and the stage README
# caption in one call.
#
# Reads 03_results/18_tf_activity/tables/:
#   hif1a_rank_cascade.csv, target_decomposition.csv, target_decomposition_summary.csv,
#   regulon_size_calibration.csv, regulon_size_spearman.csv, size_matched_null.csv
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/18_tf_activity_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview

STAGE  <- "18_tf_activity"
SCRIPT <- "02_analysis/scripts/18_tf_activity_viz.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

CFG <- FIG_CFG
TA  <- CFG$tf_activity
FIG <- CFG$figures
OKI <- unlist(CFG$colors$okabe_ito)
DIV <- unlist(CFG$colors$diverging)
REFC <- CFG$colors$reference_line %||% unname(OKI["black"])

FOCUS   <- unlist(TA$focus_tfs)
DECOMP  <- unlist(TA$decompose_tfs)
SEL_MAX <- as.integer(TA$selective_max_regulons %||% 1L)
LBL_N   <- as.integer(FIG$volcano_label_top %||% 10L)
LBL_SZ  <- as.numeric(FIG$label_size %||% 4)
PT_SZ   <- as.numeric(FIG$point_size %||% 2.4)
LN_W    <- as.numeric(FIG$line_width %||% 1.0)
# ggrepel searches for label positions at random, so every repel call here takes this seed
# and a re-render reproduces byte for byte.
RPL_SEED <- as.integer(CFG$thresholds$gsea_seed %||% 123L)

TBL <- file.path(CFG$paths$results %||% "03_results/", STAGE,
                 CFG$paths$stage_tables_subdir %||% "tables")
rd <- function(f) readr::read_csv(file.path(TBL, f), show_col_types = FALSE)

CASCADE <- rd("hif1a_rank_cascade.csv")
DEC     <- rd("target_decomposition.csv")
DECSUM  <- rd("target_decomposition_summary.csv")
CAL     <- rd("regulon_size_calibration.csv")
SPEAR   <- rd("regulon_size_spearman.csv")
NULLS   <- rd("size_matched_null.csv")

PRIMARY_LABEL <- CAL$population[1]

# One colour and one point shape per focus TF, reused across figures so a factor keeps its
# identity between panels, and named from the config palette. `okabe_ito$yellow` is left
# out because it washes into a white panel at these point and text sizes. Every series is
# directly labelled too, so shape carries the disambiguation where two warm hues land close.
TF_PALETTE_KEYS <- c(HIF1A = "vermillion", NFKB1 = "blue", STAT3 = "bluish_green",
                     CREB1 = "reddish_purple", ATF3 = "sky_blue", HSF1 = "orange",
                     EPAS1 = "black")
tf_colors <- c(vapply(TF_PALETTE_KEYS, function(k) unname(OKI[[k]]), character(1)),
               REL = unname(DIV[["up"]]))
tf_colors <- tf_colors[FOCUS]
tf_shapes <- setNames(c(16, 17, 15, 18, 8, 4, 6, 3)[seq_along(FOCUS)], FOCUS)

# =============================================================================
# FIGURE 1: the rank cascade
# =============================================================================

message("[fig 1] rank cascade ...")

casc <- CASCADE %>%
  filter(population == PRIMARY_LABEL) %>%
  mutate(variant = factor(variant, levels = c("unsigned_geneset", "signed", "unsigned",
                                              "literature_signed", "alias_recovered")),
         method  = factor(method, levels = c("fgsea", "ulm", "mlm", "consensus"))) %>%
  arrange(variant, method) %>%
  mutate(configuration = factor(configuration, levels = unique(configuration)),
         tf = factor(tf, levels = FOCUS))

casc_end <- casc %>% group_by(tf) %>%
  filter(as.integer(configuration) == max(as.integer(configuration))) %>% ungroup()

p1 <- ggplot(casc, aes(x = configuration, y = rank, colour = tf, group = tf)) +
  geom_line(linewidth = LN_W) +
  geom_point(aes(shape = tf), size = PT_SZ, stroke = 0.9) +
  geom_text_repel(data = casc_end, aes(label = tf), size = LBL_SZ,
                  nudge_x = 0.35, hjust = 0, direction = "y", segment.size = 0.3,
                  min.segment.length = 0, seed = RPL_SEED, show.legend = FALSE) +
  # Log scale keeps the top of the ranking legible, reversed so rank 1 sits at the top.
  scale_y_continuous(transform = scales::compose_trans("log10", "reverse"),
                     breaks = c(1, 3, 10, 30, 100, 300),
                     expand = expansion(mult = c(0.06, 0.06))) +
  scale_x_discrete(expand = expansion(add = c(0.5, 2.8))) +
  scale_colour_manual(values = tf_colors, name = "Factor") +
  scale_shape_manual(values = tf_shapes, guide = "none") +
  labs(title = "TF rank across network variant and estimator",
       subtitle = "Rank by descending activity within a configuration; rank 1 is the most activated",
       x = "Network variant / estimator", y = "Rank (log scale, inverted)") +
  project_theme(config = CFG) +
  # Rotated tick labels need left margin, or the canvas clips the first one.
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "none",
        plot.margin = margin(8, 12, 8, 58))

fig1_tbl <- casc %>%
  select(population, tf, variant, method, configuration, regulon_size, score, padj,
         rank, n_tfs_scored, pct_rank) %>%
  arrange(tf, configuration)

save_overview(
  p1, STAGE, "tf_rank_cascade", table = fig1_tbl,
  finding = paste(
    "HIF1A's inferred-activity rank on the sorted-Treg synovial-fluid-versus-paired-blood",
    "contrast sits in the top twelve in twelve of the thirteen network-by-estimator",
    "configurations, its one remaining placement being rank 42 of 388 under the",
    "literature-signed network scored multivariately, which makes it the steadiest of the eight.",
    "The same axes move its neighbours much further: NFKB1 from rank 3 to rank 138 under signed",
    "MLM and rank 252 under unsigned MLM, REL from rank 2 to rank 298 under literature-signed",
    "MLM, so the rank instability the mouse anchor documented for its HIF1a result falls here on",
    "the NF-kB and AP-1 members."),
  script = SCRIPT, fn = "save_overview",
  config_kv = paste0("tf_activity.network_variants=[signed, unsigned, literature_signed, ",
                     "alias_recovered]; tf_activity.methods=[ulm, mlm, consensus]; ",
                     "thresholds.gsea_min_size=", CFG$thresholds$gsea_min_size),
  input = "03_results/18_tf_activity/tables/hif1a_rank_cascade.csv",
  how_to_read = paste(
    "One line per factor, labelled at the right edge and carrying its own point shape so two",
    "lines of similar hue stay separable. The y axis is rank by descending activity within a",
    "configuration, inverted on a log scale, so higher means more activated and rank 1 is at the",
    "top. ULM scores each regulon on its own; MLM fits every regulon jointly, so a factor whose",
    "targets are shared with other regulons loses rank there. `signed` uses CollecTRI's recorded",
    "per-edge mode of regulation, `unsigned` forces every edge positive, `literature_signed`",
    "keeps only evidence-signed edges, `alias_recovered` adds targets resolved to the pre-2019",
    "symbol this matrix carries. The leftmost column is the committed unsigned-regulon fgsea",
    "rank. Denominators differ between configurations and sit in the source table.",
    "Annotation tier: an inferred activity is a statistic over target-gene expression, and",
    "nothing here pools with the donor-pseudobulk claim spine."),
  config = CFG, wide = TRUE
)

# =============================================================================
# FIGURE 2: target contribution against target promiscuity
# =============================================================================

message("[fig 2] target promiscuity ...")

dec <- DEC %>%
  filter(tf %in% DECOMP) %>%
  mutate(tf = factor(tf, levels = DECOMP),
         claimed = ifelse(n_regulons <= SEL_MAX,
                          "claimed by this regulon alone", "claimed by other regulons too"))

lab2 <- dec %>% group_by(tf) %>% slice_max(contrib, n = LBL_N, with_ties = FALSE) %>% ungroup()

share <- DECSUM %>%
  group_by(tf) %>%
  summarise(pct_selective = sum(pct_of_total_contrib[promiscuity_band == "<=1"]),
            pct_over_25 = sum(pct_of_total_contrib[promiscuity_band == ">25"]), .groups = "drop") %>%
  mutate(tf = factor(tf, levels = DECOMP),
         txt = sprintf("%.2f%% of the signed total from targets this regulon alone claims\n%.0f%% from targets in more than 25 regulons",
                       pct_selective, pct_over_25))

p2 <- ggplot(dec, aes(x = n_regulons, y = contrib)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = REFC) +
  geom_point(aes(colour = claimed), size = PT_SZ, alpha = 0.85) +
  geom_text_repel(data = lab2, aes(label = target), size = LBL_SZ,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                  box.padding = 0.35, seed = RPL_SEED, show.legend = FALSE) +
  geom_text(data = share, aes(x = 1, y = -Inf, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = -0.3, size = LBL_SZ * 0.85, colour = "grey25") +
  facet_wrap(~ tf, nrow = 1) +
  scale_x_continuous(trans = "log10", breaks = c(1, 2, 5, 10, 25, 50, 100, 200)) +
  scale_colour_manual(values = c("claimed by this regulon alone" = unname(DIV["up"]),
                                 "claimed by other regulons too" = unname(OKI["sky_blue"])),
                      name = NULL) +
  labs(title = "Target contribution against target promiscuity",
       subtitle = "Signed contribution = sign(mode of regulation) x moderated t; positive is synovial-fluid-side",
       x = "CollecTRI regulons containing the target (log scale)",
       y = "Signed contribution") +
  project_theme(config = CFG) +
  theme(legend.position = "bottom")

fig2_tbl <- dec %>%
  select(population, tf, target, mor, sign_decision, stat, contrib, n_regulons,
         n_other_regulons, selective, promiscuity_band, avg_expr, padj_gene) %>%
  arrange(tf, desc(contrib))

save_overview(
  p2, STAGE, "tf_target_promiscuity", table = fig2_tbl,
  finding = paste(
    "The synovial-fluid-side contribution to HIF1A's CollecTRI-ULM score is carried by targets",
    "that many other regulons also contain: the 27 of 293 targets HIF1A alone claims sum to 0.14%",
    "of its signed total, while the 73 targets sitting in more than 25 regulons carry 35%.",
    "NFKB1 decomposes the same way (2 exclusive targets, 0.07% of its signed total; 30% from",
    "targets in more than 25 regulons), so joint ownership of the high-t genes is a property",
    "the two regulons share and it bounds both of them equally."),
  script = SCRIPT, fn = "save_overview",
  config_kv = paste0("tf_activity.decompose_tfs=[", paste(DECOMP, collapse = ", "),
                     "]; tf_activity.selective_max_regulons=", SEL_MAX,
                     "; tf_activity.primary_population=", TA$primary_population,
                     "; figures.volcano_label_top=", LBL_N),
  input = "03_results/18_tf_activity/tables/target_decomposition.csv",
  how_to_read = paste(
    "One point per regulon target, faceted by factor. The x axis counts how many CollecTRI",
    "regulons contain that target, on a log scale, so points to the right are jointly owned and",
    "points at x = 1 belong to this regulon alone. The y axis is the target's signed",
    "contribution, its moderated t multiplied by the sign of the edge, so positive means the",
    "target moves with the synovial-fluid side and the zero rule separates the directions. Orange",
    "marks the exclusively-claimed targets. The ten largest positive contributors per facet are",
    "named. The in-panel text gives each factor's share of its signed total from",
    "exclusively-claimed targets and from targets in more than 25 regulons. Annotation tier: a",
    "contribution is arithmetic on the committed ranked list and carries no separate test."),
  config = CFG, wide = TRUE
)

# =============================================================================
# FIGURE 3: activity against regulon size, with the size-matched null
# =============================================================================

message("[fig 3] regulon-size calibration ...")

STAT_LABELS <- c(collectri_ulm_score = "CollecTRI-ULM activity score",
                 unsigned_geneset_fgsea_nes = "Unsigned-regulon fgsea NES")

cal_long <- bind_rows(
  CAL %>% transmute(tf, focus_tf, regulon_size,
                    statistic = "collectri_ulm_score",
                    score = ulm_score, size_expected = ulm_score_size_expected),
  CAL %>% transmute(tf, focus_tf, regulon_size = fgsea_set_size,
                    statistic = "unsigned_geneset_fgsea_nes",
                    score = fgsea_nes, size_expected = fgsea_nes_size_expected)
) %>%
  filter(is.finite(regulon_size), is.finite(score)) %>%
  mutate(statistic = factor(STAT_LABELS[statistic], levels = unname(STAT_LABELS)))

rho <- SPEAR %>% filter(population == PRIMARY_LABEL)
rho_perm <- SPEAR %>% filter(grepl("permuted", population))
rho_lab <- tibble(
  statistic = factor(unname(STAT_LABELS), levels = unname(STAT_LABELS)),
  label = c(sprintf("Spearman rho = %.2f over %d factors\n%.2f on the same contrast with gene labels permuted",
                    rho$spearman_ulm_score_vs_size, rho$n_tfs_ulm,
                    rho_perm$spearman_ulm_score_vs_size),
            sprintf("Spearman rho = %.2f over %d sets",
                    rho$spearman_fgsea_nes_vs_size, rho$n_tfs_fgsea)))

nulls <- NULLS %>%
  filter(null_match == "size_and_expression") %>%
  mutate(statistic = factor(STAT_LABELS[statistic], levels = unname(STAT_LABELS))) %>%
  filter(!is.na(statistic), is.finite(regulon_size))

focus_pts <- cal_long %>% filter(focus_tf)

p3 <- ggplot(cal_long, aes(x = regulon_size, y = score)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = REFC) +
  geom_point(colour = "grey72", size = PT_SZ * 0.55, alpha = 0.7) +
  # Dashed and neutral, so the fitted expectation reads as a reference rather than a series.
  geom_line(aes(y = size_expected), linewidth = LN_W, colour = REFC, linetype = "22") +
  geom_segment(data = nulls, aes(x = regulon_size, xend = regulon_size,
                                 y = null_q95, yend = obs),
               inherit.aes = FALSE, linewidth = 0.5, colour = REFC) +
  geom_point(data = nulls, aes(x = regulon_size, y = null_q95), inherit.aes = FALSE,
             shape = 2, size = PT_SZ, stroke = 0.8, colour = REFC) +
  geom_point(data = focus_pts, aes(colour = tf, shape = tf), size = PT_SZ * 1.6,
             stroke = 1.1) +
  geom_text_repel(data = focus_pts, aes(label = tf, colour = tf), size = LBL_SZ,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                  box.padding = 0.7, point.padding = 0.4, force = 8,
                  seed = RPL_SEED, show.legend = FALSE) +
  scale_shape_manual(values = tf_shapes, guide = "none") +
  geom_text(data = rho_lab, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE,
            hjust = 1.02, vjust = -0.25, size = LBL_SZ * 0.85, colour = "grey25") +
  facet_wrap(~ statistic, nrow = 1, scales = "free_y") +
  scale_x_continuous(trans = "log10", breaks = c(5, 10, 25, 50, 100, 250, 500, 1000)) +
  scale_colour_manual(values = tf_colors, guide = "none") +
  labs(title = "Activity against regulon size on the niche contrast",
       subtitle = "Positive is synovial-fluid-side; both statistics computed on the same ranked list",
       x = "Targets tested in the regulon (log scale)", y = "Activity") +
  project_theme(config = CFG)

fig3_tbl <- cal_long %>%
  left_join(nulls %>% select(tf, statistic, null_mean, null_sd, null_q95, obs,
                             pct_of_null, p_empirical, z_vs_null),
            by = c("tf", "statistic")) %>%
  arrange(statistic, desc(score))

save_overview(
  p3, STAGE, "tf_activity_vs_regulon_size", table = fig3_tbl,
  finding = paste(
    "Inferred activity on this contrast rises with regulon size across every factor tested:",
    sprintf("Spearman rho = %.2f between size and CollecTRI-ULM score over %d factors and %.2f between size and unsigned-regulon fgsea NES over %d sets,",
            rho$spearman_ulm_score_vs_size, rho$n_tfs_ulm,
            rho$spearman_fgsea_nes_vs_size, rho$n_tfs_fgsea),
    sprintf("falling to %.2f when the gene labels of the same ranked list are permuted,",
            rho_perm$spearman_ulm_score_vs_size),
    "which places the size dependence in the breadth of the synovial-fluid-side shift that a",
    "bigger regulon samples more thoroughly, and every large-regulon factor in the headline",
    "table sits on that gradient."),
  script = SCRIPT, fn = "save_overview",
  config_kv = paste0("tf_activity.null_draws=", TA$null_draws,
                     "; tf_activity.null_expression_deciles=", TA$null_expression_deciles,
                     "; thresholds.gsea_min_size=", CFG$thresholds$gsea_min_size,
                     "; thresholds.gsea_max_size=", CFG$thresholds$gsea_max_size),
  input = paste("03_results/18_tf_activity/tables/regulon_size_calibration.csv +",
                "regulon_size_spearman.csv + size_matched_null.csv"),
  how_to_read = paste(
    "One grey point per factor, faceted by statistic, x = the factor's targets present in the",
    "ranked list on a log scale. The dashed dark-grey curve is the size-conditional expectation",
    "fitted over the real regulons themselves, so a point above it is more active than its size",
    "alone accounts for. Coloured labelled points are the headline-table factors, each with its",
    "own shape. The open triangle below each is the 95th percentile of random regulons matched to",
    "that factor on size, on repressing-edge fraction and on average-expression decile",
    "composition, and the stalk spans that percentile to the observed value, so a short stalk",
    "means a factor barely beats a matched bag of genes. In-panel text gives the",
    "size-versus-activity Spearman correlation, and for the left facet the same correlation after",
    "the ranked list's gene labels are permuted. The facets use free y axes because the statistics",
    "differ in scale, so only position relative to the curve and the triangle is comparable",
    "between them. The unsigned-regulon facet omits the regulons above the sweep's size cap.",
    "Annotation tier."),
  config = CFG, wide = TRUE
)

message("\n[DONE] 18_tf_activity_viz complete.")
