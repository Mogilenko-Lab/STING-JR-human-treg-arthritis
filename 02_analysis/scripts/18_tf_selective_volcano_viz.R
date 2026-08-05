#!/usr/bin/env Rscript
# 18_tf_selective_volcano_viz.R: VIZ. One figure — the exclusively-claimed regulon targets
# placed on the volcano of the contrast they were scored on.
# =============================================================================
#   tf_selective_targets_volcano   the Treg SF-vs-PB donor-pseudobulk volcano with only the
#                                  targets no other CollecTRI regulon claims named
#
# WHY THIS FIGURE EXISTS, SEPARATELY FROM tf_selective_targets:
# `tf_target_promiscuity` reports that HIF1A's 27 exclusively-claimed targets sum to 0.14% of
# its signed contribution, and that number reads as "these genes do nothing". They in fact hold
# 15% of the regulon's signed total in MAGNITUDE and cancel, 13 going up on the synovial-fluid
# side against 14 going down. `tf_selective_targets` shows that cancellation in the regulon's
# own arithmetic; this figure shows the same genes in the coordinates a reader already trusts —
# effect size against evidence on the committed differential-expression table — so the split is
# visible as position in the contrast rather than as a property of the decomposition.
#
# The volcano grammar is the canonical RNAseq-toolkit plotter
# (01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R): raw p on the y axis for
# resolution, FDR for the significance decision, four significance categories by colour. It is
# called with top_n = 0 so it draws NO labels of its own, and the only named genes on the canvas
# are the exclusively-claimed targets. project_theme() is applied afterwards because
# project_theme is a complete theme and must be the last one added.
#
# Reads:
#   03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv        the committed limma-voom table
#   03_results/18_tf_activity/tables/target_decomposition.csv which targets are exclusive
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/18_tf_selective_volcano_viz.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_overview
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")  # create_standard_volcano

STAGE  <- "18_tf_activity"
SCRIPT <- "02_analysis/scripts/18_tf_selective_volcano_viz.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

CFG  <- FIG_CFG
TA   <- CFG$tf_activity
FIG  <- CFG$figures
THR  <- CFG$thresholds
OKI  <- unlist(CFG$colors$okabe_ito)
DIV  <- unlist(CFG$colors$diverging)

DECOMP  <- unlist(TA$decompose_tfs)
SEL_MAX <- as.integer(TA$selective_max_regulons %||% 1L)
POP     <- TA$primary_population %||% "treg"
LBL_SZ  <- as.numeric(FIG$label_size %||% 4)
PT_SZ   <- as.numeric(FIG$point_size %||% 2.4)
RPL_SEED <- as.integer(THR$gsea_seed %||% 123L)

# The volcano's decision rule is the compartment's committed DE rule, not a per-figure choice.
FDR_CUT   <- as.numeric(THR$de_fdr %||% 0.05)
LOGFC_CUT <- as.numeric(THR$de_logfc %||% 1.0)

RES <- CFG$paths$results %||% "03_results/"
TBL <- file.path(RES, STAGE, CFG$paths$stage_tables_subdir %||% "tables")

DE  <- readr::read_csv(file.path(RES, "03_pseudobulk", "tables",
                                sprintf("de_SFvsPB_%s.csv", POP)), show_col_types = FALSE)
DEC <- readr::read_csv(file.path(TBL, "target_decomposition.csv"), show_col_types = FALSE)

# `POP` is the file-name tag ("treg"); the display name is the label the compute stage recorded.
POP_LABEL <- DEC$population[1]

# --- the exclusively-claimed targets, per decomposed factor -------------------------------
sel <- DEC %>%
  filter(tf %in% DECOMP, n_regulons <= SEL_MAX) %>%
  mutate(tf = factor(tf, levels = DECOMP)) %>%
  select(tf, target, mor, sign_decision, stat, contrib, n_regulons)

# A target claimed by two of the decomposed factors would be plotted twice at one position and
# the shape legend would be a lie. They are disjoint here; stop if that ever changes.
stopifnot(!any(duplicated(sel$target)))

# --- the volcano substrate ----------------------------------------------------------------
# create_standard_volcano() reads limma's column names off the frame and takes gene ids from
# ROWNAMES, so the committed table is renamed rather than the plotter adapted. The DE table is
# already one row per HGNC symbol (13,999 rows, no duplicate symbol), which is what makes
# rownames a safe key here; assert it rather than assume it.
stopifnot(!any(duplicated(DE$gene_symbol)), !any(is.na(DE$gene_symbol)))

de_v <- DE %>%
  transmute(gene_symbol, logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
  as.data.frame()
rownames(de_v) <- de_v$gene_symbol
de_v$gene_symbol <- NULL

# Every exclusive target must be on the volcano, or the figure silently under-reports the set.
missing_targets <- setdiff(as.character(sel$target), rownames(de_v))
stopifnot(length(missing_targets) == 0)

hl <- sel %>%
  mutate(logFC      = de_v[as.character(target), "logFC"],
         P.Value    = de_v[as.character(target), "P.Value"],
         adj.P.Val  = de_v[as.character(target), "adj.P.Val"],
         neglog10p  = -log10(P.Value))

# --- per-factor counts for the in-panel text ----------------------------------------------
# Aggregation for display over a committed table: no statistic is computed here.
# The split is counted on log2FC, not on signed contribution, because log2FC is the axis this
# figure plots. The two disagree for any target carried on a repressing edge, where the recorded
# sign flips the contribution relative to the contrast — TM9SF4 is the only such target here, and
# it is why 13/14 by contribution reads as 12/15 by fold change.
hl_txt <- hl %>%
  group_by(tf) %>%
  summarise(n = n(),
            n_up = sum(logFC > 0), n_down = sum(logFC < 0),
            n_fdr = sum(adj.P.Val <= FDR_CUT),
            n_fdr_up = sum(adj.P.Val <= FDR_CUT & logFC > 0),
            n_fdr_down = sum(adj.P.Val <= FDR_CUT & logFC < 0),
            .groups = "drop") %>%
  arrange(tf) %>%
  transmute(line = sprintf("%s claims %d targets no other regulon does: %d up and %d down in this contrast; %s",
                           tf, n, n_up, n_down,
                           # A "splitting 0 up and 0 down" clause is vacuous, so it is only
                           # written when some target actually clears the cut.
                           ifelse(n_fdr == 0,
                                  sprintf("none clears FDR <= %.2g", FDR_CUT),
                                  sprintf("%d clear FDR <= %.2g, splitting %d up and %d down",
                                          n_fdr, FDR_CUT, n_fdr_up, n_fdr_down))))

# --- the figure ---------------------------------------------------------------------------
message("[fig] exclusively-claimed targets on the contrast volcano ...")

# top_n = 0 suppresses the plotter's own label layers, so the named genes are exactly the
# highlighted set and nothing competes with them for space.
p <- create_standard_volcano(
  de_v,
  decision_by     = "fdr",
  p_cutoff        = FDR_CUT,
  fc_cutoff       = LOGFC_CUT,
  top_n           = 0,
  annotate_counts = TRUE,
  title    = "Exclusively-claimed regulon targets on the contrast they were scored on",
  subtitle = paste(c(sprintf("Sorted %s, synovial fluid versus paired blood, donor pseudobulk (limma-voom)",
                             POP_LABEL),
                     hl_txt$line), collapse = "\n")
)

# The highlight is an unfilled black RING, not a coloured point. The plotter's top significance
# category is already vermillion and a burnt-orange fill sat indistinguishably on top of it — the
# highlighted genes vanished into the cloud. A ring reads against every one of the four category
# colours and leaves the point inside it visible, so a named gene still shows its own category.
# Labels are black on a white halo for the same reason: the dense part of this cloud is coloured.
p <- p +
  geom_point(data = hl, aes(x = logFC, y = neglog10p, shape = tf), inherit.aes = FALSE,
             fill = NA, colour = "black", size = PT_SZ * 2.1, stroke = 1.3) +
  geom_text_repel(data = hl, aes(x = logFC, y = neglog10p, label = target), inherit.aes = FALSE,
                  colour = "black", size = LBL_SZ, fontface = "bold",
                  bg.color = "white", bg.r = 0.18,
                  max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4,
                  segment.colour = "black", box.padding = 0.5, point.padding = 0.4, force = 12,
                  seed = RPL_SEED, show.legend = FALSE) +
  scale_shape_manual(values = setNames(c(21, 24)[seq_along(DECOMP)], DECOMP),
                     name = "Claimed by this regulon alone") +
  # project_theme is a COMPLETE theme, so it must come after the plotter's own theming or it
  # would be reset by it; anything element-level goes after this line.
  project_theme(config = CFG) +
  theme(legend.position = "bottom", legend.box = "vertical")

fig_tbl <- hl %>%
  transmute(population = POP_LABEL, tf, target, mor, sign_decision,
            stat, contrib, n_regulons, logFC, P.Value, adj.P.Val,
            past_fdr = adj.P.Val <= FDR_CUT,
            past_fdr_and_logfc = adj.P.Val <= FDR_CUT & abs(logFC) >= LOGFC_CUT) %>%
  arrange(tf, desc(contrib))

# Counts for the caption, read off the same frame the figure draws. `pct_abs`/`pct_net` come from
# the committed decomposition so the volcano's account of the share matches the decomposition's.
cnt <- hl %>% group_by(tf) %>%
  summarise(n = n(), n_up = sum(logFC > 0), n_down = sum(logFC < 0),
            n_fdr = sum(adj.P.Val <= FDR_CUT),
            n_fdr_up = sum(adj.P.Val <= FDR_CUT & logFC > 0),
            n_fdr_down = sum(adj.P.Val <= FDR_CUT & logFC < 0),
            n_contrib_up = sum(contrib > 0), n_contrib_down = sum(contrib < 0),
            .groups = "drop")
cn <- function(f, col) cnt[[col]][as.character(cnt$tf) == f]
tot <- DEC %>% filter(tf %in% DECOMP) %>% group_by(tf) %>%
  summarise(pct_abs = 100 * sum(abs(contrib[n_regulons <= SEL_MAX])) / sum(contrib),
            pct_net = 100 * sum(contrib[n_regulons <= SEL_MAX]) / sum(contrib),
            .groups = "drop")
tt <- function(f, col) tot[[col]][as.character(tot$tf) == f]
flipped <- hl %>% filter(mor < 0)

save_overview(
  p, STAGE, "tf_selective_targets_volcano", table = fig_tbl,
  finding = paste(
    sprintf(paste("Placed on the contrast they were scored on, the %d targets HIF1A alone claims",
                  "sit on both sides: %d go up in synovial fluid and %d go down, %d clear FDR %.2g,",
                  "and those %d split %d up against %d down, so the set carries evidence and",
                  "carries no direction."),
            cn("HIF1A", "n"), cn("HIF1A", "n_up"), cn("HIF1A", "n_down"),
            cn("HIF1A", "n_fdr"), FDR_CUT, cn("HIF1A", "n_fdr"),
            cn("HIF1A", "n_fdr_up"), cn("HIF1A", "n_fdr_down")),
    sprintf(paste("That is what makes its %.2f%% share of HIF1A's signed contribution a",
                  "cancellation: in magnitude the same %d targets are %.0f%% of that total."),
            tt("HIF1A", "pct_net"), cn("HIF1A", "n"), tt("HIF1A", "pct_abs")),
    sprintf(paste("NFKB1's %d exclusively-claimed targets split one each way, and the FDR cut is",
                  "cleared by neither."),
            cn("NFKB1", "n"))),
  script = SCRIPT, fn = "save_overview",
  config_kv = paste0("thresholds.de_fdr=", FDR_CUT, "; thresholds.de_logfc=", LOGFC_CUT,
                     "; tf_activity.decompose_tfs=[", paste(DECOMP, collapse = ", "),
                     "]; tf_activity.selective_max_regulons=", SEL_MAX,
                     "; tf_activity.primary_population=", POP),
  input = paste0("03_results/03_pseudobulk/tables/de_SFvsPB_", POP, ".csv, ",
                 "03_results/18_tf_activity/tables/target_decomposition.csv"),
  how_to_read = paste(
    "The standard volcano of the committed donor-pseudobulk contrast: x is log2 fold change,",
    "synovial fluid over paired blood; y is raw p on a -log10 scale; colour is the significance",
    "category, decided on FDR while the axis keeps raw p for resolution. The dashed horizontal rule",
    "is the raw p that realises the FDR cut, the vertical rules the fold-change cut.",
    "Only the targets no other CollecTRI regulon claims are ringed and named, circles for HIF1A and",
    "triangles for NFKB1, the ring unfilled so the point inside keeps its category colour. The",
    "spread of the named genes across both halves is the point, so nothing else is labelled: count",
    "how many sit either side of zero and how many clear the rules.",
    sprintf(paste("The split here is by fold change; by signed contribution it reads %d up and %d",
                  "down, since %s sits on a repressing edge and contributes positively while going",
                  "down."),
            cn("HIF1A", "n_contrib_up"), cn("HIF1A", "n_contrib_down"),
            paste(flipped$target, collapse = ", ")),
    "Annotation tier: this restates a committed DE table and a set membership, so nothing here is",
    "a new test."),
  config = CFG, wide = TRUE, height = 9
)

message("\n[DONE] 18_tf_selective_volcano_viz complete.")
