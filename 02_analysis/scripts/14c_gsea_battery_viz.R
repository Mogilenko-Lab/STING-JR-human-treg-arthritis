#!/usr/bin/env Rscript
# 14c_gsea_battery_viz.R: VIZ. The per-database GSEA battery for the JIA niche contrast.
# =============================================================================
# One cell per (sorted population x gene-set collection) over the donor-pseudobulk
# synovial-fluid-versus-paired-blood contrast, drawn with the RNAseq-toolkit GSEA plotters so
# the panels read the same way as the mouse-anchor battery in
# mouse_anchor/03_results/06_gsea/figures/by_contrast/.
#
# Viz only: every number comes from the cached gseaResult objects 14_unbiased_enrichment.R
# wrote plus the published gsea_all.csv, and 03_results/master/ stays untouched.
#
# The panels of one cell, and the selection rule each one carries in its own subtitle:
#   dotplot            GeneRatio x -log10(pooled FDR), fill NES.   top N by pooled FDR
#   facet              the dotplot split by NES sign.              top N per direction by FDR
#   barplot            NES bars from zero.                         pooled FDR < cutoff, then |NES|
#   running_sum        running ES, hit ticks, ranked moderated t.   top N by |NES|
#   running_sum_focus  the same panel over RS_FOCUS_SETS.          named, no ranking
# The two ranking metrics disagree, so a set can sit inside one panel and outside another in
# the same cell. Each panel writes a same-stem CSV of the rows it drew and the rule that
# picked them, which makes an absence checkable against a ranking.
#
# Panels use `padj_pooled`, the Benjamini-Hochberg correction across every test asked of one
# population's ranked list. Within one collection it induces the same ordering as the
# per-database `padj`, so the printed values, point sizes and significance flag move alone,
# all in the conservative direction. The per-database value travels in every same-stem CSV
# under `padj_in_database`.
#
# Reads, read-only:
#   03_results/objects/14_gsea/<tag>__<database>.rds   cached gseaResult per cell
#   03_results/14_unbiased_enrichment/tables/gsea_all.csv         the pooled family
#   03_results/14_unbiased_enrichment/tables/gsea_<tag>_<db>.csv  the collection's roster
#   03_results/14_unbiased_enrichment/tables/geneset_alias_map.csv
#   03_results/14_unbiased_enrichment/tables/geneset_manifest.csv
#
# Writes:
#   03_results/14_unbiased_enrichment/figures/by_contrast/<pop>/<DB>/*.{pdf,png}
#   03_results/14_unbiased_enrichment/tables/by_contrast/<pop>/<DB>/*.csv
#   03_results/14_unbiased_enrichment/README.md   one caption block per collection
#
# Run from the compartment root, after 14_unbiased_enrichment.R:
#   Rscript 02_analysis/scripts/14c_gsea_battery_viz.R

# ============================================================================
# 0. Setup: style contract first, then the toolkit GSEA plotters
# ============================================================================

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_figure,
                                               # contrast_path, style_series, write_caption,
                                               # round_numeric_cols

# format_pathway_names.R is sourced ahead of the plotters that call format_pathway_name().
.RTK <- "01_modules/RNAseq-toolkit"
.GP  <- file.path(.RTK, "scripts", "GSEA", "GSEA_plotting")
source(file.path(.RTK, "scripts", "custom_minimal_theme.R"))  # custom_minimal_theme_with_grid()
source(file.path(.GP, "format_pathway_names.R"))              # format_pathway_name()
source(file.path(.GP, "gsea_plotting_utils.R"))               # smart_wrap(), get_db_plot_params()
source(file.path(.GP, "gsea_dotplot.R"))
source(file.path(.GP, "gsea_dotplot_facet.R"))
source(file.path(.GP, "gsea_barplot.R"))
source(file.path(.GP, "gsea_running_sum_plot.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(methods)
  library(enrichplot)   # gseaplot2()
  library(DOSE)         # gseaResult S4 class definition
})
options(stringsAsFactors = FALSE)

# ============================================================================
# 1. Constants, all read from analysis_config.yaml
# ============================================================================

STAGE  <- "14_unbiased_enrichment"
SCRIPT <- "02_analysis/scripts/14c_gsea_battery_viz.R"

CFG      <- FIG_CFG
RESULTS  <- CFG$paths$results %||% "03_results/"
DIR_OBJ  <- CFG$paths$objects %||% "03_results/objects/"
DIR_GSEA <- file.path(DIR_OBJ, "14_gsea")
TBL      <- file.path(RESULTS, STAGE, CFG$paths$stage_tables_subdir %||% "tables")

FDR    <- as.numeric(CFG$thresholds$gsea_fdr %||% 0.05)
MINSZ  <- as.integer(CFG$thresholds$gsea_min_size %||% 5L)
MAXSZ  <- as.integer(CFG$thresholds$gsea_max_size %||% 500L)
TOPN   <- as.integer(CFG$figures$top_n %||% 20L)
RSTOP  <- as.integer(CFG$figures$running_sum_top %||% 5L)
RSYLIM <- as.numeric(unlist(CFG$figures$running_sum_ylim %||% c(-1, 1)))
NESCAP <- as.numeric(CFG$figures$nes_cap %||% 3.5)
SPECIES <- CFG$project$species %||% "Homo sapiens"

# The config's diverging scale, with no literal fallback, so a missing key stops the run
# keeping the scale identical to every other figure in the compartment.
NEG <- CFG$colors$diverging$down
MID <- CFG$colors$diverging$neutral
POS <- CFG$colors$diverging$up
if (any(vapply(list(NEG, MID, POS), function(x) is.null(x) || !nzchar(x), logical(1))))
  stop("[14c] analysis_config.yaml colors.diverging must define down, neutral and up.")

# Facet row cap per direction, below figures.top_n, which would stack forty wrapped set
# names on one canvas. The cap and the full per-direction counts go on the panel, and every
# set stays in the same-stem CSV.
FACET_TOPN <- 10L
DOT_WRAP   <- 50L    # the toolkit dotplot default; measured against for canvas width
FACET_WRAP <- 60L
BAR_WRAP   <- 55L    # gsea_barplot applies no wrap of its own
RS_WRAP    <- 46L    # running-sum legend: the set name, above its own statistics line
TITLE_WRAP <- 62L
SUB_WRAP   <- 92L

# Sorted populations: figure directory name -> ranked-list tag. Every directory holds the
# same donor-paired synovial fluid versus peripheral blood comparison, published in each
# emitted CSV as SF_vs_PB_<pop>.
POPS <- c(Treg = "treg", Tcon = "tcon", CD8 = "cd8")
contrast_label <- function(pop) sprintf("SF_vs_PB_%s", pop)

# Panels a collection gets, by how many sets it offers: 6 and above all four, 3 to 5 the
# dotplot and the running sum, 2 and below the running sum alone. Over a two-set collection
# the four panels are the same two points drawn four ways. The running sum is kept at every
# size, since it is the one panel that shows where in the ranking a set's genes sit.
FULL_BATTERY_MIN <- 6L
DOTPLOT_MIN      <- 3L

CFG_KV <- sprintf(
  "thresholds.gsea_fdr=%.2g; thresholds.gsea_min_size=%d; thresholds.gsea_max_size=%d; figures.top_n=%d; facet_top_n_per_direction=%d; figures.running_sum_top=%d; figures.running_sum_ylim=[%.1f,%.1f]; figures.nes_cap=%.1f; colors.diverging",
  FDR, MINSZ, MAXSZ, TOPN, FACET_TOPN, RSTOP, RSYLIM[1], RSYLIM[2], NESCAP)

message("=================================================================")
message("14c_gsea_battery_viz: per-database GSEA panels for the JIA niche contrast")
message(sprintf("  fdr=%.2g  top_n=%d  facet_top_n=%d  running_sum_top=%d",
                FDR, TOPN, FACET_TOPN, RSTOP))
message("=================================================================")

# ============================================================================
# 2. Guards: the sweep must have run, and the cache must match the tables
# ============================================================================

sweep_fp <- file.path(TBL, "gsea_all.csv")
if (!file.exists(sweep_fp))
  stop("[14c] ", sweep_fp, " not found. Run 02_analysis/scripts/14_unbiased_enrichment.R first.")
if (!dir.exists(DIR_GSEA))
  stop("[14c] cached gseaResult directory not found: ", DIR_GSEA,
       ". Run 02_analysis/scripts/14_unbiased_enrichment.R first.")

SWEEP <- readr::read_csv(sweep_fp, show_col_types = FALSE, progress = FALSE)
req <- c("population", "database", "pathway_id", "nes", "pvalue", "padj",
         "padj_pooled", "set_size", "n_tests_in_db", "n_tests_pooled")
miss <- setdiff(req, colnames(SWEEP))
if (length(miss) > 0)
  stop("[14c] gsea_all.csv is missing required columns: ", paste(miss, collapse = ", "))

MANIFEST <- readr::read_csv(file.path(TBL, "geneset_manifest.csv"),
                            show_col_types = FALSE, progress = FALSE)

# The alias resolution, as the sweep recorded it: for a set id present in more than one
# collection, which copy entered the pooled family (`kept_copy`) and which did not
# (`dropped_copy`). A dropped copy carries no pooled adjusted p of its own, and this map
# is what lets a cell draw it under the adjusted p it was actually pooled under.
ALIAS <- readr::read_csv(file.path(TBL, "geneset_alias_map.csv"),
                         show_col_types = FALSE, progress = FALSE)
alias_req <- c("population", "pathway_id", "kept_copy", "dropped_copy",
               "gene_content_identical")
alias_miss <- setdiff(alias_req, colnames(ALIAS))
if (length(alias_miss) > 0)
  stop("[14c] geneset_alias_map.csv is missing required columns: ",
       paste(alias_miss, collapse = ", "))

DBS <- MANIFEST$database
DBS <- DBS[order(-MANIFEST$n_sets_offered_for_pooling)]   # largest collection first

message(sprintf("[2] sweep: %d rows, %d populations, %d collections",
                nrow(SWEEP), dplyr::n_distinct(SWEEP$population),
                dplyr::n_distinct(SWEEP$database)))
for (pop in names(POPS)) {
  s <- SWEEP[SWEEP$population == pop, ]
  message(sprintf("    %-4s pooled family %s tests, %s at pooled FDR < %.2g",
                  pop, format(s$n_tests_pooled[1], big.mark = ","),
                  format(sum(s$padj_pooled < FDR, na.rm = TRUE), big.mark = ","), FDR))
}

# ============================================================================
# 3. Rehydrating one cell into a plottable gseaResult
# ============================================================================
## The toolkit plotters need a clusterProfiler `gseaResult` whose @result carries
## ID/Description/NES/p.adjust/setSize/core_enrichment and whose @geneList and @geneSets
## drive the running-sum curve. The cached object holds all of it, so the only surgery is
## the adjusted-p column: @result$p.adjust and @result$qvalue take the published
## `padj_pooled`. Every other statistic stays the cached object's own, which is what keeps a
## panel in agreement with gsea_all.csv.
##
## A cell draws its whole collection, alias copies included. A set id present in two
## collections is pooled once, under whichever collection the config declared first, and the
## copy in the other collection carries no pooled adjusted p of its own. Those copies are
## drawn with the adjusted p of the copy that was pooled, which the alias map licenses by
## verifying the gene content is identical, so the two copies are one hypothesis. That
## widens two collections, project_frozen from 1 set drawn to 7 and HSR_lens from 1 to 2,
## and `pooled_under` records the borrow on the panel and in every same-stem CSV. The pooled
## family, its size and every adjusted p in gsea_all.csv stay as the sweep left them.
##
## The borrow covers the adjusted p alone. The NES stays the alias copy's own, from the
## fgsea run behind this collection's cached object, so it can differ from the canonical
## copy's in the third decimal: HALLMARK_HYPOXIA in Treg is 2.2553 here and 2.2563 under
## Hallmark, permutation noise on the same ranking and the same genes.
## geneset_alias_map.csv publishes that gap per set as `abs_nes_difference`. Copying the
## canonical NES across would put the panel at odds with the object it was drawn from, which
## is the disagreement the freshness gate below exists to catch.

.CELL_CACHE <- new.env(parent = emptyenv())

#' Read one cached gseaResult and re-key its adjusted p to the pooled correction.
#' Returns NULL when the cell has no cache file or no drawable row.
as_gsearesult <- function(pop, db) {
  key <- paste(pop, db)
  if (!is.null(.CELL_CACHE[[key]])) {
    out <- .CELL_CACHE[[key]]
    if (isTRUE(out$.absent)) return(NULL)
    return(out$g)
  }
  tag <- POPS[[pop]]
  fp  <- file.path(DIR_GSEA, sprintf("%s__%s.rds", tag, db))
  if (!file.exists(fp)) {
    .CELL_CACHE[[key]] <- list(.absent = TRUE)
    return(NULL)
  }
  cached <- readRDS(fp)

  # gsea_all.csv holds the pooled family, so an alias copy is absent from it altogether.
  # The collection's own per-database table lists every set the collection was scored with,
  # alias copies included and flagged, so that table is the roster a cell draws from.
  per_db_fp <- file.path(TBL, sprintf("gsea_%s_%s.csv", tag, db))
  if (!file.exists(per_db_fp)) {
    .CELL_CACHE[[key]] <- list(.absent = TRUE)
    return(NULL)
  }
  roster <- readr::read_csv(per_db_fp, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::transmute(pathway_id, nes_published = nes,
                     is_alias = if ("is_pooled_alias" %in% names(.)) is_pooled_alias else FALSE)

  # This population's pooled values, keyed by (id, collection), so a borrowed adjusted p is
  # read from the published table.
  pop_pooled <- SWEEP %>%
    dplyr::filter(population == pop, !is.na(padj_pooled)) %>%
    dplyr::select(pathway_id, database, padj_pooled)
  # Which of this collection's ids were pooled elsewhere, and under which collection.
  borrow <- ALIAS %>%
    dplyr::filter(population == pop, dropped_copy == db,
                  is.na(gene_content_identical) | gene_content_identical) %>%
    dplyr::select(pathway_id, pooled_under = kept_copy)

  pooled <- roster %>%
    dplyr::left_join(borrow, by = "pathway_id") %>%
    # Which collection this row's pooled value is read from: its own, or for an alias copy
    # the identical-content canonical copy's.
    dplyr::mutate(pooled_under = dplyr::if_else(is_alias & !is.na(pooled_under),
                                                pooled_under, db)) %>%
    dplyr::left_join(pop_pooled, by = c("pathway_id", "pooled_under" = "database")) %>%
    dplyr::filter(!is.na(padj_pooled))

  r <- cached@result
  r <- r[r$ID %in% pooled$pathway_id, , drop = FALSE]
  if (nrow(r) == 0) {
    .CELL_CACHE[[key]] <- list(.absent = TRUE)
    return(NULL)
  }
  m <- match(r$ID, pooled$pathway_id)

  # Freshness gate: the cache and the published table have to be from one run. Drift here
  # means the tables were regenerated against a different object, and every quotable number
  # on the panels would be stale.
  dmax <- max(abs(r$NES - pooled$nes_published[m]), na.rm = TRUE)
  if (dmax > 1e-6)
    stop(sprintf(paste0("[14c] %s x %s: the cached gseaResult disagrees with gsea_all.csv on NES ",
                        "(max |difference| %.3g). The object cache and the published tables are ",
                        "from different runs. Re-run 02_analysis/scripts/14_unbiased_enrichment.R ",
                        "before drawing anything from them."), pop, db, dmax))

  res <- data.frame(
    ID              = r$ID,
    Description     = r$ID,       # raw id, so format_pathway_name() renders it identically everywhere
    setSize         = as.integer(r$setSize),
    enrichmentScore = as.numeric(r$enrichmentScore),
    NES             = as.numeric(r$NES),
    pvalue          = as.numeric(r$pvalue),
    p.adjust        = as.numeric(pooled$padj_pooled[m]),
    qvalue          = as.numeric(pooled$padj_pooled[m]),
    rank            = NA_integer_,
    leading_edge    = NA_character_,
    core_enrichment = as.character(r$core_enrichment),
    padj_in_database = as.numeric(r$p.adjust),
    # Where this row's pooled adjusted p came from. A name other than `db` marks an alias
    # copy drawn with the identical-content canonical copy's value.
    pooled_under     = as.character(pooled$pooled_under[m]),
    stringsAsFactors = FALSE
  )
  rownames(res) <- res$ID

  sets <- lapply(cached@geneSets, function(g) intersect(as.character(g), names(cached@geneList)))
  sets <- sets[intersect(res$ID, names(sets))]

  g <- methods::new("gseaResult",
    result      = res,
    organism    = SPECIES,
    setType     = "UNKNOWN",
    geneSets    = sets,
    geneList    = cached@geneList,
    keytype     = "UNKNOWN",
    permScores  = matrix(),
    params      = list(pvalueCutoff = 1, eps = 0, pAdjustMethod = "BH",
                       exponent = 1, minGSSize = MINSZ, maxGSSize = MAXSZ),
    gene2Symbol = character(),
    readable    = FALSE)
  .CELL_CACHE[[key]] <- list(.absent = FALSE, g = g)
  g
}

# ============================================================================
# 4. On-panel text and canvas sizing
# ============================================================================

wrap_text <- function(x, width) paste(strwrap(as.character(x), width = width), collapse = "\n")

SIGN_LINE <- "NES > 0 toward synovial fluid, NES < 0 toward paired blood."

# The toolkit plotters title the size legend "-log10(q-value)", which reads as
# clusterProfiler's own q column, while on these panels that column holds the sweep-wide
# pooled adjusted p. The legend is retitled to say so and match the subtitle. The
# override.aes is the toolkit's own, carried through here so the legend bubbles stay filled.
SIZE_LAB   <- expression(-log[10]("pooled FDR"))
SIZE_GUIDE <- guides(size = guide_legend(
  title = SIZE_LAB,
  override.aes = list(shape = 16, fill = "black", color = "black")))

cell_title <- function(pop, db)
  wrap_text(sprintf("%s synovial fluid versus paired blood, %s", pop, db), TITLE_WRAP)

#' Render a p-value for on-figure text: fixed above 0.001, else one-digit scientific.
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p >= 0.001, sprintf("%.3f", p), sprintf("%.0e", p)))

#' The display names a panel will carry, wrapped the way the toolkit plotters wrap them.
#' Measurement only; the plotters do their own wrapping from the same algorithm.
display_labels <- function(ids, width) {
  nm <- format_pathway_name(ids, use_formatting = TRUE, strip_prefix = TRUE)
  vapply(nm, smart_wrap, character(1), width = width, USE.NAMES = FALSE)
}

#' Size the canvas to the labels it will carry: width from the longest wrapped line, height
#' from the total line count. On a fixed 8.5 x 6.5 canvas the long-named collections come out
#' unreadable, the Reactome barplot squeezing its panel to nothing and clipping its title.
panel_geometry <- function(labels, legend_width = 2.4, min_height = 6.5,
                           max_width = 16, max_height = 18) {
  parts     <- strsplit(as.character(labels), "\n", fixed = TRUE)
  n_lines   <- sum(lengths(parts))
  max_chars <- suppressWarnings(max(nchar(unlist(parts)), 0L))
  w <- 0.075 * max_chars + 3.4 + legend_width          # label column + panel + legend
  h <- 0.30 * n_lines + 2.4                            # rows + title/subtitle/axis chrome
  c(width  = min(max(w, CFG$figures$width %||% 8.5), max_width),
    height = min(max(h, min_height), max_height))
}

#' A titled panel carrying one legible statement, for a cell where the panel's own selection
#' rule leaves nothing to draw. save_figure(void = TRUE) strips the axis chrome and grid at
#' save time, which keeps the style contract to its single entry point.
empty_state_panel <- function(ttl, sub, msg) {
  ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
    geom_blank() +
    annotate("label", x = 0, y = 0, label = msg, lineheight = 1.15,
             size = (CFG$figures$label_size %||% 4) + 1, fontface = "bold",
             colour = MID, fill = NEG, label.size = 0, label.r = unit(3, "pt"),
             label.padding = unit(0.8, "lines")) +
    labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
    project_theme(config = CFG)
}

# ============================================================================
# 5. The same-stem source table for one panel
# ============================================================================
## Each panel owes the numbers behind it at a readable size: the rows it drew, in draw
## order, with the rule that picked them and both adjusted-p columns. That table is also
## what makes the selection rules checkable against each other inside one cell.

panel_table <- function(g, ids, pop, db, panel, rule) {
  r <- g@result[match(ids, g@result$ID), , drop = FALSE]
  n_le <- vapply(strsplit(as.character(r$core_enrichment), "/", fixed = TRUE),
                 function(x) length(x[nzchar(x)]), integer(1))
  tibble::tibble(
    population        = pop,
    contrast          = contrast_label(pop),
    database          = db,
    panel             = panel,
    selection_rule    = rule,
    rank_in_panel     = seq_along(ids),
    pathway_id        = r$ID,
    display_name      = format_pathway_name(r$ID, use_formatting = TRUE, strip_prefix = TRUE),
    direction         = ifelse(r$NES > 0, "up", "down"),
    nes               = r$NES,
    enrichment_score  = r$enrichmentScore,
    pvalue            = r$pvalue,
    padj_in_database  = r$padj_in_database,
    padj_pooled       = r$p.adjust,
    # Where padj_pooled came from; a value other than `database` marks an alias copy.
    pooled_under      = r$pooled_under,
    set_size          = r$setSize,
    leading_edge_size = n_le,
    gene_ratio        = n_le / r$setSize)
}

write_panel_table <- function(tbl, pop, db, stem) {
  d <- file.path(contrast_path(STAGE, pop, "tables", config = CFG), db)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(round_numeric_cols(tbl), file.path(d, sprintf("%s.csv", stem)))
}

# ============================================================================
# 6. Per-cell emitter
# ============================================================================

EMITTED <- list()   # one row per panel written, for the run summary and the README

# The named pair a collection holding both should draw on its own panel, in this order.
RS_FOCUS_SETS <- c("HSR_core", "HALLMARK_HYPOXIA")

#' Draw one running-sum panel for a named set of ids and write its same-stem table.
#'
#' Every number goes on the face. The toolkit plotter carries no p-value table
#' (`pvalue_table = FALSE`), so a curve from it locates a set's genes and leaves strength and
#' significance to the CSV, and for a set resting on six genes that is the difference between
#' a result and a shape. The legend label therefore carries NES, the pooled adjusted p and
#' how many of the set's genes reached the ranked list, on a second line so the label column
#' stays narrow. This matches the running sums published by 05_scoring and 10_hsr_lens.
emit_running_sum <- function(g, ids, pop, db, stem, rule, ttl, note) {
  r <- g@result[match(ids, g@result$ID), , drop = FALSE]
  disp <- display_labels(ids, RS_WRAP)
  # An alias copy says so on its own legend line, so a borrowed adjusted p reads as borrowed.
  borrowed <- !is.na(r$pooled_under) & r$pooled_under != db
  labs <- stats::setNames(
    sprintf("%s\nNES %+.2f,  pooled FDR %s,  %d genes ranked%s",
            disp, r$NES, fmt_p(r$p.adjust), as.integer(r$setSize),
            ifelse(borrowed, sprintf(",  pooled under %s", r$pooled_under), "")),
    ids)

  # gseaplot2 draws its legend from @result$Description, so the display names go onto a copy.
  # format_pathway_name() runs once per name, and the panels above already formatted the
  # raw-id Description themselves.
  g_rs <- g
  g_rs@result$Description <- format_pathway_name(g_rs@result$ID, use_formatting = TRUE,
                                                 strip_prefix = TRUE)
  # max_name_length sits past the longest label, so the wrapping above owns the line breaks
  # and the plotter's own strwrap leaves them alone.
  p_rs <- tryCatch(
    gsea_running_sum_plot(g_rs, gene_set_ids = ids, labels = labs, title = ttl,
                          max_name_length = max(nchar(labs)) + 1L),
    error = function(e) {
      message(sprintf("  [14c] %s skipped (%s / %s): %s", stem, pop, db,
                      conditionMessage(e)))
      NULL
    })
  if (is.null(p_rs)) return(invisible(NULL))

  p_rs <- style_series(p_rs, ylim = RSYLIM, config = CFG)
  # The legend goes underneath the stack. style_series collects it outside-right, which suits
  # a two-curve figure; here five labels of a name plus a statistics line took three fifths
  # of the canvas and squeezed the enrichment curve into the left third. `&` reaches every
  # panel without indexing into the patchwork, so the collected guide moves as one.
  #
  # The per-panel themes and the plot-level theme both have to say "bottom": patchwork places
  # a collected guide from the plot-level theme in plot_annotation(), so setting it per-panel
  # with `&` alone leaves the legend on the right.
  LEGEND_BELOW <- ggplot2::theme(legend.position = "bottom",
                                 legend.justification = "left",
                                 legend.direction = "vertical",
                                 legend.key.width = unit(1.6, "lines"))
  p_rs <- p_rs & LEGEND_BELOW
  p_rs <- p_rs + patchwork::plot_annotation(
    caption = wrap_text(sprintf("%s Selection: %s. Running enrichment score clamped to [%.0f, %.0f], the range shared by every running sum of this battery",
                                SIGN_LINE, rule, RSYLIM[1], RSYLIM[2]), SUB_WRAP),
    theme = project_theme(config = CFG) + LEGEND_BELOW)
  # Height grows with the legend, a stack of two-line entries under the panels.
  chars <- max(nchar(unlist(strsplit(labs, "\n", fixed = TRUE))))
  geo <- c(width  = min(max(0.075 * chars + 3.0, 11), 15),
           height = 7.2 + 0.62 * length(ids))
  save_figure(p_rs, STAGE, file.path(db, stem), contrast = pop, config = CFG,
              width = geo[["width"]], height = geo[["height"]])
  write_panel_table(panel_table(g, ids, pop, db, stem, rule), pop, db, stem)
  note(stem, rule, geo[["width"]], geo[["height"]], length(ids))
  invisible(TRUE)
}

emit_cell <- function(pop, db) {
  g <- as_gsearesult(pop, db)
  if (is.null(g)) {
    message(sprintf("  [14c] SKIP %s / %s: no cached object or no set in the pooled family", pop, db))
    return(invisible(NULL))
  }

  n_sets <- nrow(g@result)
  n_sig  <- sum(g@result$p.adjust < FDR, na.rm = TRUE)
  n_up   <- sum(g@result$NES > 0, na.rm = TRUE)
  n_dn   <- sum(g@result$NES < 0, na.rm = TRUE)

  want_dot  <- n_sets >= DOTPLOT_MIN
  want_full <- n_sets >= FULL_BATTERY_MIN
  message(sprintf("  [14c] %-5s / %-16s %5d sets, %4d at pooled FDR < %.2g -> %s",
                  pop, db, n_sets, n_sig, FDR,
                  if (want_full) "dotplot, facet, barplot, running_sum"
                  else if (want_dot) "dotplot, running_sum" else "running_sum"))

  ttl <- cell_title(pop, db)

  # save_figure's same-stem purge stops at the stage directory, so the <DB>/ subdirectory is
  # cleared by hand and the run owns this cell's namespace.
  fig_dir <- file.path(contrast_path(STAGE, pop, "figures", config = CFG), db)
  tab_dir <- file.path(contrast_path(STAGE, pop, "tables",  config = CFG), db)
  for (d in c(fig_dir, tab_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    stale <- list.files(d, pattern = "\\.(png|pdf|csv)$", full.names = TRUE)
    if (length(stale)) file.remove(stale)
  }

  note <- function(panel, rule, w, h, n_rows)
    EMITTED[[length(EMITTED) + 1L]] <<- tibble::tibble(
      population = pop, database = db, panel = panel, selection_rule = rule,
      n_rows_drawn = n_rows, n_sets_in_cell = n_sets, width_in = w, height_in = h)

  # ---- 1. dotplot: top N by pooled FDR ------------------------------------
  if (want_dot) {
    # A collection smaller than the cap gets a rule naming the whole collection, since "top
    # 20" would claim a selection that never happened.
    rule <- if (n_sets <= TOPN) sprintf("every one of the %d sets in this collection", n_sets)
            else sprintf("top %d by pooled FDR", TOPN)
    ids  <- g@result$ID[order(g@result$p.adjust)][seq_len(min(TOPN, n_sets))]
    geo  <- panel_geometry(display_labels(ids, DOT_WRAP), legend_width = 2.6)
    sub  <- wrap_text(sprintf("%s Selection: %s. Black outline: pooled FDR < %.2g",
                              SIGN_LINE, rule, FDR), SUB_WRAP)
    p <- gsea_dotplot(g, filterBy = "p.adjust", showCategory = TOPN, padj_cutoff = FDR,
                      title = ttl, wrap_width = DOT_WRAP,
                      neg_color = NEG, mid_color = MID, pos_color = POS,
                      nes_limits = c(-NESCAP, NESCAP)) +
      labs(subtitle = sub) +
      project_theme(config = CFG) +
      SIZE_GUIDE
    if (n_sig == 0)
      p <- p + annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
                        label = sprintf("No set reaches pooled FDR < %.2g in this cell", FDR),
                        size = (CFG$figures$label_size %||% 4) + 1, fontface = "bold",
                        colour = "white", fill = NEG, label.size = 0,
                        label.r = unit(2, "pt"), label.padding = unit(0.4, "lines"))
    save_figure(p, STAGE, file.path(db, "dotplot"), contrast = pop, config = CFG,
                width = geo[["width"]], height = geo[["height"]])
    write_panel_table(panel_table(g, ids, pop, db, "dotplot", rule), pop, db, "dotplot")
    note("dotplot", rule, geo[["width"]], geo[["height"]], length(ids))
  }

  # ---- 2. facet: top N per direction by pooled FDR -------------------------
  if (want_full) {
    rule <- sprintf("top %d per direction by pooled FDR", FACET_TOPN)
    pick <- function(sel) {
      idx <- which(sel)
      idx <- idx[order(g@result$p.adjust[idx])]
      g@result$ID[idx[seq_len(min(FACET_TOPN, length(idx)))]]
    }
    ids <- c(pick(g@result$NES > 0), pick(g@result$NES < 0))
    geo <- panel_geometry(display_labels(ids, FACET_WRAP), legend_width = 2.6,
                          min_height = 7.5)
    sub <- wrap_text(sprintf("%s Selection: %s, out of %d up and %d down in this cell",
                             SIGN_LINE, rule, n_up, n_dn), SUB_WRAP)
    p <- gsea_dotplot_facet(g, showCategory = FACET_TOPN, padj_cutoff = FDR, title = ttl,
                            wrap_width = FACET_WRAP,
                            neg_color = NEG, mid_color = MID, pos_color = POS,
                            nes_limits = c(-NESCAP, NESCAP)) +
      labs(subtitle = sub) +
      project_theme(config = CFG) +
      SIZE_GUIDE
    save_figure(p, STAGE, file.path(db, "facet"), contrast = pop, config = CFG,
                width = geo[["width"]], height = geo[["height"]])
    write_panel_table(panel_table(g, ids, pop, db, "facet", rule), pop, db, "facet")
    note("facet", rule, geo[["width"]], geo[["height"]], length(ids))
  }

  # ---- 3. barplot: significant only, then top N by |NES| -------------------
  if (want_full) {
    sig_idx <- which(g@result$p.adjust < FDR)
    ids <- g@result$ID[sig_idx[order(abs(g@result$NES[sig_idx]), decreasing = TRUE)]]
    ids <- ids[seq_len(min(TOPN, length(ids)))]
    rule <- if (n_sig == 0)
              sprintf("sets at pooled FDR < %.2g only, and none reaches it in this cell", FDR)
            else if (n_sig <= TOPN)
              sprintf("every one of the %d sets at pooled FDR < %.2g", n_sig, FDR)
            else sprintf("pooled FDR < %.2g only, then top %d of those %d by |NES|",
                         FDR, TOPN, n_sig)
    sub <- wrap_text(sprintf("%s Selection: %s", SIGN_LINE, rule), SUB_WRAP)
    if (length(ids) == 0) {
      p   <- empty_state_panel(ttl, sub, sprintf(
        "No set reaches pooled FDR < %.2g in this cell.\nThis panel draws significant sets only.\nThe dotplot and facet panels show the full ranking.", FDR))
      geo <- c(width = CFG$figures$width %||% 8.5, height = CFG$figures$height %||% 6.5)
    } else {
      geo <- panel_geometry(display_labels(ids, BAR_WRAP), legend_width = 1.8)
      p <- gsea_barplot(g, padj_cutoff = FDR, top_n = TOPN, title = ttl,
                        neg_color = NEG, mid_color = MID, pos_color = POS,
                        nes_limits = c(-NESCAP, NESCAP)) +
        # gsea_barplot applies no wrap of its own, and an unwrapped Reactome or GO name
        # squeezes the bar panel to nothing. Wrapping keeps every character on the canvas.
        scale_x_discrete(labels = function(l) stringr::str_wrap(l, width = BAR_WRAP)) +
        labs(subtitle = sub) +
        project_theme(config = CFG)
    }
    save_figure(p, STAGE, file.path(db, "barplot"), contrast = pop, config = CFG,
                width = geo[["width"]], height = geo[["height"]],
                void = (length(ids) == 0))
    if (length(ids) > 0)
      write_panel_table(panel_table(g, ids, pop, db, "barplot", rule), pop, db, "barplot")
    note("barplot", rule, geo[["width"]], geo[["height"]], length(ids))
  }

  # ---- 4. running sum(s) ----------------------------------------------------
  have <- g@result$ID %in% names(g@geneSets)
  cand <- g@result[have, , drop = FALSE]
  if (nrow(cand) > 0) {
    rule <- if (nrow(cand) <= RSTOP)
              sprintf("every one of the %d set(s) in this collection", nrow(cand))
            else sprintf("top %d by |NES|", RSTOP)
    ids  <- cand$ID[order(abs(cand$NES), decreasing = TRUE)][seq_len(min(RSTOP, nrow(cand)))]
    emit_running_sum(g, ids, pop, db, "running_sum", rule, ttl, note)

    # A second panel over a named pair, where the collection holds one. Reading the curated
    # heat-shock core against HALLMARK_HYPOXIA is the comparison project_frozen exists to
    # make, and on the collection-wide panel those two are two of seven curves that a |NES|
    # ranking may leave out. Named explicitly, so the pair is stable across
    # populations, and it carries nothing the collection panel lacks.
    focus <- intersect(RS_FOCUS_SETS, cand$ID)
    if (length(focus) >= 2)
      emit_running_sum(
        g, focus, pop, db, "running_sum_focus",
        sprintf("named pair, no ranking: %s", paste(focus, collapse = " and ")),
        wrap_text(sprintf("%s synovial fluid versus paired blood, %s against %s",
                          pop, focus[1], focus[2]), TITLE_WRAP),
        note)
  } else {
    message(sprintf("  [14c] running_sum skipped (%s / %s): no set carries a membership list",
                    pop, db))
  }

  invisible(TRUE)
}

# ============================================================================
# 7. Main loop
# ============================================================================

message(sprintf("[7] drawing %d populations x %d collections", length(POPS), length(DBS)))
for (pop in names(POPS)) {
  message(sprintf("[7] == %s ==", pop))
  for (db in DBS) emit_cell(pop, db)
}

EMIT <- dplyr::bind_rows(EMITTED)
if (nrow(EMIT) == 0) stop("[14c] no panel was emitted at all.")

# ============================================================================
# 8. Per-collection README captions
# ============================================================================
## One caption block per collection, covering that collection's cells across all three
## populations. Each names the selection rule its panels use, since that rule governs every
## absence a reader sees, and each says what the collection contains, since a set enriching
## is a statement about gene content.

DB_CONTENT <- c(
  Hallmark = paste(
    "Fifty broad MSigDB Hallmark programs, the collection a bench reader can name.",
    "HALLMARK_HYPOXIA and both interferon-response sets are members, so this is the",
    "one panel of the battery where the hypoxia and interferon readings sit on a",
    "single axis with everything else."),
  KEGG = paste(
    "KEGG canonical pathways from MSigDB C2. The configuration asks for CP:KEGG and",
    "msigdbr 26 resolves that to the CP:KEGG_LEGACY subcollection; the resolved name",
    "is recorded in geneset_manifest.csv."),
  Reactome = paste(
    "Reactome canonical pathways from MSigDB C2. Set names here are the longest in",
    "the battery and are wrapped in full, so the label column is wide."),
  WikiPathways = "WikiPathways canonical pathways from MSigDB C2.",
  GO_BP = paste(
    "Gene Ontology biological process terms from MSigDB C5, the largest collection in",
    "the battery. The gap between a per-database adjusted p and the pooled adjusted p",
    "is smallest here, because this collection supplies about half of the pooled family."),
  GO_MF = "Gene Ontology molecular function terms from MSigDB C5.",
  GO_CC = "Gene Ontology cellular component terms from MSigDB C5.",
  TF_Targets = paste(
    "CollecTRI regulons, one unsigned gene set per transcription factor with",
    "activating and repressing targets pooled, so a set is that factor's",
    "transcriptional neighbourhood. A set enriching says the factor's targets move",
    "with one side of the contrast; it is a statement about target-gene expression",
    "and the factor's own activity is untested."),
  project_frozen = paste(
    "The seven frozen curated lists this compartment owns, all seven drawn. Six are",
    "re-pins of MSigDB Hallmark sets with identical gene content, so they are pooled",
    "once under Hallmark — geneset_manifest.csv records them as",
    "n_sets_aliased_out_of_pooling = 6 — and each is drawn here with the adjusted p it",
    "was pooled under, marked `pooled under Hallmark` on the panel and in the",
    "pooled_under column of every same-stem CSV. A set drawn twice is tested once.",
    "The seventh is HSR_core, the curated heat-shock-response lens, held",
    "independent of the mouse anchor and general to proteotoxic stress. This is also the",
    "one collection that holds both HSR_core and HALLMARK_HYPOXIA, so it carries the",
    "extra running_sum_focus panel over that named pair: the inflamed niche imposes",
    "proteotoxic and low-oxygen stress together and cross-sectional human data leave",
    "them entangled, so the pair is named and both curves are drawn, which a |NES|",
    "ranking would leave to chance."),
  mouse_projection = paste(
    "The three mouse-derived up arms projected onto human symbols: WT_heat_up",
    "WT_heat_up, KO_heat_up and Interaction_up. They are ordinary members of the sweep,",
    "each with the same standing as any other set, and WT_heat_up doubles as the",
    "reproduction check against the published targeted result on the same ranked list."),
  sting_axes = paste(
    "The two frozen axes from the SAVI positive-control compartment, sting_specific_up",
    "and ifn_only_up, which separate STING-attributable content from generic type-I",
    "interferon content."),
  MitoPathways = paste(
    "MitoPathways 3.0 from MitoCarta, 119 sets over 1,035 unique human symbols. The",
    "human list is the SOURCE here and the mouse-anchor copy of this same database is the",
    "homologene conversion of it, which is the opposite direction to every other",
    "cross-species asset in this compartment: the sets scored here are the less derived",
    "of the two. The catalogue is hierarchical and set names carry their depth, so",
    "MITOPATHWAYS_Metabolism.Detoxification.ROS_and_glutathione_metabolism is a child of",
    "MITOPATHWAYS_Metabolism.Detoxification and its genes are a subset of the parent's.",
    "Parent and child are scored as separate sets and both enter the pooled family, so a",
    "nested pair enriching together is one observation counted twice. The upstream",
    "catalogue holds 149 sets: 30 of 1 to 4 genes fall under the 5-gene floor the toolkit",
    "applies when it writes the file, and they are small nested leaves whose genes all",
    "survive in their parents, so the unique-symbol count is 1,035 either way."),
  HSR_lens = paste(
    "The curated heat-shock-response lens, frozen in this compartment from human MSigDB",
    "sets: HSR_core is the 56-symbol cleaned cytosolic subset, HSR_sensitivity the",
    "176-symbol full union. Both are drawn. HSR_core is also the seventh file of",
    "project_frozen, which is where it is pooled, so it appears here with that adjusted p,",
    "marked `pooled under project_frozen`; HSR_sensitivity is the one new hypothesis this",
    "collection contributes to the pooled family. The lens is independent",
    "of the mouse anchor and general to proteotoxic stress, since HSF1 also fires on",
    "oxidative, proteasome and metal stress."),
  TCR_activation = paste(
    "The curated TCR and immediate-early T-cell activation lens, 66 human symbols",
    "spanning TCR-proximal signalling, early costimulation, immediate-early",
    "transcription factors and activation effectors, with FOXP3 absent because it marks",
    "Treg lineage identity. The curated panel is defined in human",
    "symbols and the mouse-anchor list is the ortholog conversion of it, so the set scored",
    "here carries no ortholog step. Held as the activation-pole comparator against the",
    "heat-shock-response lens.")
)
DB_CONTENT_GENERIC <- paste(
  "Read each set's direction on its own. The panels rank sets within a cell; which",
  "regulator drives any of them is untested here.")

#' Which panels this collection got, and why any were left out.
panels_note <- function(db) {
  got <- sort(unique(EMIT$panel[EMIT$database == db]))
  n   <- max(EMIT$n_sets_in_cell[EMIT$database == db], na.rm = TRUE)
  if (all(c("dotplot", "facet", "barplot", "running_sum") %in% got))
    return(sprintf("All four panels are drawn for this collection (%d sets in the pooled family).", n))
  sprintf(paste0("Panels drawn for this collection: %s. The collection offers %d set(s) to the ",
                 "pooled family, and over a collection that small a dotplot, a direction split ",
                 "and a significant-only barplot are the same few points drawn several ways, so ",
                 "the redundant panels are left out. The running sum is kept at every size ",
                 "because it is the only panel that shows where in the ranking a set's genes sit."),
          paste(got, collapse = ", "), n)
}

for (db in DBS) {
  if (!db %in% EMIT$database) next
  body <- if (db %in% names(DB_CONTENT)) unname(DB_CONTENT[db]) else DB_CONTENT_GENERIC

  # A worked example of the two ranking metrics disagreeing, so the reader meets it in the
  # caption, on one panel of one cell.
  rank_note <- ""
  if (identical(db, "Hallmark")) {
    h <- SWEEP %>% dplyr::filter(population == "Treg", database == "Hallmark",
                                 !is.na(padj_pooled))
    if (nrow(h) > 0 && "HALLMARK_HYPOXIA" %in% h$pathway_id) {
      r_nes <- which(h$pathway_id[order(-abs(h$nes))] == "HALLMARK_HYPOXIA")[1]
      r_p   <- which(h$pathway_id[order(h$padj_pooled)] == "HALLMARK_HYPOXIA")[1]
      rank_note <- sprintf(paste0(
        " Worked example of the two rankings diverging: in the Treg Hallmark cell ",
        "HALLMARK_HYPOXIA (NES %+.4f, pooled FDR %s, %d genes in the ranked list) is %d by ",
        "|NES| and %d by pooled adjusted p, so it sits inside the dotplot's top %d and ",
        "outside the running sum's top %d. An absence from one panel is a statement about ",
        "that panel's ranking metric."),
        h$nes[h$pathway_id == "HALLMARK_HYPOXIA"][1],
        fmt_p(h$padj_pooled[h$pathway_id == "HALLMARK_HYPOXIA"][1]),
        as.integer(h$set_size[h$pathway_id == "HALLMARK_HYPOXIA"][1]),
        r_nes, r_p, TOPN, RSTOP)
    }
  }

  n_cells <- dplyr::n_distinct(EMIT$population[EMIT$database == db])
  write_caption(
    stage    = STAGE,
    filename = sprintf("figures/by_contrast/&lt;population&gt;/%s/*.png", db),
    finding  = sprintf(paste0(
      "%s GSEA of the donor-pseudobulk synovial-fluid-versus-paired-blood contrast in %d ",
      "sorted populations, drawn with the RNAseq-toolkit plotters on the cached gseaResult ",
      "from 03_results/objects/14_gsea/ with the adjusted p re-keyed to the sweep-wide pooled ",
      "correction published in gsea_all.csv. %s %s"),
      db, n_cells, body, panels_note(db)),
    script    = SCRIPT,
    fn        = "gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot",
    config_kv = CFG_KV,
    input     = sprintf("03_results/objects/14_gsea/{treg,tcon,cd8}__%s.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv", db),
    # Glyphs, sign convention and the pooled-correction note are stated once, in the
    # battery-level caption, so each collection's block carries only what differs.
    how_to_read = paste0(sprintf(paste0(
      "Selection rules, which govern every absence: dotplot top %d by pooled adjusted p; ",
      "facet top %d per direction by pooled adjusted p; barplot sets at pooled FDR < %.2g only, ",
      "then top %d of those by |NES|; running_sum top %d by |NES|; running_sum_focus, where it ",
      "exists, a named pair and no ranking at all. Glyphs, the sign convention ",
      "and the pooled correction are described once in the `figures/by_contrast/ (per-database GSEA battery)` section of this README. Each panel ",
      "writes its own same-stem CSV under tables/by_contrast/&lt;population&gt;/%s/ listing the ",
      "rows it drew, in draw order, with the rule that picked them, the per-collection ",
      "adjusted p under padj_in_database, and under pooled_under the collection each pooled ",
      "adjusted p was read from."),
      TOPN, FACET_TOPN, FDR, TOPN, RSTOP, db),
      rank_note,
      " A set enriching says its gene content moves with one side of this contrast. Correlative. ",
      "Claim tier: L3 (enrichment statistics).")
  )
}

# ============================================================================
# 9. Battery-level README caption
# ============================================================================

n_fig  <- nrow(EMIT)
skipped <- EMIT %>% dplyr::group_by(database) %>%
  dplyr::summarise(n_panels = dplyr::n_distinct(panel),
                   n_sets = max(n_sets_in_cell), .groups = "drop") %>%
  dplyr::filter(n_panels < 4)
skipped_txt <- if (nrow(skipped) == 0) "Every collection carries all four panels." else
  sprintf(paste0("%d collection%s too small for the full battery and carr%s fewer panels: %s. ",
                 "The omission is a redundancy judgement recorded per collection above, and no ",
                 "statistic is withheld: every set of every collection is in gsea_all.csv."),
          nrow(skipped), if (nrow(skipped) == 1) " is" else "s are",
          if (nrow(skipped) == 1) "ies" else "y",
          paste(sprintf("%s (%d set(s), %d panel type(s))", skipped$database, skipped$n_sets,
                        skipped$n_panels), collapse = "; "))

treg <- SWEEP %>% dplyr::filter(population == "Treg")
wt   <- treg %>% dplyr::filter(pathway_id == "WT_heat_up")
hyp  <- treg %>% dplyr::filter(database == "Hallmark", pathway_id == "HALLMARK_HYPOXIA")

write_caption(
  stage    = STAGE,
  filename = "figures/by_contrast/ (per-database GSEA battery)",
  finding  = sprintf(paste0(
    "The full browse surface for the unbiased sweep: %d panels across %d sorted populations and ",
    "%d gene-set collections, one directory per (population, collection) cell. Read against the ",
    "whole family, the Treg contrast carries %s pooled-significant sets out of %s tests, so the ",
    "mouse-derived WT_heat_up arm (NES %+.4f, pooled FDR %s, %d genes in the ranked list) ",
    "and HALLMARK_HYPOXIA (NES %+.4f, pooled FDR %s, %d genes) both sit inside a broad ",
    "co-enrichment shared with many other sets. %s"),
    n_fig, dplyr::n_distinct(EMIT$population), dplyr::n_distinct(EMIT$database),
    format(sum(treg$padj_pooled < FDR, na.rm = TRUE), big.mark = ","),
    format(treg$n_tests_pooled[1], big.mark = ","),
    wt$nes[1], fmt_p(wt$padj_pooled[1]), as.integer(wt$set_size[1]),
    hyp$nes[1], fmt_p(hyp$padj_pooled[1]), as.integer(hyp$set_size[1]),
    skipped_txt),
  script    = SCRIPT,
  fn        = "emit_cell",
  config_kv = CFG_KV,
  input     = "03_results/objects/14_gsea/*.rds + 03_results/14_unbiased_enrichment/tables/{gsea_all,gsea_&lt;population&gt;_&lt;collection&gt;,geneset_alias_map,geneset_manifest}.csv",
  how_to_read = sprintf(paste0(
    "Layout. figures/by_contrast/&lt;population&gt;/&lt;COLLECTION&gt;/{dotplot,facet,barplot,",
    "running_sum}.{pdf,png}, with the rows behind each panel in the mirrored path under ",
    "tables/by_contrast/. Population directories are Treg, Tcon and CD8, and the contrast inside ",
    "every one of them is the same donor-paired synovial fluid versus peripheral blood ",
    "comparison, published in each CSV as SF_vs_PB_&lt;population&gt;. ",
    "Where to start. The three Hallmark dotplots, one per population: fifty named programs on a ",
    "top-%d axis with both hypoxia and interferon among them. ",
    "Glyphs, shared by every cell of the battery. dotplot: x = GeneRatio (leading-edge genes ",
    "divided by set size), point size = -log10(pooled adjusted p), fill = NES with orange %s ",
    "positive and blue %s negative and the fill squished at plus or minus %.1f, black outline = ",
    "pooled FDR < %.2g. The dotplot selects by adjusted p and orders its y-axis by GeneRatio ",
    "descending, so vertical position there is a gene-ratio ranking. facet: the same dotplot ",
    "split into an NES > 0 block and an NES < 0 block. barplot: NES bars from zero, ordered by ",
    "NES. running_sum: three stacked panels, the running enrichment score with its leading-edge ",
    "peak on top, gene-hit ticks at each member's rank in the middle, and the ranked moderated t ",
    "at the bottom, with the score clamped to [%.0f, %.0f] so curves stay comparable between ",
    "collections, and each curve's NES, pooled adjusted p and effective size in the legend beneath ",
    "so the panel can be read from its own face. Sign. NES > 0 means the set's genes concentrate on ",
    "the synovial-fluid side of the ranking and NES < 0 on the paired-blood side. ",
    "Adjusted p. Every panel uses the Benjamini-Hochberg correction across the whole family of ",
    "tests asked of one population's ranked list, which is stricter than a single-collection ",
    "correction; the per-collection value travels in each same-stem CSV under padj_in_database. ",
    "A set configured into two collections is pooled once and drawn in both, the second time with ",
    "the adjusted p of the copy that was pooled, labelled `pooled under &lt;collection&gt;` on the ",
    "panel and recorded in the pooled_under column; the gene content of the two copies is verified ",
    "identical in geneset_alias_map.csv, so this is one hypothesis shown twice and tested once. ",
    "Ranking. The panel types rank by different metrics, adjusted p for dotplot and facet ",
    "and |NES| for barplot and running_sum, so read an absence against the rule named in that ",
    "panel's own subtitle before reading it as a null; running_sum_focus ranks nothing and draws a ",
    "named pair. ",
    "This is a browse surface, wide on purpose and privileging no set, and the claim spine stays ",
    "the donor-pseudobulk effect sizes. Correlative: enrichment describes where a set's genes sit ",
    "in a ranking. Claim tier: L3 (enrichment statistics), and no row of this battery reaches an ",
    "effect-size accumulator."),
    TOPN, POS, NEG, NESCAP, FDR, RSYLIM[1], RSYLIM[2])
)

# ============================================================================
# 10. Run summary and structural asserts
# ============================================================================

fig_root <- file.path(RESULTS, STAGE, "figures", CFG$figures$by_contrast_dir %||% "by_contrast")
tab_root <- file.path(RESULTS, STAGE, "tables",  CFG$figures$by_contrast_dir %||% "by_contrast")
figs <- list.files(fig_root, pattern = "\\.(pdf|png)$", recursive = TRUE, full.names = TRUE)
tabs <- list.files(tab_root, pattern = "\\.csv$",       recursive = TRUE, full.names = TRUE)
bytes <- sum(file.size(figs), na.rm = TRUE) + sum(file.size(tabs), na.rm = TRUE)

message("")
message(sprintf("14c_gsea_battery_viz complete."))
message(sprintf("  %d panels, %d figure files, %d source tables, %.1f MB total",
                nrow(EMIT), length(figs), length(tabs), bytes / 1024^2))
print(EMIT %>% dplyr::count(database, panel) %>%
        tidyr::pivot_wider(names_from = panel, values_from = n, values_fill = 0L),
      n = 40)

stopifnot("no figure file produced" = length(figs) >= 1)
stopifnot("panel count does not match figure files" = length(figs) == 2L * nrow(EMIT))
