#!/usr/bin/env Rscript
# =============================================================================
# 15_coresh_search_viz.R  --  VIZ (no statistics)
# =============================================================================
# Two plain overview panels for the co-regulation search of the public human GEO
# compendium. Reads only tables written by 15_coresh_search.R; computes nothing.
#
#   figures/_overview/coresh_pctvar_top_hits.{pdf,png}
#       which public human datasets the JIA niche up-arm co-varies in, per query
#   figures/_overview/coresh_module_nes.{pdf,png}
#       where the modules those datasets define sit in each population's ranked list
#
# Run from the compartment root, AFTER 15_coresh_search.R:
#   Rscript 02_analysis/scripts/15_coresh_search_viz.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

source("02_analysis/helpers/figure_style.R")

STAGE  <- "15_coresh_search"
SCRIPT <- "02_analysis/scripts/15_coresh_search_viz.R"
TDIR   <- file.path("03_results", STAGE, "tables")

`%|N|%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

TOP_HITS_SHOWN <- 10L                      # bars per query panel; the search returns ~44k rows
NES_CAP  <- as.numeric(FIG_CFG$figures$nes_cap %|N|% 3.5)
GSEA_FDR <- as.numeric(FIG_CFG$thresholds$gsea_fdr %|N|% 0.05)
COL_UP   <- FIG_CFG$colors$diverging$up      %|N|% "#B35806"
COL_DOWN <- FIG_CFG$colors$diverging$down    %|N|% "#2166AC"
COL_MID  <- FIG_CFG$colors$diverging$neutral %|N|% "#F7F7F7"

POP_LABEL <- c(treg = "Treg", tcon = "Tcon", cd8 = "CD8")
GATE_LABEL <- c(fdr_only = "FDR only", fdr_logfc = "FDR + log2FC")

for (f in c("coresh_hits.csv", "coresh_derived_gsea.csv", "coresh_derived_annotation.csv",
            "coresh_gsea_summary.csv"))
  if (!file.exists(file.path(TDIR, f)))
    stop("[15_viz] ", f, " not found in ", TDIR, " -- run 15_coresh_search.R first.")

hits <- data.table::fread(file.path(TDIR, "coresh_hits.csv"))
gsea <- data.table::fread(file.path(TDIR, "coresh_derived_gsea.csv"))
ann  <- data.table::fread(file.path(TDIR, "coresh_derived_annotation.csv"))
summ <- data.table::fread(file.path(TDIR, "coresh_gsea_summary.csv"))

if (exists("purge_figures")) {
  purge_figures(STAGE, "coresh_pctvar_top_hits", overview = TRUE, config = FIG_CFG)
  purge_figures(STAGE, "coresh_module_nes",      overview = TRUE, config = FIG_CFG)
}

## ---------------------------------------------------------------------------
## Panel 1 -- the search result: top-ranked public datasets per query
## ---------------------------------------------------------------------------

prepare_pctvar_table <- function(hits, gsea, n_show = TOP_HITS_SHOWN) {
  d <- hits[rank <= n_show]
  became <- unique(gsea[, .(query_name, gse)])[, became_module := TRUE]
  d <- merge(d, became, by = c("query_name", "gse"), all.x = TRUE)
  d[is.na(became_module), became_module := FALSE]
  ## The query size is constant within a panel, so it belongs on the strip, not
  ## repeated on every bar where it costs the label its last characters.
  panel_lab <- function(pop, gate, q) sprintf("%s — %s   (query: %s Entrez ids)",
                                              POP_LABEL[pop], GATE_LABEL[gate],
                                              format(n_query_genes(q), big.mark = ","))
  lev <- unique(panel_lab(rep(c("treg", "tcon", "cd8"), each = 2),
                          rep(c("fdr_only", "fdr_logfc"), 3),
                          paste0(rep(c("treg", "tcon", "cd8"), each = 2), "_up_",
                                 rep(c("fdr_only", "fdr_logfc"), 3))))
  d[, `:=`(
    panel = factor(panel_lab(population, gate, query_name), levels = lev),
    bar_label = sprintf("%s  ·  %s matched", gse, format(size, big.mark = ",")))]
  ## Label offset is per panel: facets have free x scales, so one global nudge
  ## would push labels off the small-pctVar panels and crowd the large ones.
  d[, label_x := pctVar + max(pctVar) * 0.02, by = panel]
  d[order(panel, rank)]
}

PROV <- data.table::fread(file.path(TDIR, "coresh_query_provenance.csv"))
n_query_genes <- local({
  lut <- stats::setNames(PROV$n_unique_entrez, PROV$query_name)
  function(q) unname(lut[q])
})

create_pctvar_overview <- function(d) {
  d <- copy(d)
  ## Bar identity rides on the y AXIS, not on a geom_text at the bar tip: the top
  ## bar of a free-scale panel reaches the panel edge, so an in-panel label is the
  ## one that gets clipped. An axis label gets its own gutter and cannot truncate.
  d[, row_id := paste(panel, sprintf("%02d", rank), sep = "|")]
  lut <- stats::setNames(d$bar_label, d$row_id)
  ggplot(d, aes(x = pctVar, y = stats::reorder(row_id, -rank), fill = became_module)) +
    geom_col(width = 0.72) +
    facet_wrap(~ panel, scales = "free", ncol = 2) +
    scale_fill_manual(values = c(`TRUE` = COL_UP, `FALSE` = "grey70"),
                      labels = c(`TRUE` = "became a module", `FALSE` = "ranked only"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    scale_y_discrete(labels = function(x) unname(lut[x])) +
    labs(x = "pctVar  (% of that dataset's variance explained by the query genes)",
         y = NULL,
         title = "Where the JIA synovial up-arm co-varies across ~44,000 public human GEO datasets",
         subtitle = "Top 10 datasets per query, ranked by co-regulation score") +
    project_theme(config = FIG_CFG) +
    theme(legend.position = "bottom")
}

PCT <- prepare_pctvar_table(hits, gsea)
pct_tbl <- PCT[, .(population, gate, query_name, rank, gse, gpl, pctVar,
                   query_genes_in_dataset = size, became_module)]

save_overview(
  create_pctvar_overview(PCT), STAGE, "coresh_pctvar_top_hits", table = pct_tbl,
  finding = sprintf(
    paste("Across the %s public human datasets searched, the JIA synovial-fluid-versus-blood up-arm",
          "reaches a best co-regulation score of pctVar %.1f%% under the stringent FDR-plus-log2FC",
          "gate (%s) and %.1f%% under the relaxed FDR-only gate (%s); the relaxed gate scores higher",
          "throughout on a query roughly three times larger, so pctVar is comparable only within a",
          "gate and the two rankings are read as separate searches rather than as one league table."),
    format(PROV$n_datasets_searched[1], big.mark = ","),
    max(PCT[gate == "fdr_logfc"]$pctVar), PCT[gate == "fdr_logfc"][which.max(pctVar)]$gse,
    max(PCT[gate == "fdr_only"]$pctVar),  PCT[gate == "fdr_only"][which.max(pctVar)]$gse),
  script = SCRIPT, fn = "create_pctvar_overview",
  config_kv = sprintf("coresh.top_n_hits=%s; coresh.species=human; coresh.pvalues=%s; n_hits_shown=%d",
                      FIG_CFG$coresh$top_n_hits %|N|% 5, FIG_CFG$coresh$pvalues %|N|% FALSE, TOP_HITS_SHOWN),
  input = "03_results/15_coresh_search/tables/coresh_hits.csv",
  how_to_read = paste(
    "One panel per query: a sorted population's up arm at one significance gate. Each bar",
    "is one public human GEO dataset; length is pctVar, the share of that dataset's",
    "variance the query genes jointly explain — a co-regulation score, unsigned and always",
    "positive, higher meaning tighter co-movement. The row label gives the accession and",
    "how many query Entrez ids that dataset measures, because pctVar is normalised by that",
    "count. Orange bars became modules; grey were ranked only. Panels have independent x",
    "ranges: pctVar is not comparable across query sizes, and the relaxed gate's query is",
    "far larger. Exploratory tier: a ranking of public datasets, not a test."),
  config = FIG_CFG, width = 15, height = 12)

## ---------------------------------------------------------------------------
## Panel 2 -- the modules those datasets define, scored back on the ranked lists
## ---------------------------------------------------------------------------

prepare_module_table <- function(gsea, ann) {
  terms <- unique(ann[, .(set_name, terms_aligned_with_query)])
  terms[, short_terms := {
    t <- strsplit(terms_aligned_with_query, " / ", fixed = TRUE)
    vapply(t, function(x) paste(utils::head(x[nzchar(x)], 2), collapse = ", "), character(1))
  }]
  d <- merge(gsea, terms, by.x = "pathway_id", by.y = "set_name", all.x = TRUE)
  d[is.na(short_terms), short_terms := ""]
  d[, module_label := sprintf("%s  ·  %s, %s%s", gse, POP_LABEL[seed_population],
                              GATE_LABEL[gate],
                              ifelse(nzchar(short_terms), sprintf("  ·  %s", short_terms), ""))]
  d[, `:=`(nes_clamped = pmax(pmin(nes, NES_CAP), -NES_CAP),
           neglog_padj = -log10(pmax(padj, 1e-300)),
           scored_in   = factor(POP_LABEL[population], levels = unname(POP_LABEL)))]
  d[order(seed_population, gate, hit_rank)]
}

create_module_nes_dotplot <- function(d) {
  ord <- unique(d[order(seed_population, gate, -hit_rank)]$module_label)
  d[, module_label := factor(module_label, levels = ord)]
  ggplot(d, aes(x = scored_in, y = module_label)) +
    geom_point(aes(size = neglog_padj, fill = nes_clamped,
                   colour = seeded_from_this_population), shape = 21, stroke = 1.1) +
    scale_fill_gradient2(low = COL_DOWN, mid = COL_MID, high = COL_UP, midpoint = 0,
                         limits = c(-NES_CAP, NES_CAP), name = "NES") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "seeded from this population",
                                   `FALSE` = "seeded elsewhere"), name = NULL) +
    scale_size_continuous(range = c(3, 11), name = expression(-log[10]~FDR)) +
    labs(x = "scored on this population's SF-vs-PB ranked list", y = NULL,
         title = "Co-regulation modules mined from public human data,\nscored back on the JIA niche contrast",
         subtitle = "Rows are modules, named by the public dataset they were mined from") +
    project_theme(config = FIG_CFG) +
    theme(legend.position = "right")
}

MOD <- prepare_module_table(gsea, ann)
mod_tbl <- MOD[, .(pathway_id, scored_in_population = population, contrast, nes, pvalue, padj,
                   set_size, direction, seed_population, gate, gse,
                   seeded_from_this_population, n_genes, n_seed_genes, frac_seed_genes,
                   pctVar, hit_rank, terms_aligned_with_query)]

treg_sum <- summ[population == "treg"]
save_overview(
  create_module_nes_dotplot(MOD), STAGE, "coresh_module_nes", table = mod_tbl,
  finding = sprintf(
    paste("In the Treg synovial-fluid-versus-blood ranking the %d modules split both ways",
          "(NES %.2f to %.2f; %d up and %d down at FDR < %.2f), and which way a module goes tracks",
          "how much of it is the query that seeded it rather than anything it newly recruited:",
          "Spearman rho = %.2f between NES and seed fraction, on a median seed fraction of only",
          "%.0f%%. The enrichment is therefore reporting seed content back to itself as much as it",
          "is reporting public biology, which is why nothing here is read as evidence."),
    treg_sum$n_modules, treg_sum$nes_min, treg_sum$nes_max,
    treg_sum$n_sig_up, treg_sum$n_sig_down, GSEA_FDR,
    treg_sum$spearman_nes_vs_frac_seed, 100 * treg_sum$median_frac_seed),
  script = SCRIPT, fn = "create_module_nes_dotplot",
  config_kv = sprintf("gsea_min_size=%s; gsea_max_size=%s; gsea_seed=%s; gsea_nperm=%s; nes_cap=%s; engine=clusterProfiler::GSEA(by=fgsea)",
                      FIG_CFG$thresholds$gsea_min_size %|N|% 5, FIG_CFG$thresholds$gsea_max_size %|N|% 500,
                      FIG_CFG$thresholds$gsea_seed %|N|% 123, FIG_CFG$thresholds$gsea_nperm %|N|% 100000, NES_CAP),
  input = "03_results/15_coresh_search/tables/coresh_derived_gsea.csv + coresh_derived_annotation.csv + coresh_gsea_summary.csv",
  how_to_read = paste(
    "Each row is one co-regulation module: the genes loading most strongly onto the query",
    "direction inside one public human GEO dataset. The row label gives that dataset's",
    "accession, the population and gate whose up arm seeded it, and the compendium's own",
    "sample-metadata terms tracking the query axis there — a descriptor of the public",
    "dataset, not of the JIA data. Columns are the three sorted populations, each scored",
    "on its own synovial-fluid-versus-blood ranked list. Fill is NES, orange positive and",
    "blue negative, clamped. Size is -log10 FDR. A black outline marks the circular cells,",
    "where a module is scored on the list that seeded it; grey outlines are the",
    "informative comparison. Exploratory tier."),
  config = FIG_CFG, width = 15, height = 13)

message("[15_viz] wrote 2 overview panel(s) for ", STAGE)
