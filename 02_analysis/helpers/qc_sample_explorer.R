# qc_sample_explorer.R — human_treg_arthritis · 02_analysis/helpers
# =============================================================================
# build_qc_sample_explorer()
#
# Reusable, repurposable helper: renders an interactive per-GSM QC explorer for
# the sorted JIA Treg/Tcon/CD8 libraries (GSE160097). Each of the 40 GSMs is
# aggregated to a SINGLE point so the whole cohort — and the one dropped sample —
# is legible at a glance. Carved out as a standalone helper (mirroring
# 02_analysis/helpers/runsum_widget.R) so ANY chunk / ANY qmd can source this file
# and call the function instead of copying plotly code. It reads only a published
# per-cell QC table and never recomputes anything a downstream stage depends on
# (compute never plots; viz never computes).
#
# @param path A length-1 character path to the per-cell QC metrics table
#   `03_results/01_qc/tables/qc_metrics_per_cell.csv`. Expected columns:
#   gsm, donor, tissue, population, population_short, total_counts,
#   n_genes_by_counts, pct_counts_mt, pct_counts_ribo, doublet_score,
#   predicted_doublet, mad_outlier, low_genes, excluded_gsm, pass_qc.
# @return A self-contained plotly htmlwidget (survives `embed-resources: true`),
#   or invisibly NULL (with a Quarto callout-warning emitted) if the table is
#   missing.
#
# AGGREGATION — one row per GSM: n_cells (barcodes), n_kept (pass_qc == TRUE),
# frac_kept, and the medians of total_counts / n_genes_by_counts / pct_counts_mt.
# donor / tissue / population / population_short carry through unchanged (constant
# within a GSM). A GSM is flagged EXCLUDED when any of its cells is marked
# `excluded_gsm` (the whole-library hard drop).
#
# ENCODING — x = median genes/cell (LOG axis), y = median %mt, point SIZE = n_cells
# (area-scaled), COLOUR = sorted population (Treg / Tcon / CD8; same Okabe-Ito
# palette as the running-sum widget + NES forest). Rich per-GSM hover: gsm, donor,
# tissue, population, n_cells, frac_kept, median counts / genes / %mt, and a
# kept/EXCLUDED flag. The excluded near-empty library (GSM4859852) is over-drawn
# with a black cross and a pinned annotation so it stands out at far-left
# low-median-genes as the UMI/gene outlier it is.
#
# LAYOUT — compact fixed height (~420). Following the runsum_widget.R fix, the
# legend and any control live in SEPARATE regions: the population legend is a
# horizontal strip along the BOTTOM (below the x-axis title), leaving the plot
# area and top margin uncluttered so no control overlaps it.
# =============================================================================

#' The one population palette, read from analysis_config.yaml::colors.populations, so this
#' widget shows Treg, Tcon and CD8 in the same hues as the static figures and cannot drift
#' from them. Resolved here rather than sourced from the figure-style shim because this
#' helper is designed to be sourced on its own from any qmd, and the shim pulls the whole
#' plotting contract lib with it. Path is compartment-root relative, matching how a qmd
#' sources this file.
read_population_palette <- function(
    path = "02_analysis/config/analysis_config.yaml") {
  cfg   <- yaml::read_yaml(path)
  okabe <- cfg$colors$okabe_ito
  named <- cfg$colors$populations
  if (is.null(named) || length(named) == 0)
    stop("colors.populations absent from ", path,
         " — the population palette has one home and this is it.", call. = FALSE)
  vapply(named, function(hue) as.character(okabe[[hue]]), character(1))
}

build_qc_sample_explorer <- function(path) {
  suppressPackageStartupMessages({ library(plotly); library(dplyr) })

  if (length(path) != 1 || !file.exists(path)) {
    cat(paste0(
      "::: {.callout-warning}\n**Interactive per-sample QC explorer not available.**\n\n",
      "Expected the per-cell QC table at `",
      path, "`.\n:::\n\n"))
    return(invisible(NULL))
  }

  pop_col <- read_population_palette()

  cells <- read.csv(path, stringsAsFactors = FALSE)
  cells$excluded_gsm <- as.logical(cells$excluded_gsm)
  cells$pass_qc      <- as.logical(cells$pass_qc)

  # ---- Aggregate per GSM to one point each. ----
  g <- cells %>%
    group_by(gsm) %>%
    summarise(
      donor            = donor[1],
      tissue           = tissue[1],
      population       = population[1],
      population_short = population_short[1],
      n_cells          = dplyr::n(),
      n_kept           = sum(pass_qc),
      frac_kept        = mean(pass_qc),
      med_counts       = median(total_counts),
      med_genes        = median(n_genes_by_counts),
      med_mt           = median(pct_counts_mt),
      excluded         = any(excluded_gsm),
      .groups = "drop") %>%
    mutate(
      population_short = factor(population_short, levels = names(pop_col)),
      donor_s          = sub("JIA_patient_", "JIA patient ", donor),
      tissue_s         = c(synovial_fluid = "synovial fluid",
                           peripheral_blood = "peripheral blood")[tissue],
      flag             = ifelse(excluded, "EXCLUDED", "kept"))

  # Area-scaled marker sizing: largest library ~= 26px, smallest still visible.
  max_px  <- 26
  sizeref <- 2 * max(g$n_cells) / (max_px^2)

  hov <- with(g, sprintf(paste0(
    "<b>%s</b>  (%s)<br>",
    "%s &middot; %s<br>",
    "population: %s<br>",
    "n_cells: %s &middot; kept: %.1f%%<br>",
    "median counts: %s<br>",
    "median genes: %s<br>",
    "median %%mt: %.1f%%"),
    gsm, flag, donor_s, tissue_s, population,
    format(n_cells, big.mark = ",", trim = TRUE), 100 * frac_kept,
    format(round(med_counts), big.mark = ",", trim = TRUE),
    format(round(med_genes), big.mark = ",", trim = TRUE),
    med_mt))

  # ---- Base scatter: one coloured trace per population (for the legend). ----
  p <- plot_ly(height = 420)
  for (pop in names(pop_col)) {
    d <- g[g$population_short == pop, ]
    if (!nrow(d)) next
    p <- add_trace(p,
      data = d, x = ~med_genes, y = ~med_mt,
      type = "scatter", mode = "markers",
      name = pop, legendgroup = pop,
      marker = list(color = pop_col[[pop]], sizemode = "area",
                    size = d$n_cells, sizeref = sizeref, sizemin = 4,
                    opacity = 0.82, line = list(width = 1, color = "#FFFFFF")),
      text = hov[g$population_short == pop], hoverinfo = "text")
  }

  # ---- Over-draw the excluded near-empty library so it stands out. ----
  ex <- g[g$excluded, ]
  if (nrow(ex)) {
    p <- add_trace(p,
      x = ex$med_genes, y = ex$med_mt,
      type = "scatter", mode = "markers",
      name = "excluded (near-empty)", legendgroup = "excluded",
      marker = list(color = "rgba(0,0,0,0)", symbol = "x-thin",
                    size = 15, line = list(width = 2.4, color = "#000000")),
      text = hov[g$excluded], hoverinfo = "text")
  }

  p <- layout(p,
    height = 420, autosize = TRUE,
    # ~52px top for the title; ~78px bottom for the x-axis title + the horizontal
    # legend strip beneath it (the runsum_widget.R separate-regions fix).
    margin = list(l = 62, r = 18, t = 52, b = 78),
    title = list(text = paste0(
      "Per-sample QC — 40 sorted JIA Treg / Tcon / CD8 libraries (one point per GSM)"),
      font = list(size = 12), y = 0.98, yanchor = "top"),
    xaxis = list(title = "median genes / cell (log scale)", type = "log",
                 automargin = TRUE),
    yaxis = list(title = "median mitochondrial %", automargin = TRUE,
                 rangemode = "tozero"),
    hovermode = "closest",
    legend = list(
      title = list(text = "Sorted population:"),
      orientation = "h", x = 0.5, xanchor = "center", y = -0.22, yanchor = "top",
      font = list(size = 11), tracegroupgap = 8),
    annotations = if (nrow(ex)) lapply(seq_len(nrow(ex)), function(i) list(
      x = log10(ex$med_genes[i]), y = ex$med_mt[i],
      text = sprintf("<b>%s</b><br>excluded: near-empty library<br>(~%s barcodes, median %d genes/cell)",
                     ex$gsm[i], format(ex$n_cells[i], big.mark = ",", trim = TRUE),
                     round(ex$med_genes[i])),
      showarrow = TRUE, arrowhead = 2, arrowsize = 1, ax = 55, ay = -42,
      font = list(size = 10, color = "#000000"),
      bgcolor = "rgba(255,255,255,0.82)", bordercolor = "#000000", borderwidth = 1)) else list())

  config(p, displayModeBar = TRUE, displaylogo = FALSE)
}
