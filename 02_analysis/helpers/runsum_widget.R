# runsum_widget.R — human_treg_arthritis · 02_analysis/helpers
# =============================================================================
# build_runsum_widget()
#
# Reusable, repurposable helper: renders the interactive fgsea running-enrichment
# widget for the mouse WT_heat (39 °C Treg) signature scored in JIA SF-vs-PB
# donor-pseudobulk ranked lists, one running-ES line per sorted T-cell population
# (Treg / Tcon / CD8) and gene-set arm (up / down). Carved out of the inline
# `{r runsum-interactive}` chunk of 05_gonogo_review.qmd so ANY chunk / ANY qmd
# can source this file and call the function instead of copying plotly code.
#
# @param source Either (a) a length-1 character path to the tables directory
#   holding the six `runsum_interactive_{treg,tcon,cd8}_WT_heat_{up,down}.csv`
#   files, or (b) a named list / character vector of the six file paths keyed
#   `<pop>_<arm>` (e.g. `treg_up`, `cd8_down`). A directory is the common case;
#   the explicit-six form lets a caller point at relocated / renamed tables.
# @return A self-contained plotly htmlwidget (survives `embed-resources: true`),
#   or invisibly NULL (with a Quarto callout-warning emitted) if any of the six
#   source tables is missing.
#
# Each CSV has columns: rank, gene, stat, running_es, hit, leading_edge,
# gene_set, population, contrast. `running_es` is precomputed by the stage-05
# viz script to overlay the static running-sum figures exactly — this helper
# plots it verbatim and never recomputes (compute never plots; viz never
# computes).
#
# Renders: one running-ES curve per population; a per-set hit-rug (stacked ticks
# below that set's curves) marking where the mouse-anchor core signature genes
# land; hover on the rug showing gene symbol + SF-vs-PB signed moderated-t stat +
# leading-edge flag; an up/down gene-set dropdown; population toggling via the
# legend (legendgroup ties each population's curve + rug together).
#
# LAYOUT — kept deliberately IDENTICAL to the integration report's parallel copy
# at integration/03_results/final_report/helpers.R::build_runsum_widget()
# (dropdown top-LEFT, legend as a horizontal strip along the BOTTOM, compact
# height, title pinned top). The two are intentional parallel copies: they live
# in separate repos (this treg gitlink vs the integration layer), so a cross-repo
# source() is not appropriate — when you change the layout in one, mirror it in
# the other to keep them in sync.
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

build_runsum_widget <- function(source) {
  suppressPackageStartupMessages({ library(plotly) })

  pop_lab <- c(treg = "Treg", tcon = "Tcon", cd8 = "CD8")
  pop_col <- read_population_palette()
  sets    <- c(up = "WT_heat_up", down = "WT_heat_down")

  # ---- Resolve the six source paths from either a dir or an explicit set. ----
  keys <- c("treg_up", "treg_down", "tcon_up", "tcon_down", "cd8_up", "cd8_down")
  if (is.list(source) || (is.character(source) && length(source) > 1)) {
    src <- unlist(source)
    if (!all(keys %in% names(src)))
      stop("build_runsum_widget(): explicit paths must be named <pop>_<arm> (",
           paste(keys, collapse = ", "), ")")
    paths <- src[keys]
  } else {
    tdir  <- as.character(source)
    paths <- setNames(
      file.path(tdir, sprintf("runsum_interactive_%s.csv",
                              sub("_(up|down)$", "_WT_heat_\\1", keys))),
      keys)
  }

  if (!all(file.exists(paths))) {
    cat(paste0(
      "::: {.callout-warning}\n**Interactive running-sum not yet available.**\n\n",
      "Expected six `runsum_interactive_*` tables under ",
      "`human_treg_arthritis/03_results/05_scoring/tables/`.\n:::\n\n"))
    return(invisible(NULL))
  }

  # ---- Read the six tables into one long frame (data logic preserved). ----
  rs <- do.call(rbind, lapply(names(pop_lab), function(pt) {
    do.call(rbind, lapply(names(sets), function(sd) {
      d <- read.csv(paths[[paste0(pt, "_", sd)]], stringsAsFactors = FALSE)
      d$hit <- as.logical(d$hit); d$leading_edge <- as.logical(d$leading_edge)
      d$set <- sd; d
    }))
  }))

  # per-set rug band (below that set's own curves), one row per population.
  # Tight spacing (small offset + step) keeps the three rug rows compact so the
  # whole widget stays short without losing the per-population separation.
  rug_y <- lapply(names(sets), function(sd) {
    sub <- rs[rs$set == sd, ]
    rng <- range(sub$running_es); span <- max(diff(rng), 1e-3)
    base <- rng[1] - 0.04 * span; step <- 0.03 * span
    c(Treg = base, Tcon = base - step, CD8 = base - 2 * step)
  })
  names(rug_y) <- names(sets)

  p <- plot_ly(height = 400)
  for (sd in names(sets)) {
    vis <- (sd == "up")
    for (pt in names(pop_lab)) {
      pop <- pop_lab[[pt]]
      d <- rs[rs$population == pt & rs$set == sd, ]; d <- d[order(d$rank), ]
      p <- add_trace(p, x = d$rank, y = d$running_es, type = "scattergl", mode = "lines",
                     line = list(color = pop_col[[pop]], width = 2),
                     name = pop, legendgroup = pop, visible = vis,
                     showlegend = vis, hoverinfo = "skip")
      h <- d[d$hit, ]
      p <- add_trace(p, x = h$rank, y = rep(rug_y[[sd]][[pop]], nrow(h)),
                     type = "scatter", mode = "markers",
                     marker = list(color = pop_col[[pop]], symbol = "line-ns-open",
                                   size = 11, line = list(width = 1.6, color = pop_col[[pop]])),
                     name = paste(pop, "core genes"), legendgroup = pop,
                     visible = vis, showlegend = FALSE, hoverinfo = "text",
                     text = sprintf(paste0("<b>%s</b> (%s core gene)<br>",
                                           "population: %s<br>rank: %d<br>",
                                           "SF-vs-PB signed moderated t: %.2f<br>leading edge: %s"),
                                    h$gene, sets[[sd]], pop, h$rank, h$stat,
                                    ifelse(h$leading_edge, "yes", "no")))
    }
  }

  # up/down selector: 12 traces in order [up: L,R,L,R,L,R][down: L,R,L,R,L,R]
  sl_up   <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, rep(FALSE, 6))
  sl_down <- c(rep(FALSE, 6), TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)

  p <- layout(p,
    # Compact fixed height (~400), bumped slightly from the original so the title +
    # dropdown band above and the bottom legend strip never crowd the plot area.
    height = 400,
    autosize = TRUE,
    # Margins carve dedicated bands: ~64px top for title + dropdown, ~80px bottom
    # for the x-axis title plus the horizontal legend strip below it.
    margin = list(l = 58, r = 16, t = 64, b = 80),
    title = list(text = paste0(
      "Interactive running-sum — mouse WT_heat core gene set across sorted ",
      "Treg / Tcon / CD8 (SF-vs-PB)"),
      font = list(size = 12), y = 0.98, yanchor = "top"),
    xaxis = list(title = "rank in SF-vs-PB signed moderated-t list (1 = most SF-up)", automargin = TRUE),
    yaxis = list(title = "running enrichment score", automargin = TRUE),
    hovermode = "closest",
    # Legend as a horizontal strip along the BOTTOM (below the x-axis title) so it
    # never collides with the dropdown in the cramped top margin. groupclick keeps
    # each population's curve + hit-rug toggling together.
    legend = list(
      title = list(text = "Population:"),
      orientation = "h", groupclick = "togglegroup",
      x = 0.5, xanchor = "center", y = -0.22, yanchor = "top",
      font = list(size = 11), tracegroupgap = 4),
    # Dropdown pinned top-LEFT, in its own band under the title — a separate region
    # from the bottom legend, so the two never overlap.
    updatemenus = list(list(type = "dropdown", direction = "down",
      x = 0, y = 1.12, xanchor = "left", yanchor = "top", showactive = TRUE,
      buttons = list(
        list(method = "update", label = "Gene set: WT_heat_up",
             args = list(list(visible = c(rep(TRUE, 6), rep(FALSE, 6)), showlegend = sl_up))),
        list(method = "update", label = "Gene set: WT_heat_down",
             args = list(list(visible = c(rep(FALSE, 6), rep(TRUE, 6)), showlegend = sl_down)))))))
  config(p, displayModeBar = TRUE, displaylogo = FALSE)
}
