#!/usr/bin/env Rscript
# render.R — render a breakpoint review notebook (.qmd) to gfm + self-contained html.
# Adapted from mouse_anchor/02_analysis/notebooks/render.R. Uses rmarkdown + system
# pandoc (no Quarto CLI needed). These notebooks are R chunks that READ the stage's
# published CSV tables and embed the already-rendered figure PNGs — they never compute
# anything a later stage depends on (steering gates only).
#
# Usage (from the compartment root):
#   Rscript 02_analysis/notebooks/render.R [notebook.qmd] [gfm|html|both]
# Defaults: 05_gonogo_review/05_gonogo_review.qmd, both.

args       <- commandArgs(trailingOnly = TRUE)
default_nb <- "02_analysis/notebooks/05_gonogo_review/05_gonogo_review.qmd"
nb   <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_nb
what <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "both"
stopifnot("notebook not found" = file.exists(nb),
          "target must be gfm|html|both" = what %in% c("gfm", "html", "both"))

if (what %in% c("gfm", "both"))
  rmarkdown::render(nb, output_format = rmarkdown::github_document(html_preview = FALSE), quiet = TRUE)
if (what %in% c("html", "both"))
  rmarkdown::render(nb, output_format = rmarkdown::html_document(self_contained = TRUE,
                                                                toc = TRUE, toc_float = TRUE), quiet = TRUE)
message("rendered ", nb, "  ->  ", what)
