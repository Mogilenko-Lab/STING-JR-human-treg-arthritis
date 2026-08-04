## 02_analysis/helpers/figure_style.R — per-project figure-style shim.
## =====================================================================
## ONE import per viz script. Delegates to the SciAgent-toolkit contract
## lib (02_analysis/helpers/figure-style/figure_helpers.R, symlinked by
## `sciagent activate`) and loads the project analysis_config.yaml.
##
## Usage in any viz script:
##   source("02_analysis/helpers/figure_style.R")
##   p <- ggplot(...) + project_theme(config = FIG_CFG)
##   save_overview(p, "04_gsea", "name", table = df, ..., config = FIG_CFG)

## ---------------------------------------------------------------------------
## 1. Resolve the symlinked contract lib and source it (graceful fallback).
## ---------------------------------------------------------------------------
.HELPERS_LIB <- file.path(
    dirname(sys.frame(1L)$ofile %||% "02_analysis/helpers"),
    "figure-style",
    "figure_helpers.R"
)

# Portable %||% for use before the lib is sourced.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.FIGURE_STYLE_LOADED <- FALSE
if (file.exists(.HELPERS_LIB)) {
    source(.HELPERS_LIB)
    .FIGURE_STYLE_LOADED <- TRUE
} else {
    warning(
        "[figure_style] toolkit lib not found at: ", .HELPERS_LIB,
        "\nRun `sciagent activate` in this repo to link lib/figure-style/. ",
        "Falling back to minimal stubs.", call. = FALSE
    )
    # --- MINIMAL FALLBACK stubs (keep scripts from hard-crashing) -----------
    load_figure_config <- function(path = "02_analysis/config/analysis_config.yaml") {
        if (requireNamespace("yaml", quietly = TRUE)) yaml::read_yaml(path)
        else { warning("[figure_style] yaml unavailable; using empty config."); list() }
    }
    # `variant` accepted-but-ignored (drop-in compat with old call sites).
    project_theme <- function(base_size = NULL, legend = TRUE, variant = NULL,
                               config = NULL, ...) {
        if (!requireNamespace("ggplot2", quietly = TRUE))
            stop("[figure_style] ggplot2 required for project_theme().")
        ggplot2::theme_minimal(base_size = base_size %||% 14)
    }
    set_paper_style <- function(...) project_theme(...)
}

## ---------------------------------------------------------------------------
## 2. Load the project config once (FIG_CFG is the stable project-wide handle).
## ---------------------------------------------------------------------------
FIG_CFG <- tryCatch(
    load_figure_config("02_analysis/config/analysis_config.yaml"),
    error = function(e) {
        warning("[figure_style] Could not load analysis_config.yaml: ", conditionMessage(e),
                call. = FALSE)
        list()
    }
)

## ---------------------------------------------------------------------------
## 3. The two categorical axes (project-level; mirrors config.py's accessors).
## ---------------------------------------------------------------------------
## `colors.populations` and `colors.tissue` each map a level to an `okabe_ito` KEY rather
## than a hex literal, so a hue is resolved here and no viz script carries its own. The
## two blocks must stay disjoint: several figures draw both axes, and one puts a tissue
## panel and a cell-state panel in the same row, where a shared hue would make them
## indistinguishable. `.resolve_named_hues()` is shared so both resolve identically.
.resolve_named_hues <- function(config, block) {
    okabe <- config$colors$okabe_ito
    named <- config$colors[[block]]
    if (is.null(named) || length(named) == 0)
        stop("[figure_style] analysis_config.yaml::colors.", block, " is absent — every ",
             "figure that colours by ", block, " reads it, so a missing block is a config ",
             "error rather than something to fall back from.", call. = FALSE)
    bad <- setdiff(unlist(named, use.names = FALSE), names(okabe))
    if (length(bad))
        stop("[figure_style] colors.", block, " names ", paste(bad, collapse = ", "),
             ", which is not a key of colors.okabe_ito (",
             paste(names(okabe), collapse = ", "),
             "). A level names a palette entry, never a hex literal.", call. = FALSE)
    vapply(named, function(hue) as.character(okabe[[hue]]), character(1))
}

## Named character vector {population -> hex}. Keyed by BOTH spellings a script may hold:
## the short figure label (`Treg`) and the long FACS gate name (`CD4_Treg`), so either
## resolves to the same hue. `scale_*_manual()` ignores names absent from the data, so the
## extra aliases never reach a legend.
population_colors <- function(config = FIG_CFG) {
    short <- .resolve_named_hues(config, "populations")
    long  <- c(CD4_Treg = "Treg", CD4_Tcon = "Tcon", CD8 = "CD8")
    long  <- long[long %in% names(short)]
    c(short, setNames(unname(short[long]), names(long))[!names(long) %in% names(short)])
}

## Named character vector {tissue level -> hex}, warm for the inflamed joint and cool for
## paired blood. Keys are checked against `design.tissue_levels`, so renaming a tissue
## level cannot silently leave its colour behind under the old name.
tissue_colors <- function(config = FIG_CFG) {
    resolved <- .resolve_named_hues(config, "tissue")
    levels   <- unlist(config$design$tissue_levels, use.names = FALSE)
    if (length(levels) && !setequal(names(resolved), levels))
        stop("[figure_style] colors.tissue keys (", paste(names(resolved), collapse = ", "),
             ") do not match design.tissue_levels (", paste(levels, collapse = ", "),
             "). The two axes of this compartment are declared together; renaming a ",
             "tissue level means renaming its colour.", call. = FALSE)
    resolved
}
