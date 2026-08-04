#!/usr/bin/env Rscript
# 14_sweep_named_sets.R: COMPUTE (no plotting)
# =============================================================================
# The closing figure of this compartment's narrative asks a single question: scored
# against every set in eleven collections on the same donor-level ranked lists, with
# no set privileged, where do hypoxia and cGAS-STING land? The only real editorial
# decision in that figure is WHICH named sets it draws, so this script makes that
# decision explicit, writes it to a committed table with a reason per row, and hands
# the table to the viz script. Nothing is selected silently inside a plotting call.
#
# WHY A SCRIPT AND NOT A HARDCODED VECTOR IN THE VIZ. The obvious selection is a trap.
# A case-insensitive substring search for "STING" over the Treg sweep returns eight
# rows, and two of them are false positives that match the substring inside the word
# "RE-STING": REACTOME_PHASE_4_RESTING_MEMBRANE_POTENTIAL and
# GOBP_REGULATION_OF_RESTING_MEMBRANE_POTENTIAL are membrane-potential terms with no
# relation to the cGAS-STING axis. Both are recorded here with the reason they were
# dropped, and the audit is asserted rather than trusted, so a future collection that
# adds another "resting" term fails loudly instead of drawing a spurious point.
#
# The genuine cGAS-STING family in this sweep has SIX members, not the four a reader
# might expect. Two of them are regulation-of terms with opposite sign, and the
# positive-regulation term runs negative. All six are drawn: a family where the
# annotation's own sign structure disagrees with itself is part of the honest picture.
#
# Output (03_results/14_unbiased_enrichment/tables/):
#   sweep_named_sets.csv         the selection decision, one row per set considered,
#                                including the two rows that were excluded and why
#   sweep_named_sets_stats.csv   population x named set: NES, FDR, set size, and the
#                                set's rank inside its population's whole sweep
#   sweep_setsize_baseline.csv   population x set-size band: how often a set of that
#                                size reaches pooled significance at all
#
# Run from the compartment root, AFTER 14_unbiased_enrichment.R:
#   Rscript 02_analysis/scripts/14_sweep_named_sets.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, round_numeric_cols

STAGE <- "14_unbiased_enrichment"
TDIR  <- file.path(FIG_CFG$paths$results %||% "03_results/", STAGE,
                   FIG_CFG$paths$stage_tables_subdir %||% "tables")
FDR   <- as.numeric(FIG_CFG$thresholds$gsea_fdr %||% 0.05)
UE    <- FIG_CFG$unbiased_enrichment
POP_LEVELS <- c("Treg", "Tcon", "CD8")

sweep <- readr::read_csv(file.path(TDIR, "gsea_all.csv"),
                         show_col_types = FALSE, progress = FALSE)
stopifnot("gsea_all.csv must carry all three sorted populations" =
            setequal(unique(sweep$population), POP_LEVELS))

# =============================================================================
# 1. The cGAS-STING family, established by audit rather than by assumption
# =============================================================================
## The two documented false positives. Named as literals so the exclusion is a
## committed decision a reviewer can check, and asserted below so a renamed or
## newly added membrane-potential term cannot slip past.
STING_FALSE_POSITIVES <- c(
  "REACTOME_PHASE_4_RESTING_MEMBRANE_POTENTIAL",
  "GOBP_REGULATION_OF_RESTING_MEMBRANE_POTENTIAL")

treg_ids  <- sweep$pathway_id[sweep$population == "Treg"]
naive_hit <- sort(grep("STING", treg_ids, ignore.case = TRUE, value = TRUE))
sting_ids <- setdiff(naive_hit, STING_FALSE_POSITIVES)

message(sprintf("[named_sets] substring audit: %d ids contain 'STING', %d are membrane-potential false positives, %d genuine",
                length(naive_hit), length(intersect(naive_hit, STING_FALSE_POSITIVES)),
                length(sting_ids)))
stopifnot(
  "the two documented false positives must both still be present, or the exclusion list is stale" =
    all(STING_FALSE_POSITIVES %in% naive_hit),
  "every remaining STING substring hit must name the cGAS-STING axis; a new 'resting' term needs adding to STING_FALSE_POSITIVES" =
    all(grepl("CGAS_STING|STING_PATHWAY|STING_MEDIATED|STING_SPECIFIC",
              sting_ids, ignore.case = TRUE)),
  "the cGAS-STING family in this sweep is six sets" = length(sting_ids) == 6L)

# =============================================================================
# 2. The comparison sets, and one display label per row
# =============================================================================
## `thread` groups the rows into the four things this figure is asked to compare, and
## `thread_order` fixes their top-to-bottom order on the panel so a regenerated figure
## cannot silently reshuffle. Labels turn MSigDB underscores into spaces; the full
## identifier stays in `pathway_id`, and nothing is truncated.
## Only the all-caps MSigDB-style identifiers are prettified; the frozen list names
## (WT_heat_up, sting_specific_up) are the names those signatures are known by and stay
## verbatim, because renaming a signature in a display label is how a checkable name
## drifts into an unchecked one.
lbl <- function(id) ifelse(grepl("^[A-Z0-9_]+$", id), gsub("_", " ", id), id)

selection <- dplyr::bind_rows(
  tibble::tibble(
    pathway_id = c("WT_heat_up", "KO_heat_up", "Interaction_up"),
    display_label = pathway_id,
    thread = "mouse 39 °C-derived arms", thread_order = 1L,
    why_included = c(
      paste("the primary projected up arm, and the reason this compartment ran the",
            "sweep at all; drawn so its position can be read against every other set tested"),
      paste("the cGAS-knockout comparator for the same 39 °C contrast; drawn beside the",
            "primary arm because it beats it on pooled FDR in all three populations and",
            "shares most of its genes, so the two are close to one observation twice"),
      paste("the interaction arm of the same mouse contrast, 6 of 7 genes in the ranked",
            "list; drawn so the reader sees the size at which this sweep stops resolving"))),
  tibble::tibble(
    pathway_id = "HALLMARK_HYPOXIA", display_label = lbl(pathway_id),
    thread = "hypoxia", thread_order = 2L,
    why_included = paste("the curated versioned hypoxia set this compartment has used",
                         "throughout; the hypoxia half of the question this figure closes")),
  tibble::tibble(
    pathway_id = "ifn_only_up", display_label = pathway_id,
    thread = "type-I interferon", thread_order = 3L,
    why_included = paste("the generic type-I interferon axis frozen from the SAVI",
                         "compartment; carried so a cGAS-STING result can be read against",
                         "the interferon response it would be confused with")),
  tibble::tibble(
    pathway_id = sting_ids, display_label = lbl(pathway_id),
    thread = "cGAS-STING", thread_order = 4L,
    why_included = paste("a member of the complete cGAS-STING family in this sweep;",
                         "all six are drawn, including the two regulation-of terms whose",
                         "signs disagree, so no member is chosen for its result"))
) |>
  dplyr::mutate(included = TRUE, why_excluded = NA_character_)

excluded <- tibble::tibble(
  pathway_id = STING_FALSE_POSITIVES,
  display_label = lbl(pathway_id),
  thread = "cGAS-STING", thread_order = 4L,
  why_included = NA_character_, included = FALSE,
  why_excluded = paste("substring false positive: 'STING' matches inside 'RESTING'.",
                       "A membrane-potential term with no relation to the cGAS-STING axis"))

named <- dplyr::bind_rows(selection, excluded)

## A set is scored in a population only if enough of its genes reach that population's
## ranked list to clear gsea_min_size, so presence is a per-population fact and the
## count is recorded per row. GOBP_POSITIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY
## carries 5 testable genes in Treg and Tcon and fewer in CD8, so it was never tested
## there. The figure has to say "not tested" for that cell; drawing nothing and letting
## a reader infer a null would be the one reading the data cannot support.
tested_in <- sweep |> dplyr::group_by(pathway_id) |>
  dplyr::summarise(n_populations_tested = dplyr::n_distinct(population),
                   populations_tested = paste(sort(unique(population)), collapse = "|"),
                   .groups = "drop")
named <- named |> dplyr::left_join(tested_in, by = "pathway_id") |>
  dplyr::mutate(n_populations_tested = tidyr::replace_na(n_populations_tested, 0L),
                populations_tested = tidyr::replace_na(populations_tested, ""))
stopifnot("no set may be listed twice" = !any(duplicated(named$pathway_id)),
          "every selected set must at least be tested in Treg, the primary compartment" =
            all(grepl("Treg", named$populations_tested[named$included])))
n_partial <- sum(named$included & named$n_populations_tested < length(POP_LEVELS))
if (n_partial)
  message(sprintf("[named_sets] %d drawn set(s) were not tested in every population: %s",
                  n_partial,
                  paste(sprintf("%s (%s only)",
                                named$pathway_id[named$included & named$n_populations_tested < 3],
                                named$populations_tested[named$included & named$n_populations_tested < 3]),
                        collapse = "; ")))

## Record which collection each set came from, so a reader can see that the family
## spans GO_BP, Reactome, WikiPathways and a frozen list rather than one source.
coll <- sweep |> dplyr::filter(population == "Treg") |>
  dplyr::select(pathway_id, source_collection = database) |> dplyr::distinct()
named <- named |> dplyr::left_join(coll, by = "pathway_id")

## The KO comparator's gene overlap with the primary arm, read off the same frozen
## lists the sweep scored, so the "close to one observation twice" claim in
## why_included carries a number a reader can check.
mp_files <- unlist(UE$mouse_projection$files)
read_set <- function(p) unique(trimws(readLines(p, warn = FALSE)))
wt_f <- grep("WT_heat_up", mp_files, value = TRUE)[1]
ko_f <- grep("KO_heat_up", mp_files, value = TRUE)[1]
n_shared <- NA_integer_
if (length(wt_f) && length(ko_f) && file.exists(wt_f) && file.exists(ko_f)) {
  wt <- read_set(wt_f); ko <- read_set(ko_f)
  n_shared <- length(intersect(wt, ko))
  message(sprintf("[named_sets] WT_heat_up %d genes, KO_heat_up %d genes, %d shared",
                  length(wt), length(ko), n_shared))
} else {
  warning("[named_sets] mouse projection signature files unreachable; gene overlap left NA")
}
named <- named |>
  dplyr::mutate(n_genes_shared_with_WT_heat_up =
                  ifelse(pathway_id == "KO_heat_up", n_shared, NA_integer_)) |>
  dplyr::arrange(thread_order, dplyr::desc(included), pathway_id)

readr::write_csv(named, file.path(TDIR, "sweep_named_sets.csv"))
message(sprintf("  sweep_named_sets.csv: %d rows (%d drawn, %d excluded)",
                nrow(named), sum(named$included), sum(!named$included)))

# =============================================================================
# 3. Where each named set sits inside its own population's whole sweep
# =============================================================================
## Two orderings are computed and kept SEPARATE, because they answer different
## questions and mixing them in one sentence misreads both.
##   rank_padj_pooled  ascending pooled FDR, ties broken by nominal p then |NES|.
##       Benjamini-Hochberg produces exact FDR ties (HALLMARK_HYPOXIA ties
##       HALLMARK_APOPTOSIS in Treg), so an explicit tiebreak is what makes this
##       rank reproducible at all.
##   rank_nes_signed   descending signed NES, so rank 1 is the most
##       synovial-fluid-shifted set in that population.
## `rank_nes_abs` is carried too because a ranking on |NES| answers "largest shift in
## either direction" and lands a positive set much lower once the paired-blood-side
## translation and ribosome sets are counted.
ranked <- sweep |>
  dplyr::group_by(population) |>
  dplyr::mutate(
    n_tests_pooled_check = dplyr::n(),
    n_sig_pooled = sum(padj_pooled < FDR),
    ## order(order(v)) is the position each row takes in the sort defined by v, which is
    ## what a "rank under this exact ordering" means once ties are broken explicitly.
    rank_padj_pooled = order(order(padj_pooled, pvalue, -abs(nes))),
    rank_nes_signed  = order(order(-nes)),
    rank_nes_abs     = order(order(-abs(nes)))) |>
  dplyr::ungroup()
stopifnot("the pooled test count in the file must equal the rows present" =
            all(ranked$n_tests_pooled_check == ranked$n_tests_pooled))

## The grid is completed on purpose: every drawn set gets a row in every population,
## and a cell the sweep never scored carries tested = FALSE with NA statistics rather
## than vanishing. A missing row and a null result look identical on a panel, and only
## one of them is a result.
totals <- ranked |> dplyr::group_by(population) |>
  dplyr::summarise(n_tests_pooled = dplyr::first(n_tests_pooled),
                   n_sig_pooled = dplyr::first(n_sig_pooled), .groups = "drop")

stats <- tidyr::expand_grid(
    population = POP_LEVELS,
    selection |> dplyr::select(pathway_id, display_label, thread, thread_order)) |>
  dplyr::left_join(
    ranked |> dplyr::select(population, pathway_id, database, direction, nes, pvalue,
                            padj, padj_pooled, set_size, leading_edge_size,
                            rank_padj_pooled, rank_nes_signed, rank_nes_abs),
    by = c("population", "pathway_id")) |>
  dplyr::left_join(totals, by = "population") |>
  dplyr::transmute(
    population = factor(population, levels = POP_LEVELS), thread, thread_order,
    pathway_id, display_label, database, direction, nes, pvalue, padj, padj_pooled,
    tested = !is.na(nes),
    significant_pooled = !is.na(padj_pooled) & padj_pooled < FDR,
    set_size, leading_edge_size,
    rank_padj_pooled = as.integer(rank_padj_pooled),
    rank_nes_signed = as.integer(rank_nes_signed),
    rank_nes_abs = as.integer(rank_nes_abs),
    n_tests_pooled, n_sig_pooled) |>
  dplyr::arrange(population, thread_order, dplyr::desc(tested), padj_pooled)

stopifnot("every drawn set needs a row in every population" =
            nrow(stats) == nrow(selection) * length(POP_LEVELS))
readr::write_csv(round_numeric_cols(stats), file.path(TDIR, "sweep_named_sets_stats.csv"))
message(sprintf("  sweep_named_sets_stats.csv: %d rows", nrow(stats)))

# =============================================================================
# 4. Set-size baseline: how often a set of a given size reaches significance
# =============================================================================
## The hypoxia-versus-cGAS-STING comparison is confounded with how many genes MSigDB
## assigns to each label, and a figure that leaves this out reads as a biological
## verdict when it is substantially an artifact of set size. Two of the bands below
## are chosen to be the ones the figure's own points fall in: 10 to 22 spans all six
## cGAS-STING sets bar the 5-gene positive-regulation term, and 130 to 150 spans
## HALLMARK_HYPOXIA and both mouse arms. The rest give the trend across the range.
BAND_LO <- c(5, 10, 23, 50, 100, 130, 151, 250)
BAND_HI <- c(9, 22, 49, 99, 129, 150, 249, 500)
band_of <- function(n) {
  i <- vapply(n, function(x) {
    h <- which(x >= BAND_LO & x <= BAND_HI)
    if (length(h)) h[1] else NA_integer_
  }, integer(1))
  i
}

baseline <- sweep |>
  dplyr::mutate(band = band_of(set_size)) |>
  dplyr::filter(!is.na(band)) |>
  dplyr::group_by(population, band) |>
  dplyr::summarise(n_sets_tested = dplyr::n(),
                   n_pooled_significant = sum(padj_pooled < FDR),
                   .groups = "drop") |>
  dplyr::mutate(population = factor(population, levels = POP_LEVELS),
                size_low = BAND_LO[band], size_high = BAND_HI[band],
                band_label = sprintf("%d to %d genes", size_low, size_high),
                frac_pooled_significant = n_pooled_significant / n_sets_tested) |>
  dplyr::select(population, band_label, size_low, size_high, n_sets_tested,
                n_pooled_significant, frac_pooled_significant) |>
  dplyr::arrange(population, size_low)

## The two bands the caption cites, echoed to the log so a run that moved them is
## visible without opening the file.
for (b in c("10 to 22 genes", "130 to 150 genes")) {
  r <- baseline |> dplyr::filter(population == "Treg", band_label == b)
  message(sprintf("  Treg %s: %d of %d reach pooled FDR < %.2g (%.1f%%)",
                  b, r$n_pooled_significant[1], r$n_sets_tested[1], FDR,
                  100 * r$frac_pooled_significant[1]))
}
readr::write_csv(round_numeric_cols(baseline), file.path(TDIR, "sweep_setsize_baseline.csv"))
message(sprintf("  sweep_setsize_baseline.csv: %d rows", nrow(baseline)))

## Per-population totals, echoed for the same reason: "hypoxia comes out" means less
## when a large share of everything tested comes out.
for (p in POP_LEVELS) {
  s <- sweep |> dplyr::filter(population == p)
  message(sprintf("  %s: %d of %d tests reach pooled FDR < %.2g (%.1f%%)",
                  p, sum(s$padj_pooled < FDR), nrow(s), FDR,
                  100 * mean(s$padj_pooled < FDR)))
}
message("[named_sets] done")
