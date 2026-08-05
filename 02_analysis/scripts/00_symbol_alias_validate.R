#!/usr/bin/env Rscript
# 00_symbol_alias_validate.R: VALIDATION of the committed symbol-alias map.
# =============================================================================
# What this answers. The map built by 00_symbol_alias_map.R changes the effective size of
# gene sets this compartment has already published numbers for, so before any stage is
# re-run the map is checked against the ranked lists on disk: how many genes each named
# set recovers, which pairs did it, and whether the recovery agrees with the independent
# audit of 2026-08-05. A disagreement is REPORTED here rather than reconciled — the audit
# and this script are two reads of the same disk and either can be wrong, so the
# disagreement is the finding.
#
# It also runs the migration test for the lift. 18_tf_activity.R carried a private copy of
# the alias machinery and published two ledger tables from it; the shared helper must
# reproduce both. That check is the cheapest proof the lift changed nothing, and it runs
# here rather than by re-running that stage, whose expensive random-regulon nulls have
# nothing to do with alias resolution.
#
# NOTHING PUBLISHED IS TOUCHED. Every output lands under 03_results/_scratch/, and no
# stage table or figure is regenerated. The alias fix moves set sizes in stages 05 and 09
# through 16, and the frozen mouse contract those stages consume is itself being corrected
# upstream, so regenerating them now would only mean regenerating them again.
#
# Reads, read-only:
#   00_data/references/symbol_alias/symbol_alias_map.csv       the committed map
#   03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv
#   03_results/03_pseudobulk/tables/gene_symbols.csv           the post-QC vocabulary
#   03_results/00_build/tables/reference_feature_symbols.csv   the CellRanger feature union
#   03_results/18_tf_activity/tables/{alias_recovery,symbol_vocabulary_check}.csv
#   the five STING families from msigdbr, the frozen lenses, the projected mouse arms
#
# Writes 03_results/_scratch/symbol_alias/:
#   geneset_symbol_ledger.csv    the full per-set x population ledger, all seven buckets
#   audit_reproduction.csv       expected versus reproduced, with a mismatch column
#   alias_pairs_applied.csv      every pair applied, per set and population
#   migration_18_tf_activity.csv the lift's regression result
#
# Run from the compartment root, after 00_symbol_alias_map.R:
#   Rscript 02_analysis/scripts/00_symbol_alias_validate.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")   # FIG_CFG
source("02_analysis/helpers/symbol_alias.R")   # build_alias_map, resolve_sets, symbol_ledger

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

CFG <- FIG_CFG
SA  <- CFG$symbol_alias
UE  <- CFG$unbiased_enrichment
OUT <- file.path(CFG$paths$scratch %||% "03_results/_scratch/", "symbol_alias")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
POPS <- c(treg = "Treg", tcon = "Tcon", cd8 = "CD8")

emit <- function(df, fname) {
  p <- file.path(OUT, fname)
  readr::write_csv(df, p)
  message(sprintf("  [SAVE] %-32s %5d rows x %2d cols", fname, nrow(df), ncol(df)))
  invisible(p)
}

message("=================================================================")
message("00_symbol_alias_validate — the map against the ranked lists on disk")
message("=================================================================")

# ============================================================================
# 1. The three nested vocabulary layers, and the map
# ============================================================================
# Each layer answers a different question about a missing gene, and the whole point of the
# ledger is that they are not the same question. Absent from layer 1 is a statement about
# the reference; present in 1 and absent from 3 is a statement about detection; present in
# 3 and absent from 4 is a statement about this contrast's power.

MAP <- readr::read_csv(SA$map_path, show_col_types = FALSE)
stopifnot(nrow(MAP) > 0, all(MAP$resolution %in% ALIAS_RESOLUTIONS))
MATRIX_SYMBOLS <- unique(readr::read_csv(SA$matrix_vocabulary,
                                         show_col_types = FALSE)$gene_symbol)
if (!file.exists(SA$reference_feature_symbols))
  stop("[00_val] ", SA$reference_feature_symbols, " is missing. Without the CellRanger ",
       "feature union, 'genuinely absent' cannot be told apart from 'never detected', ",
       "which is the defect this ledger exists to keep apart. Re-run 00_build_anndata.py.")
REF_FEATURES <- unique(readr::read_csv(SA$reference_feature_symbols,
                                       show_col_types = FALSE)$gene_symbol)

read_ranked <- function(tag) {
  p <- gsub("{population}", tag, SA$ranked_list, fixed = TRUE)
  d <- utils::read.table(p, sep = "\t", header = FALSE,
                         col.names = c("symbol", "stat"), stringsAsFactors = FALSE)
  d <- d[!duplicated(d$symbol), , drop = FALSE]
  stats::setNames(d$stat, d$symbol)
}
RANKED <- lapply(names(POPS), read_ranked); names(RANKED) <- names(POPS)

message(sprintf("[1] layers: reference union %d > post-QC matrix %d > ranked %s",
                length(REF_FEATURES), length(MATRIX_SYMBOLS),
                paste(sprintf("%s %d", POPS, lengths(RANKED)), collapse = " / ")))
message(sprintf("    map: %d accepted, %d flagged for review, %d rejected",
                sum(MAP$resolution == "accepted"),
                sum(MAP$resolution == "flagged_for_review"),
                sum(!MAP$resolution %in% c("accepted", "flagged_for_review"))))

# ============================================================================
# 2. The sets under test
# ============================================================================
# The five STING families are fetched live from msigdbr, because they are the sets whose
# published matched counts the audit itemises and the point is to reproduce those from the
# reference rather than from a cache.

read_txt <- function(p) {
  g <- trimws(readLines(p, warn = FALSE)); unique(g[nzchar(g)])
}
SIG_DIR <- CFG$paths$signature_contract %||% "../mouse_anchor/03_results/human_projection/"

STING_IDS <- c(
  "GOBP_CGAS_STING_SIGNALING_PATHWAY",
  "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",
  "WP_STING_PATHWAY_IN_KAWASAKILIKE_DISEASE_AND_COVID19",
  "GOBP_NEGATIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY")

message("[2] fetching the five STING families from msigdbr ",
        as.character(packageVersion("msigdbr")), " ...")
msig <- msigdbr::msigdbr(species = CFG$project$species %||% "Homo sapiens")
sym_col <- if ("gene_symbol" %in% colnames(msig)) "gene_symbol" else "human_gene_symbol"
all_msig <- lapply(split(as.character(msig[[sym_col]]), as.character(msig$gs_name)), unique)
missing_ids <- setdiff(STING_IDS, names(all_msig))
if (length(missing_ids))
  stop("[00_val] MSigDB no longer carries: ", paste(missing_ids, collapse = ", "))
SETS <- all_msig[STING_IDS]

arm_files <- c(
  file.path(SIG_DIR, "signatures/WT_heat", c("WT_heat_up.txt", "WT_heat_down.txt")),
  file.path(SIG_DIR, "signatures/KO_heat", c("KO_heat_up.txt", "KO_heat_down.txt")),
  file.path(SIG_DIR, "signatures/Interaction",
            c("Interaction_up.txt", "Interaction_fdrOnly_up.txt")))
lens_files <- c(
  unlist(UE$project_frozen$files), unlist(UE$hsr_lens$files),
  unlist(UE$tcr_activation$files), unlist(UE$sting_axes$files),
  file.path(CFG$paths$references %||% "00_data/references/", "etreg_GSE161426",
            c("eTreg_up.txt", "eTreg_down.txt")))
for (p in unique(c(arm_files, lens_files))) {
  if (!file.exists(p)) { message("    absent, skipped: ", p); next }
  SETS[[sub("\\.txt$", "", basename(p))]] <- read_txt(p)
}
message(sprintf("[2] %d sets under test", length(SETS)))

# ============================================================================
# 3. The ledger, per population
# ============================================================================

LEDGER <- bind_rows(lapply(names(POPS), function(tag) {
  g <- names(RANKED[[tag]])
  led <- symbol_ledger(SETS, MAP, ranked_vocabulary = g,
                       matrix_vocabulary = MATRIX_SYMBOLS,
                       reference_vocabulary = REF_FEATURES)
  assert_ledger_closes(led, sprintf("ledger for %s", POPS[[tag]]))
  res <- resolve_sets(SETS, g, MAP)
  # A resolved set must never lose a gene and must never gain one twice.
  stopifnot(all(vapply(res$sets, function(s) !any(duplicated(s)), logical(1))))
  led$population <- POPS[[tag]]
  led$set_size_resolved_recomputed <- vapply(
    names(SETS), function(nm) length(intersect(res$sets[[nm]], g)), integer(1))
  led
})) %>% relocate(population)

# The map may only ever ADD. If resolve_sets() and the ledger disagree on the resolved
# size, one of the two is applying a pair the other is not.
stopifnot("resolve_sets and symbol_ledger must agree on every resolved set size" =
            all(LEDGER$set_size_resolved == LEDGER$set_size_resolved_recomputed))
emit(LEDGER %>% select(-set_size_resolved_recomputed), "geneset_symbol_ledger.csv")

PAIRS <- bind_rows(lapply(names(POPS), function(tag) {
  r <- resolve_sets(SETS, names(RANKED[[tag]]), MAP)
  if (!nrow(r$applied)) return(NULL)
  r$applied %>% mutate(population = POPS[[tag]]) %>% relocate(population)
}))
emit(PAIRS, "alias_pairs_applied.csv")

collapses <- bind_rows(lapply(names(POPS), function(tag) {
  r <- resolve_sets(SETS, names(RANKED[[tag]]), MAP)
  bind_rows(
    if (nrow(r$collapsed)) r$collapsed %>% mutate(population = POPS[[tag]],
                                                  issue = "duplicate_collapse"),
    if (nrow(r$many_to_one)) r$many_to_one %>% mutate(population = POPS[[tag]],
                                                      issue = "many_to_one"))
}))
if (nrow(collapses)) {
  message("[3] duplicate collapses and many-to-one resolutions, reported not merged:")
  print(as.data.frame(collapses), row.names = FALSE)
} else message("[3] no duplicate collapse and no many-to-one resolution in any set.")

# ============================================================================
# 4. Reproduction against the 2026-08-05 audit
# ============================================================================
# The audited numbers, transcribed once, as the check to reproduce. They were measured
# against the ranked lists on disk BEFORE the mouse-side ortholog fix lands, so once the
# projected arms change, the arm rows here move for that reason and this table is what
# says which rows moved and why.
#
# `expected_matched` is the published exact-match count and MUST NOT move: alias
# resolution only ever adds. `expected_alias` is what the audit says it adds.

AUDIT <- tribble(
  ~gene_set,                                                     ~population, ~expected_matched, ~expected_alias,
  "GOBP_CGAS_STING_SIGNALING_PATHWAY",                           "Treg", 22L, 3L,
  "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES",  "Treg", 10L, 3L,
  "WP_STING_PATHWAY_IN_KAWASAKILIKE_DISEASE_AND_COVID19",        "Treg", 16L, 2L,
  "GOBP_NEGATIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY",    "Treg", 13L, 1L,
  "GOBP_POSITIVE_REGULATION_OF_CGAS_STING_SIGNALING_PATHWAY",    "Treg",  5L, 1L,
  "WT_heat_up",         "Treg", 119L, 0L,  "WT_heat_up",         "Tcon", 130L, 0L,  "WT_heat_up",         "CD8", 113L, 0L,
  "WT_heat_down",       "Treg",  56L, 2L,  "WT_heat_down",       "Tcon",  61L, 2L,  "WT_heat_down",       "CD8",  57L, 2L,
  "KO_heat_up",         "Treg", 132L, 3L,  "KO_heat_up",         "Tcon", 143L, 3L,  "KO_heat_up",         "CD8", 125L, 3L,
  "KO_heat_down",       "Treg",  66L, 2L,  "KO_heat_down",       "Tcon",  72L, 2L,  "KO_heat_down",       "CD8",  70L, 2L,
  "Interaction_up",     "Treg",   6L, 0L,
  "Interaction_fdrOnly_up", "Treg", 17L, 0L,
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",   "Treg", 103L, 7L,
  "HALLMARK_HYPOXIA",                     "Treg", 139L, 4L,
  "HALLMARK_IL2_STAT5_SIGNALING",         "Treg", 167L, 4L,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",   "Treg",  91L, 2L,
  "HALLMARK_INFLAMMATORY_RESPONSE",       "Treg", 141L, 2L,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",     "Treg", 175L, 1L,
  "HSR_core",                             "Treg",  43L, 1L,
  "HSR_sensitivity",                      "Treg", 137L, 3L,
  "TCR_activation",                       "Treg",  63L, 1L,
  # The audit reports +2 here. One of the two is MIR4435-2HG->MIR4435-1HG, withheld for
  # review, so the expectation under the exclusion is +1 and 13 rather than 14. A result
  # of +2 would mean the exclusion is not being applied.
  "sting_specific_up",                    "Treg",  12L, 1L,
  "ifn_only_up",                          "Treg",  60L, 1L,
  "eTreg_up",                             "Treg", 167L, 0L,
  "eTreg_down",                           "Treg", 146L, 0L)

REPRO <- AUDIT %>%
  left_join(LEDGER %>% select(population, gene_set, n_unique_set_genes, n_matched,
                              n_matched_via_alias, n_alias_flagged_for_review,
                              set_size_resolved, alias_pairs_applied, alias_pairs_flagged),
            by = c("gene_set", "population")) %>%
  mutate(expected_new = .data$expected_matched + .data$expected_alias,
         matched_mismatch = .data$n_matched != .data$expected_matched,
         alias_mismatch = .data$n_matched_via_alias != .data$expected_alias,
         mismatch = .data$matched_mismatch | .data$alias_mismatch) %>%
  select(population, gene_set, n_unique_set_genes, expected_matched, n_matched,
         matched_mismatch, expected_alias, n_matched_via_alias, alias_mismatch,
         expected_new, set_size_resolved, n_alias_flagged_for_review, mismatch,
         alias_pairs_applied, alias_pairs_flagged)
emit(REPRO, "audit_reproduction.csv")

message("\n[4] reproduction against the audit:")
print(as.data.frame(REPRO %>% select(population, gene_set, expected_matched, n_matched,
                                     expected_alias, n_matched_via_alias, expected_new,
                                     set_size_resolved, mismatch)), row.names = FALSE)
n_bad <- sum(REPRO$mismatch, na.rm = TRUE)
if (n_bad == 0) {
  message(sprintf("[4] all %d audited rows reproduce exactly.", nrow(REPRO)))
} else {
  message(sprintf("[4] %d of %d audited rows DISAGREE — reported, not reconciled:",
                  n_bad, nrow(REPRO)))
  print(as.data.frame(REPRO %>% filter(.data$mismatch)), row.names = FALSE)
}

# ============================================================================
# 5. The migration test: the lift must not have changed behaviour
# ============================================================================
# 18_tf_activity.R published its two ledger tables from a private copy of this machinery.
# The shared helper is re-run here on that stage's own inputs and the two tables are
# rebuilt and compared to what is committed. Content must be identical apart from the
# deliberate network_symbol -> reference_symbol rename and the new explicit rejection
# classes, which the private copy returned silently as NA.
#
# This mirrors sections 2 and 3 of that stage; the expensive random-regulon nulls have
# nothing to do with alias resolution, so the stage itself is not re-run.

message("\n[5] migration test against the committed 18_tf_activity tables ...")
TF_TBL   <- "03_results/18_tf_activity/tables"
NET_PATH <- UE$tf_network$path
TA       <- CFG$tf_activity
FOCUS    <- unlist(TA$focus_tfs)

net_raw <- readr::read_csv(NET_PATH, show_col_types = FALSE, progress = FALSE)
net_base <- net_raw %>%
  transmute(source, target, mor = as.numeric(.data[[TA$mor_col %||% "weight"]]),
            sign_decision = .data[[TA$sign_decision_col %||% "sign_decision"]]) %>%
  distinct(source, target, .keep_all = TRUE)

mig_maps <- lapply(names(POPS), function(tag) {
  universe <- names(RANKED[[tag]])
  build_alias_map(sort(setdiff(unique(net_base$target), universe)), universe,
                  db = org.Hs.eg.db::org.Hs.eg.db)
})
names(mig_maps) <- names(POPS)

new_recovery <- bind_rows(lapply(names(POPS), function(tag) {
  amap <- mig_maps[[tag]]
  acc <- amap %>% filter(.data$resolution == "accepted")
  rej <- amap %>% filter(.data$resolution == "rejected_symbol_belongs_to_another_gene")
  bind_rows(
    net_base %>% filter(!.data$target %in% names(RANKED[[tag]])) %>%
      inner_join(acc, by = c("target" = "reference_symbol")) %>%
      transmute(population = POPS[[tag]], tf = .data$source, reference_symbol = .data$target,
                matrix_symbol = .data$matrix_symbol, .data$mor,
                resolution = "accepted", focus_tf = .data$source %in% FOCUS),
    net_base %>% filter(.data$target %in% rej$reference_symbol) %>%
      transmute(population = POPS[[tag]], tf = .data$source, reference_symbol = .data$target,
                matrix_symbol = rej$matrix_symbol[match(.data$target, rej$reference_symbol)],
                .data$mor, resolution = "rejected_symbol_belongs_to_another_gene",
                focus_tf = .data$source %in% FOCUS))
}))

old_recovery <- readr::read_csv(file.path(TF_TBL, "alias_recovery.csv"),
                                show_col_types = FALSE) %>%
  rename(reference_symbol = network_symbol)
key <- function(d) sort(with(d, paste(population, tf, reference_symbol, matrix_symbol,
                                      mor, resolution, focus_tf)))
rec_identical <- identical(key(old_recovery), key(new_recovery))

new_vocab <- bind_rows(lapply(names(POPS), function(tag) {
  universe <- names(RANKED[[tag]])
  acc <- mig_maps[[tag]] %>% filter(.data$resolution == "accepted")
  bind_rows(lapply(FOCUS, function(tf) {
    t <- unique(net_base$target[net_base$source == tf])
    unmatched <- setdiff(t, universe)
    rec <- intersect(unmatched, acc$reference_symbol)
    tibble(population = POPS[[tag]], tf = tf,
           n_targets_in_network = length(t),
           n_matched = length(intersect(t, universe)),
           n_unmatched = length(unmatched),
           n_expression_filtered = length(intersect(unmatched, MATRIX_SYMBOLS)),
           n_absent_from_count_matrix = length(setdiff(unmatched, MATRIX_SYMBOLS)),
           n_alias_recoverable = length(rec),
           alias_recoverable_symbols = paste(sort(rec), collapse = "/"))
  }))
}))
old_vocab <- readr::read_csv(file.path(TF_TBL, "symbol_vocabulary_check.csv"),
                             show_col_types = FALSE)
vocab_cols <- intersect(names(old_vocab), names(new_vocab))
vocab_identical <- isTRUE(all.equal(
  as.data.frame(old_vocab[order(old_vocab$population, old_vocab$tf), vocab_cols]),
  as.data.frame(new_vocab[order(new_vocab$population, new_vocab$tf), vocab_cols]),
  check.attributes = FALSE))

MIG <- tibble(
  table = c("alias_recovery.csv", "symbol_vocabulary_check.csv"),
  n_rows_committed = c(nrow(old_recovery), nrow(old_vocab)),
  n_rows_shared_helper = c(nrow(new_recovery), nrow(new_vocab)),
  n_accepted_committed = c(sum(old_recovery$resolution == "accepted"), NA_integer_),
  n_accepted_shared_helper = c(sum(new_recovery$resolution == "accepted"), NA_integer_),
  content_identical = c(rec_identical, vocab_identical),
  compared_on = c("population/tf/reference_symbol/matrix_symbol/mor/resolution/focus_tf",
                  paste(vocab_cols, collapse = ",")))
emit(MIG, "migration_18_tf_activity.csv")
print(as.data.frame(MIG %>% select(-compared_on)), row.names = FALSE)
if (all(MIG$content_identical)) {
  message("[5] migration test PASSES: the shared helper reproduces both committed tables.")
} else {
  message("[5] migration test FAILS — the lift changed behaviour. Differences:")
  if (!rec_identical) {
    print(head(setdiff(key(new_recovery), key(old_recovery)), 20))
    print(head(setdiff(key(old_recovery), key(new_recovery)), 20))
  }
  if (!vocab_identical)
    print(all.equal(as.data.frame(old_vocab[order(old_vocab$population, old_vocab$tf),
                                            vocab_cols]),
                    as.data.frame(new_vocab[order(new_vocab$population, new_vocab$tf),
                                            vocab_cols]), check.attributes = FALSE))
}

message(sprintf("\n[DONE] validation complete: %d audited rows, %d disagreements, migration %s.",
                nrow(REPRO), n_bad, if (all(MIG$content_identical)) "passes" else "FAILS"))
