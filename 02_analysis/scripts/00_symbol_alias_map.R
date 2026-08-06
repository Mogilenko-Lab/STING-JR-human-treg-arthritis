#!/usr/bin/env Rscript
# 00_symbol_alias_map.R: ASSET BUILDER. The committed reference-to-matrix symbol map.
# =============================================================================
# GSE160097 was quantified against a CellRanger hg19 reference, so this compartment's
# count matrix is frozen to that build's HGNC vintage — cGAS is MB21D1, STING is
# TMEM173, MARCHF5 is MARCH5, MRE11 is MRE11A. Reference gene sets ship current symbols
# and every match in this compartment is an exact string match, so a renamed gene leaves
# a set silently and the loss reads as biological absence. This script resolves that once
# and writes the answer down.
#
# WHY A COMMITTED CSV. The map is a property of (matrix vocabulary x org.Hs.eg.db release)
# and of nothing else, so it is computed once here. Python cannot reach org.Hs.eg.db, and a
# subprocess per call would load the database each time and leave nothing to inspect.
# Hardcoding the pairs in config would go stale without a diff. A committed CSV is auditable,
# is read identically by R and Python, and moves when someone regenerates it on purpose.
#
# Rebuild when the matrix vocabulary changes, when org.Hs.eg.db moves, or when a
# collection is added to `unbiased_enrichment:`. The provenance file records what the map was
# built from, which makes a stale map visible.
#
# Reads, read-only:
#   02_analysis/config/analysis_config.yaml               symbol_alias + unbiased_enrichment
#   03_results/03_pseudobulk/tables/gene_symbols.csv      the vocabulary resolved INTO
#   the 11 declared reference collections, plus every frozen *.txt gene list this
#   compartment reads (the projected mouse arms in both directions, the curated lenses,
#   the reference interferon axes, the eTreg sets)
#
# Writes (committed assets, NOT stage results):
#   00_data/references/symbol_alias/symbol_alias_map.csv
#   00_data/references/symbol_alias/symbol_alias_provenance.csv
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/00_symbol_alias_map.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

source("02_analysis/helpers/figure_style.R")          # FIG_CFG
source("02_analysis/helpers/source_hash_manifest.R")  # source_sha256
source("02_analysis/helpers/symbol_alias.R")          # build_alias_map

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

CFG <- FIG_CFG
SA  <- CFG$symbol_alias
UE  <- CFG$unbiased_enrichment
if (is.null(SA))
  stop("[00_alias] analysis_config.yaml has no `symbol_alias:` block — add it before running.")
if (is.null(UE))
  stop("[00_alias] analysis_config.yaml has no `unbiased_enrichment:` block.")

SPECIES  <- CFG$project$species %||% "Homo sapiens"
MAP_OUT  <- SA$map_path
PROV_OUT <- SA$provenance_path
VOCAB_CSV <- SA$matrix_vocabulary
FLAGGED  <- unlist(SA$flagged_for_review) %||% character()
EXTRA_DIRS <- unlist(SA$extra_reference_dirs) %||% character()
dir.create(dirname(MAP_OUT), recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("00_symbol_alias_map — reference symbols -> this matrix's vocabulary")
message("=================================================================")

# ============================================================================
# 1. The vocabulary the map resolves INTO
# ============================================================================
# gene_symbols.csv is the seam map 03a_pseudobulk_export.py emits, so it holds every
# symbol the count matrix carried BEFORE filterByExpr. Building the map against it rather
# than against a ranked list keeps the map a property of the matrix: a pair whose target
# is later dropped for low expression is then reported as expression-filtered, which is a
# power statement about a contrast.

stopifnot(file.exists(VOCAB_CSV))
MATRIX_SYMBOLS <- unique(readr::read_csv(VOCAB_CSV, show_col_types = FALSE)$gene_symbol)
MATRIX_SYMBOLS <- MATRIX_SYMBOLS[!is.na(MATRIX_SYMBOLS) & nzchar(MATRIX_SYMBOLS)]
message(sprintf("[1] matrix vocabulary: %d symbols from %s",
                length(MATRIX_SYMBOLS), VOCAB_CSV))
for (p in c("MB21D1", "TMEM173", "MARCH5", "MRE11A", "CGAS", "STING1"))
  message(sprintf("    %-8s in matrix vocabulary: %s", p, p %in% MATRIX_SYMBOLS))

# ============================================================================
# 2. Every reference symbol this compartment consumes
# ============================================================================
# The union is deliberately wider than any single stage needs, because the map has to
# answer for every consumer without being rebuilt per stage.

sources <- list()
note <- function(label, syms) {
  syms <- unique(syms[!is.na(syms) & nzchar(syms)])
  sources[[label]] <<- syms
  message(sprintf("    %-28s %6d unique symbols", label, length(syms)))
  invisible(syms)
}

message("[2] collecting reference symbols ...")

# --- MSigDB collections, fetched exactly as the sweep fetches them --------------------
# 14_unbiased_enrichment.R owns the canonical loader; this is the same config block and
# the same msigdbr call, kept minimal because only the SYMBOLS matter here. msigdbr 26.1.0
# renamed CP:KEGG to CP:KEGG_LEGACY, and the retry is what lets the config keep the
# mouse-anchor spelling.
fetch_msigdb <- function(category, subcategory = "") {
  has_new <- "collection" %in% names(formals(msigdbr::msigdbr))
  go <- function(sub_use) {
    if (has_new)
      msigdbr::msigdbr(species = SPECIES, collection = category,
                       subcollection = if (nzchar(sub_use)) sub_use else NULL)
    else
      msigdbr::msigdbr(species = SPECIES, category = category,
                       subcategory = if (nzchar(sub_use)) sub_use else NULL)
  }
  df <- tryCatch(go(subcategory), error = function(e) {
    alt <- sub("CP:KEGG$", "CP:KEGG_LEGACY", subcategory)
    if (nzchar(subcategory) && alt != subcategory &&
        grepl("[Uu]nknown subcollection", conditionMessage(e))) return(go(alt))
    stop("[00_alias] MSigDB fetch failed for ", category, "/", subcategory, ": ",
         conditionMessage(e))
  })
  sym <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  as.character(df[[sym]])
}
for (m in UE$msigdb)
  note(paste0("msigdb:", m$name), fetch_msigdb(m$category, m$subcategory %||% ""))

# --- the CollecTRI regulon targets ----------------------------------------------------
if (!is.null(UE$tf_network) && file.exists(UE$tf_network$path)) {
  d <- readr::read_csv(UE$tf_network$path, show_col_types = FALSE, progress = FALSE)
  note("collectri:targets", as.character(d[[UE$tf_network$target_col %||% "target"]]))
} else message("    CollecTRI network absent; its targets are not covered by this map.")

# --- toolkit reference databases (MitoPathways) ---------------------------------------
for (spec in UE$custom_rds %||% list()) {
  if (!file.exists(spec$path)) { message("    missing custom_rds: ", spec$path); next }
  obj <- readRDS(spec$path)
  note(paste0("custom_rds:", spec$name), as.character(obj$T2G[[2]]))
}

# --- every frozen one-symbol-per-line list, from config and from the asset trees -------
# Config names the files the sweep scores; the extra directories catch the ones only the
# Python stages read — the projected DOWN arms above all, which no sweep collection
# declares and which stages 05, 09 and 13 all consume.
read_txt <- function(p) {
  g <- trimws(readLines(p, warn = FALSE))
  unique(g[nzchar(g)])
}
FILE_KEYS <- c("project_frozen", "mouse_projection", "sting_axes", "hsr_lens",
               "tcr_activation")
declared <- unlist(lapply(FILE_KEYS, function(k) unlist(UE[[k]]$files %||% character())))
found <- unlist(lapply(EXTRA_DIRS, function(d)
  if (dir.exists(d)) list.files(d, pattern = "\\.txt$", recursive = TRUE,
                                full.names = TRUE) else character()))
txt_files <- unique(c(declared, found))
txt_files <- txt_files[file.exists(txt_files)]
note("frozen_gene_lists", unlist(lapply(txt_files, read_txt)))
message(sprintf("    (%d frozen list files: %d declared in config, %d found under %s)",
                length(txt_files), sum(declared %in% txt_files), length(found),
                paste(EXTRA_DIRS, collapse = " ")))

REFERENCE_SYMBOLS <- sort(unique(unlist(sources, use.names = FALSE)))
message(sprintf("[2] reference universe: %d unique symbols over %d sources",
                length(REFERENCE_SYMBOLS), length(sources)))

# ============================================================================
# 3. Build the map
# ============================================================================

message("[3] resolving against org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")),
        " ...")
MAP <- build_alias_map(REFERENCE_SYMBOLS, MATRIX_SYMBOLS,
                       db = org.Hs.eg.db::org.Hs.eg.db, flagged_pairs = FLAGGED)
SUMM <- attr(MAP, "summary")

# Every flagged pair must be a pair the guard would otherwise have accepted, or the
# exclusion list is stale and reads as a decision that is no longer doing anything.
missed_flags <- setdiff(FLAGGED, paste0(MAP$reference_symbol, "->", MAP$matrix_symbol))
if (length(missed_flags))
  stop(sprintf(paste0("[00_alias] flagged_for_review names %d pair(s) this build does not ",
                      "produce as a candidate at all (%s). Either the pair is stale or the ",
                      "reference universe no longer contains it — decide, do not leave it ",
                      "silently inert."),
               length(missed_flags), paste(missed_flags, collapse = ", ")))

MAP <- MAP %>% arrange(match(.data$resolution, ALIAS_RESOLUTIONS), .data$reference_symbol)
readr::write_csv(MAP, MAP_OUT)
message(sprintf("  [SAVE] %s  %d candidate pairs", MAP_OUT, nrow(MAP)))
print(as.data.frame(MAP %>% count(resolution, name = "n_pairs")), row.names = FALSE)

message("  the four named pre-2019 probes:")
print(as.data.frame(MAP %>% filter(.data$reference_symbol %in%
                                     c("CGAS", "STING1", "MARCHF5", "MRE11"))),
      row.names = FALSE)
message("  a sample of the pairs the ownership guard rejected:")
print(as.data.frame(MAP %>%
                      filter(.data$resolution == "rejected_symbol_belongs_to_another_gene") %>%
                      transmute(pair = paste0(.data$reference_symbol, "->",
                                              .data$matrix_symbol)) %>% head(20)),
      row.names = FALSE)

# ============================================================================
# 4. Provenance — what the map was built from, so a stale map is visible
# ============================================================================
# The no-candidate counts live here, outside the map's rows. A reference symbol
# with no alias in this vocabulary is not a decision about a pair, and carrying thousands
# of such rows would defeat the one review step this asset most needs: reading the map
# cold and recognising every accepted pair as a nomenclature update.

PROV <- tibble(
  built_by = "02_analysis/scripts/00_symbol_alias_map.R",
  org_hs_eg_db_version = as.character(packageVersion("org.Hs.eg.db")),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  species = SPECIES,
  matrix_vocabulary_path = VOCAB_CSV,
  matrix_vocabulary_sha256 = source_sha256(VOCAB_CSV),
  n_matrix_symbols = length(MATRIX_SYMBOLS),
  n_reference_symbols = SUMM$n_reference_symbols,
  n_reference_sources = length(sources),
  reference_sources = paste(names(sources), collapse = "; "),
  n_frozen_list_files = length(txt_files),
  n_reference_absent_from_vocabulary = SUMM$n_absent_from_vocabulary,
  n_reference_not_in_org_db = SUMM$n_not_in_org_db,
  n_reference_no_alias_in_vocabulary = SUMM$n_no_alias_in_vocabulary,
  n_candidate_pairs = SUMM$n_candidates,
  n_accepted = SUMM$n_accepted,
  n_flagged_for_review = SUMM$n_flagged,
  n_rejected = SUMM$n_rejected,
  flagged_for_review = paste(FLAGGED, collapse = "; "))
readr::write_csv(PROV, PROV_OUT)
message(sprintf("  [SAVE] %s", PROV_OUT))
print(as.data.frame(t(PROV)))

# The map may only ever ADD a symbol to a set, so nothing it accepts may already be a
# reference symbol that matched, and every accepted target must be in the vocabulary.
stopifnot(
  "an accepted target must be in the matrix vocabulary" =
    all(MAP$matrix_symbol[MAP$resolution == "accepted"] %in% MATRIX_SYMBOLS),
  "an accepted reference symbol must be ABSENT from the matrix vocabulary" =
    !any(MAP$reference_symbol[MAP$resolution == "accepted"] %in% MATRIX_SYMBOLS),
  "one reference symbol resolves to at most one matrix symbol" =
    !any(duplicated(MAP$reference_symbol[MAP$resolution == "accepted"])),
  "the four named probes must resolve" =
    all(c("CGAS", "STING1", "MARCHF5", "MRE11") %in%
          MAP$reference_symbol[MAP$resolution == "accepted"]))

message("\n[DONE] symbol alias map built. Validate with 00_symbol_alias_validate.R.")
