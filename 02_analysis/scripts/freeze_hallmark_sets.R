#!/usr/bin/env Rscript
## freeze_hallmark_sets.R — freeze MSigDB Hallmark gene sets to the project txt idiom
## ===========================================================================
## Phase 0 reproducer. Pulls MSigDB Hallmark ("H", Homo sapiens) via msigdbr and
## freezes two sets — HALLMARK_HYPOXIA and HALLMARK_UNFOLDED_PROTEIN_RESPONSE — as
## one-symbol-per-line .txt files plus a manifest CSV under the committed-reference
## convention 00_data/references/msigdb_hallmark/ (same idiom as etreg_GSE161426/).
##
## The symbols are HGNC gene symbols in the `gene_symbol` column of the msigdbr
## v26.1.0 tidy table (new API: msigdbr(species=, collection=)). Symbols are written
## SORTED (deterministic, order-independent) so the frozen list is reproducible.
##
## Downstream: adata objects carry symbols in var['gene_symbol']; these lists route
## directly into score_cells_aucell_ucell for per-cell scoring.
##
## Run in-container from the compartment root or anywhere:
##   Rscript 02_analysis/scripts/freeze_hallmark_sets.R
## ===========================================================================

suppressPackageStartupMessages(library(msigdbr))

## --- resolve compartment root relative to THIS script (self-contained per repo) ---
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg)) normalizePath(file_arg) else normalizePath("02_analysis/scripts/freeze_hallmark_sets.R")
PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), "..", ".."))
OUT_DIR <- file.path(PROJECT_ROOT, "00_data", "references", "msigdb_hallmark")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MSIGDBR_VERSION <- as.character(packageVersion("msigdbr"))
SETS <- c("HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")

## --- pull Hallmark (H) for human; gene symbols live in `gene_symbol` ---
h <- msigdbr(species = "Homo sapiens", collection = "H")

manifest <- data.frame(
  set_name = character(), n_genes = integer(), collection = character(),
  species = character(), msigdbr_version = character(), date_frozen = character(),
  stringsAsFactors = FALSE
)

for (set_name in SETS) {
  genes <- sort(unique(h$gene_symbol[h$gs_name == set_name]))
  txt_path <- file.path(OUT_DIR, paste0(set_name, ".txt"))
  writeLines(genes, txt_path)
  cat(sprintf("%s: %d genes -> %s\n", set_name, length(genes), txt_path))
  manifest <- rbind(manifest, data.frame(
    set_name = set_name, n_genes = length(genes), collection = "H",
    species = "Homo sapiens", msigdbr_version = MSIGDBR_VERSION,
    date_frozen = "FROZEN", stringsAsFactors = FALSE
  ))
}

manifest_path <- file.path(OUT_DIR, "hallmark_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
cat(sprintf("manifest -> %s (msigdbr %s)\n", manifest_path, MSIGDBR_VERSION))
