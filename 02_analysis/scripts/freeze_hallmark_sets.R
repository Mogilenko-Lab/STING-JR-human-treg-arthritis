#!/usr/bin/env Rscript
## freeze_hallmark_sets.R — freeze MSigDB Hallmark gene sets to the project txt idiom
## ===========================================================================
## Reference reproducer. Pulls MSigDB Hallmark ("H", Homo sapiens) via msigdbr and
## freezes the sets this compartment scores or decomposes against, as one-symbol-per-line
## .txt files plus a manifest CSV and a PROVENANCE.md under the committed-reference
## convention 00_data/references/msigdb_hallmark/ (same idiom as etreg_GSE161426/).
##
## The frozen sets and their roles:
##   HALLMARK_HYPOXIA                    the gene-purge reference for the mouse WT_heat sets
##   HALLMARK_UNFOLDED_PROTEIN_RESPONSE  ER-side proteostasis
##   HALLMARK_TNFA_SIGNALING_VIA_NFKB    inflammatory signalling arm
##   HALLMARK_INTERFERON_ALPHA_RESPONSE  type-I interferon arm (cGAS/STING-adjacent)
##   HALLMARK_INFLAMMATORY_RESPONSE      broad inflammation arm
##   HALLMARK_IL2_STAT5_SIGNALING        T-cell activation / growth-signalling arm
##
## The symbols are HGNC gene symbols in the `gene_symbol` column of the msigdbr
## v26.1.0 tidy table (new API: msigdbr(species=, collection=)). Symbols are written
## SORTED (deterministic, order-independent) so the frozen list is reproducible.
##
## Each set carries a VALIDATED expected size: the script stops when the installed msigdbr
## data drifts from the release these assets were frozen against, since a size change would
## move every downstream set-intersection tally silently.
##
## Downstream: adata objects carry symbols in var['gene_symbol']; these lists route
## directly into score_cells_aucell_ucell for per-cell scoring, and into the mouse
## WT_heat decomposition arms of 11_heat_decomposition.py.
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

## set_name -> validated n_genes for MSigDB v2026.1.Hs via msigdbr 26.1.0.
SETS <- c(
  HALLMARK_HYPOXIA                   = 200L,
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE = 113L,
  HALLMARK_TNFA_SIGNALING_VIA_NFKB   = 200L,
  HALLMARK_INTERFERON_ALPHA_RESPONSE =  97L,
  HALLMARK_INFLAMMATORY_RESPONSE     = 200L,
  HALLMARK_IL2_STAT5_SIGNALING       = 199L
)

## --- pull Hallmark (H) for human; gene symbols live in `gene_symbol` ---
h <- msigdbr(species = "Homo sapiens", collection = "H")

absent <- setdiff(names(SETS), unique(h$gs_name))
if (length(absent)) {
  stop(sprintf("Hallmark set(s) absent from msigdbr %s: %s",
               MSIGDBR_VERSION, paste(absent, collapse = ", ")), call. = FALSE)
}

manifest <- data.frame(
  set_name = character(), n_genes = integer(), n_genes_expected = integer(),
  collection = character(), species = character(), msigdbr_version = character(),
  date_frozen = character(), stringsAsFactors = FALSE
)

for (set_name in names(SETS)) {
  genes <- sort(unique(h$gene_symbol[h$gs_name == set_name]))
  expected <- SETS[[set_name]]
  if (length(genes) != expected) {
    stop(sprintf(
      "%s size drift: expected %d genes, found %d. The installed msigdbr data no longer matches the release these assets were validated against.",
      set_name, expected, length(genes)), call. = FALSE)
  }
  txt_path <- file.path(OUT_DIR, paste0(set_name, ".txt"))
  writeLines(genes, txt_path)
  cat(sprintf("%s: %d genes -> %s\n", set_name, length(genes), txt_path))
  manifest <- rbind(manifest, data.frame(
    set_name = set_name, n_genes = length(genes), n_genes_expected = expected,
    collection = "H", species = "Homo sapiens", msigdbr_version = MSIGDBR_VERSION,
    date_frozen = "FROZEN", stringsAsFactors = FALSE
  ))
}

manifest_path <- file.path(OUT_DIR, "hallmark_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
cat(sprintf("manifest -> %s (msigdbr %s)\n", manifest_path, MSIGDBR_VERSION))

provenance_path <- file.path(OUT_DIR, "PROVENANCE.md")
writeLines(c(
  "# Frozen MSigDB Hallmark gene sets",
  "",
  sprintf("HGNC symbols for %d MSigDB Hallmark (`H`) collection sets, Homo sapiens, pulled from the offline `msigdbr` %s package by exact `gs_name` and written one symbol per line, sorted.",
          length(SETS), MSIGDBR_VERSION),
  "",
  "| Set | n genes | Role in this compartment |",
  "|---|---|---|",
  "| `HALLMARK_HYPOXIA` | 200 | the gene-purge reference for the mouse `WT_heat` sets |",
  "| `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 | ER-side proteostasis readout |",
  "| `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 | inflammatory-signalling decomposition arm |",
  "| `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 | type-I interferon decomposition arm (cGAS/STING-adjacent) |",
  "| `HALLMARK_INFLAMMATORY_RESPONSE` | 200 | broad-inflammation decomposition arm |",
  "| `HALLMARK_IL2_STAT5_SIGNALING` | 199 | T-cell activation / growth-signalling decomposition arm |",
  "",
  "Each set carries a validated expected size in the generating script. A size mismatch is a hard stop, because the installed msigdbr data would no longer match the release these assets were frozen against and every downstream set-intersection tally would shift silently.",
  "",
  "These sets are used WHOLE. No taxonomy refinement is applied to any of them, which for a purge or decomposition test is the conservative direction: a larger curated set removes or claims more genes from the mouse signature, so it understates rather than inflates what survives.",
  "",
  "Reproduce from the compartment root:",
  "",
  "```bash",
  "Rscript 02_analysis/scripts/freeze_hallmark_sets.R",
  "```"
), provenance_path)
cat(sprintf("provenance -> %s\n", provenance_path))
