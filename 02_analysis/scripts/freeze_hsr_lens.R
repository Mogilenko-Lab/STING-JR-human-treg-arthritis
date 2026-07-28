#!/usr/bin/env Rscript
## freeze_hsr_lens.R — freeze the mouse-anchor human HSR lens to txt lists
## ===========================================================================

source("02_analysis/helpers/source_hash_manifest.R")
## Stage 10 reproducer. Materializes the curated human heat-shock-response lens
## from the mouse-anchor single source of truth:
##   ../mouse_anchor/00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human_lens.rds
##
## This script does NOT re-pull MSigDB. It freezes the exact byte-identical lens
## used by the mouse anchor into this compartment's one-symbol-per-line .txt idiom:
##   00_data/references/temp_hsr_lens/HSR_core.txt
##   00_data/references/temp_hsr_lens/HSR_sensitivity.txt
##   00_data/references/temp_hsr_lens/hsr_lens_manifest.csv
##   00_data/references/temp_hsr_lens/PROVENANCE.md
##
## Honest ceiling: even the clean HSR core is proteotoxic-stress-general, not
## fever-specific; only the mouse 37/39 contrast can measure thermal-ness.
##
## Run in-container from the compartment root or anywhere:
##   Rscript 02_analysis/scripts/freeze_hsr_lens.R
## ===========================================================================

## --- resolve compartment root relative to THIS script (self-contained per repo) ---
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg)) {
  normalizePath(file_arg)
} else {
  normalizePath("02_analysis/scripts/freeze_hsr_lens.R")
}
PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), "..", ".."))

SOURCE_RDS_REL <- "../mouse_anchor/00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human_lens.rds"
SOURCE_RDS <- normalizePath(file.path(PROJECT_ROOT, SOURCE_RDS_REL), mustWork = FALSE)
OUT_DIR <- file.path(PROJECT_ROOT, "00_data", "references", "temp_hsr_lens")
EXPECTED_SOURCE_SHA256 <- "0046bffbdd405860b12e1686a7f2d10bac7d4ba2640cd7574ed549bc894cf487"

if (!file.exists(SOURCE_RDS)) {
  stop(sprintf(
    "frozen mouse-anchor HSR lens is absent: %s\nExpected relative to compartment root: %s",
    SOURCE_RDS, SOURCE_RDS_REL
  ), call. = FALSE)
}
source_hash <- source_sha256(SOURCE_RDS)
if (!identical(source_hash, EXPECTED_SOURCE_SHA256)) {
  stop(sprintf(
    paste0("mouse-anchor HSR lens hash drift: expected %s, observed %s.\n",
           "Review the anchor source and update this freeze pin before regenerating."),
    EXPECTED_SOURCE_SHA256, source_hash
  ), call. = FALSE)
}

lens <- readRDS(SOURCE_RDS)
required <- c("HSR_core", "HSR_sensitivity")
missing <- setdiff(required, names(lens))
if (length(missing)) {
  stop(sprintf("HSR lens missing required set(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
}
if (length(lens$HSR_core) != 56) {
  stop(sprintf("HSR_core drift: expected 56 genes, found %d", length(lens$HSR_core)), call. = FALSE)
}
if (length(lens$HSR_sensitivity) != 176) {
  stop(sprintf("HSR_sensitivity drift: expected 176 genes, found %d", length(lens$HSR_sensitivity)), call. = FALSE)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

provenance <- attr(lens, "provenance")
msigdbr_version <- "unknown"
if (!is.null(provenance)) {
  candidates <- c("msigdbr_version", "msigdb_version", "msigdbr")
  hit <- candidates[candidates %in% names(provenance)]
  if (length(hit)) {
    msigdbr_version <- as.character(provenance[[hit[1]]])[1]
  }
}

manifest <- data.frame(
  set_name = character(),
  n_genes = integer(),
  source_rds = character(),
  source_sha256 = character(),
  msigdbr_version = character(),
  date_frozen = character(),
  stringsAsFactors = FALSE
)

for (set_name in required) {
  genes <- sort(unique(as.character(lens[[set_name]])))
  txt_path <- file.path(OUT_DIR, paste0(set_name, ".txt"))
  writeLines(genes, txt_path)
  cat(sprintf("%s: %d genes -> %s\n", set_name, length(genes), txt_path))
  manifest <- rbind(manifest, data.frame(
    set_name = set_name,
    n_genes = length(genes),
    source_rds = SOURCE_RDS_REL,
    source_sha256 = source_hash,
    msigdbr_version = msigdbr_version,
    date_frozen = "FROZEN",
    stringsAsFactors = FALSE
  ))
}

manifest_path <- file.path(OUT_DIR, "hsr_lens_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

provenance_path <- file.path(OUT_DIR, "PROVENANCE.md")
writeLines(c(
  "# Frozen human HSR lens",
  "",
  "These gene sets are the curated human heat-shock-response lens used as the JIA secondary annotation lens:",
  "",
  "- `HSR_core`: 56 HGNC symbols, the cleaned cytosolic/core HSR subset.",
  "- `HSR_sensitivity`: 176 HGNC symbols, the full sensitivity union.",
  "",
  sprintf("They were frozen from the mouse-anchor RDS `%s`.", SOURCE_RDS_REL),
  sprintf("Source SHA-256: `%s`.", source_hash),
  "This compartment does not re-pull MSigDB; the anchor RDS is the single source of truth so the JIA lens is byte-identical to the anchor.",
  "",
  "Honest ceiling: even the clean HSR core is proteotoxic-stress-general (HSF1 fires on oxidative, proteasome, and metal stress), not fever-specific. Only the mouse anchor's experimental 37/39 contrast can measure thermal-ness. In JIA we carry the lens and read it correlatively; we do not decompose temperature causality from human scRNA-seq.",
  "",
  "Reproduce from the compartment root:",
  "",
  "```bash",
  "Rscript 02_analysis/scripts/freeze_hsr_lens.R",
  "```"
), provenance_path)

cat(sprintf("manifest -> %s\n", manifest_path))
cat(sprintf("provenance -> %s\n", provenance_path))
