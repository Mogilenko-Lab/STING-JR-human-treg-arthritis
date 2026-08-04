#!/usr/bin/env Rscript
## freeze_tcr_activation_lens.R — freeze the curated human TCR activation lens to a txt list
## ===========================================================================
source("02_analysis/helpers/source_hash_manifest.R")
## Reproducer for the activation-pole reference list this compartment scores in the
## unbiased sweep. Materializes the curated human TCR/IEG T-cell activation lens from
## the mouse-anchor single source of truth:
##   ../mouse_anchor/00_data/references/gene_sets/tcr_activation_lens/tcr_activation_human.rds
##
## The HUMAN asset is the source panel, not a conversion of anything. The anchor's lens
## is defined by a frozen 66-gene HUMAN symbol panel (tcr_activation_panel.csv, spanning
## TCR-proximal signalling, early costimulation, immediate-early transcription factors
## and activation effectors, with FOXP3 dropped as a lineage-identity marker rather than
## an activation marker); the anchor's mouse list is the msigdbr-native ortholog
## conversion OF that panel. So reading the human asset here involves no ortholog step
## and is one derivation closer to the curated source than the mouse anchor's own use.
##
## This script re-pulls nothing. It freezes the byte-identical anchor asset into this
## compartment's one-symbol-per-line .txt idiom, the same way freeze_hsr_lens.R does:
##   00_data/references/tcr_activation_lens/TCR_activation.txt
##   00_data/references/tcr_activation_lens/tcr_activation_manifest.csv
##   00_data/references/tcr_activation_lens/PROVENANCE.md
##
## Honest ceiling: the lens marks T-cell activation, which the synovial niche imposes
## alongside every other stress it imposes. It is carried as an independent curated
## comparator and read correlatively.
##
## Run in-container from the compartment root or anywhere:
##   Rscript 02_analysis/scripts/freeze_tcr_activation_lens.R
## ===========================================================================

## --- resolve compartment root relative to THIS script (self-contained per repo) ---
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg)) {
  normalizePath(file_arg)
} else {
  normalizePath("02_analysis/scripts/freeze_tcr_activation_lens.R")
}
PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), "..", ".."))

SOURCE_RDS_REL <- "../mouse_anchor/00_data/references/gene_sets/tcr_activation_lens/tcr_activation_human.rds"
SOURCE_RDS <- normalizePath(file.path(PROJECT_ROOT, SOURCE_RDS_REL), mustWork = FALSE)
OUT_DIR <- file.path(PROJECT_ROOT, "00_data", "references", "tcr_activation_lens")
EXPECTED_SOURCE_SHA256 <- "fa36441a283a855ac8031e6073ef89fd4f472b9223c85e2b3c41f6a3bb286cf3"

## The panel size is asserted rather than read, because the whole point of a frozen
## reference is that a silent drift in the anchor changes what every downstream
## enrichment statistic was computed against.
EXPECTED_N_GENES <- 66L
SET_NAME <- "TCR_activation"

if (!file.exists(SOURCE_RDS)) {
  stop(sprintf(
    "frozen mouse-anchor TCR activation lens is absent: %s\nExpected relative to compartment root: %s",
    SOURCE_RDS, SOURCE_RDS_REL
  ), call. = FALSE)
}
source_hash <- source_sha256(SOURCE_RDS)
if (!identical(source_hash, EXPECTED_SOURCE_SHA256)) {
  stop(sprintf(
    paste0("mouse-anchor TCR activation lens hash drift: expected %s, observed %s.\n",
           "Review the anchor source and update this freeze pin before regenerating."),
    EXPECTED_SOURCE_SHA256, source_hash
  ), call. = FALSE)
}

lens <- readRDS(SOURCE_RDS)
if (!SET_NAME %in% names(lens)) {
  stop(sprintf("TCR activation lens missing required set: %s (found: %s)",
               SET_NAME, paste(names(lens), collapse = ", ")), call. = FALSE)
}
genes <- sort(unique(as.character(lens[[SET_NAME]])))
if (length(genes) != EXPECTED_N_GENES) {
  stop(sprintf("%s drift: expected %d genes, found %d", SET_NAME, EXPECTED_N_GENES,
               length(genes)), call. = FALSE)
}
## An Ensembl-keyed or mouse-cased list would intersect this compartment's HGNC-keyed
## ranked lists at approximately zero, and fgsea reports that as an empty result rather
## than as an error. Cheap to check here, expensive to discover downstream.
if (any(grepl("^ENS[A-Z]*G[0-9]{6,}$", genes))) {
  stop("TCR activation lens carries Ensembl ids; this compartment matches on HGNC symbol.",
       call. = FALSE)
}
if (!all(genes == toupper(genes))) {
  stop(sprintf(paste0("TCR activation lens carries non-uppercase symbols (%s), which is the ",
                      "mouse-cased asset rather than the human panel."),
               paste(utils::head(genes[genes != toupper(genes)], 5), collapse = ", ")),
       call. = FALSE)
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

txt_path <- file.path(OUT_DIR, paste0(SET_NAME, ".txt"))
writeLines(genes, txt_path)
cat(sprintf("%s: %d genes -> %s\n", SET_NAME, length(genes), txt_path))

manifest <- data.frame(
  set_name = SET_NAME,
  n_genes = length(genes),
  source_rds = SOURCE_RDS_REL,
  source_sha256 = source_hash,
  msigdbr_version = msigdbr_version,
  date_frozen = "FROZEN",
  stringsAsFactors = FALSE
)
manifest_path <- file.path(OUT_DIR, "tcr_activation_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

provenance_path <- file.path(OUT_DIR, "PROVENANCE.md")
writeLines(c(
  "# Frozen human TCR activation lens",
  "",
  "One gene set, used as an independent curated comparator in the unbiased enrichment sweep:",
  "",
  sprintf("- `%s`: %d HGNC symbols spanning TCR-proximal signalling, early costimulation, immediate-early transcription factors and activation effectors.", SET_NAME, length(genes)),
  "",
  sprintf("Frozen from the mouse-anchor RDS `%s`.", SOURCE_RDS_REL),
  sprintf("Source SHA-256: `%s`.", source_hash),
  "",
  "The HUMAN list is the curated source panel. The anchor defines this lens as a frozen human-symbol panel and derives its mouse list from that panel by msigdbr-native ortholog conversion, so the list frozen here involves no ortholog step. This compartment re-pulls nothing; the anchor RDS is the single source of truth, so the list is byte-identical to the anchor's human asset.",
  "",
  "`FOXP3` is deliberately absent: it marks Treg lineage identity rather than activation.",
  "",
  "Honest ceiling: the lens marks T-cell activation, which the inflamed synovial niche imposes alongside the other stresses it imposes. Overlap of any empirical arm with this lens is a statement about shared gene content. Read correlatively.",
  "",
  "Reproduce from the compartment root:",
  "",
  "```bash",
  "Rscript 02_analysis/scripts/freeze_tcr_activation_lens.R",
  "```"
), provenance_path)

cat(sprintf("manifest -> %s\n", manifest_path))
cat(sprintf("provenance -> %s\n", provenance_path))
