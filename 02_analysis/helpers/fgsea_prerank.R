#!/usr/bin/env Rscript
# fgsea_prerank.R — pre-ranked GSEA of gene sets against ONE ranked list.
# =======================================================================
# The primary-evidence step: score the mouse WT_heat up/down sets against a
# donor-pseudobulk SF-vs-PB ranked list (signed moderated t). Called as a
# subprocess by 05_score_signatures.py so the compute stays reproducible.
#
# COMPUTE ONLY — writes tables + a `gseaResult` RDS; NEVER plots. The static
# running-sum figure (05_score_signatures_viz.R) drives the canonical toolkit
# plotter off the RDS; the report widget reads the interactive runsum tables.
#
# Engine: clusterProfiler::GSEA(by = "fgsea"). This wraps the SAME fgsea engine
# the bare fgsea() call used, but returns a real DOSE/clusterProfiler
# `gseaResult` S4 object (with @result, @geneList, @geneSets, @params) that the
# RNAseq-toolkit running-sum plotter requires. NES is therefore near-identical
# to the previous bare-fgsea result (migration 2026-07-11).
#
# Usage:
#   Rscript fgsea_prerank.R <ranked.rnk> <out.csv> <contrast_label> \
#       <min_size> <max_size> <seed> <nperm> name1=genes1.txt [name2=genes2.txt ...]
#
# ranked.rnk : 2-col TSV (symbol \t stat), no header.
# Outputs (all under dirname(out.csv), stem = basename(out.csv) w/o .csv):
#   <out.csv>                              master_gsea_table schema NES table
#   <stem>.rds                             the gseaResult S4 object (for the plotter)
#   runsum_interactive_<tag>_<setname>.csv per-set interactive-widget substrate
#     (tag = pop tag parsed from `gsea_pseudobulk_<tag>.csv`; setname = e.g. WT_heat_up)
#
# master NES CSV columns:
#   pathway_id,pathway_name,database,nes,pvalue,padj,set_size,core_enrichment,contrast,direction

suppressPackageStartupMessages({
  library(clusterProfiler)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Toolkit helper: named-list -> TERM2GENE data frame (gs_name, gene_symbol).
TOOLKIT <- "01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_processing/pathway_utils.R"
suppressMessages(source(TOOLKIT))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8)
  stop("need: <rnk> <out.csv> <contrast> <min_size> <max_size> <seed> <nperm> name=genes.txt ...")

rnk_path  <- args[1]
out_csv   <- args[2]
contrast  <- args[3]
min_size  <- as.integer(args[4])
max_size  <- as.integer(args[5])
seed      <- as.integer(args[6])
nperm     <- as.integer(args[7])
set_specs <- args[8:length(args)]

database_for_set <- function(nm) {
  base <- sub("_(up|down)$", "", nm)
  if (base %in% c("WT_heat", "WT_heat_no_hypoxia", "unassigned")) {
    return("mouse_projection")
  }
  if (base %in% c("HSR_core", "HSR_sensitivity", "hsr_curated")) {
    return("curated_hsr_reactome_go")
  }
  if (base %in% c("upr_er", "hypoxia", "nfkb_tnfa", "ifn_type_i",
                  "inflammatory", "t_activation") || grepl("^HALLMARK_", base)) {
    return("msigdb_hallmark")
  }
  "custom_gene_set"
}

# --- ranked vector (named, unique, descending) ------------------------------
rnk <- read.table(rnk_path, sep = "\t", header = FALSE,
                  col.names = c("symbol", "stat"), stringsAsFactors = FALSE)
rnk <- rnk[!is.na(rnk$stat) & nzchar(rnk$symbol), ]
rnk <- rnk[!duplicated(rnk$symbol), ]
stats <- rnk$stat
names(stats) <- rnk$symbol
stats <- sort(stats, decreasing = TRUE)

# --- pathways: name=path/to/genes.txt ---------------------------------------
pathways   <- list()
directions <- c()
databases  <- c()
for (spec in set_specs) {
  kv <- strsplit(spec, "=", fixed = TRUE)[[1]]
  nm <- kv[1]; fp <- kv[2]
  genes <- trimws(readLines(fp, warn = FALSE))
  genes <- genes[nzchar(genes)]
  pathways[[nm]]  <- unique(genes)
  directions[nm]  <- if (grepl("_up$", nm)) "up" else if (grepl("_down$", nm)) "down" else "na"
  databases[nm]   <- database_for_set(nm)
}

# --- clusterProfiler::GSEA (by = "fgsea") -> a real gseaResult S4 object -----
# by="fgsea" uses the same fgsea engine; eps=0 gives exact (multilevel) p-values;
# pvalueCutoff=1 keeps ALL sets in @result (we report every set, not just FDR-sig).
t2g <- list_to_term2gene(pathways)  # columns: gs_name, gene_symbol
set.seed(seed)
gsea <- clusterProfiler::GSEA(
  geneList      = stats,
  TERM2GENE     = t2g,
  by            = "fgsea",
  exponent      = 1,
  eps           = 0,
  minGSSize     = min_size,
  maxGSSize     = max_size,
  nPermSimple   = nperm,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  seed          = TRUE,
  verbose       = FALSE
)

res <- gsea@result

# --- master-schema NES CSV from gsea@result ---------------------------------
# Include EVERY requested set (fill NA for any dropped by size filtering) so the
# downstream effect-size table always has both up + down rows.
row_for <- function(nm) {
  hit <- res[res$ID == nm, , drop = FALSE]
  if (nrow(hit) == 1) {
    data.frame(pathway_id = nm, pathway_name = nm, database = databases[[nm]],
               nes = hit$NES, pvalue = hit$pvalue, padj = hit$p.adjust,
               set_size = hit$setSize, core_enrichment = hit$core_enrichment,
               contrast = contrast, direction = directions[[nm]],
               stringsAsFactors = FALSE)
  } else {
    data.frame(pathway_id = nm, pathway_name = nm, database = databases[[nm]],
               nes = NA_real_, pvalue = NA_real_, padj = NA_real_,
               set_size = length(intersect(pathways[[nm]], names(stats))),
               core_enrichment = "", contrast = contrast,
               direction = directions[[nm]], stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, lapply(names(pathways), row_for))

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("[fgsea_prerank] %s: %d set(s) x %d ranked genes (nperm=%d) -> %s\n",
            contrast, length(pathways), length(stats), nperm, out_csv))

# --- persist the gseaResult S4 object for the toolkit running-sum plotter ----
rds_path <- sub("\\.csv$", ".rds", out_csv)
saveRDS(gsea, rds_path)
cat(sprintf("[fgsea_prerank] gseaResult -> %s\n", rds_path))

# --- interactive-widget substrate: per-set ranked running-sum tables ---------
# Recomputes the DOSE gseaScores running ES (weighted KS, exponent from the
# fitted object) so the emitted curve is IDENTICAL to the plotted one, and flags
# leading-edge/core genes via the object's own core_enrichment. One tidy CSV per
# gene set: the exact-genes substrate collaborators hover in the report widget.
stem     <- sub("\\.csv$", "", basename(out_csv))                 # gsea_pseudobulk_treg
pop_tag  <- sub("^gsea_pseudobulk_", "", stem)                    # treg / tcon / cd8
exponent <- as.numeric(gsea@params[["exponent"]] %||% 1)
gl       <- gsea@geneList                                         # named, sorted decreasing
N        <- length(gl)
gnames   <- names(gl)

for (nm in names(pathways)) {
  gs   <- intersect(pathways[[nm]], gnames)
  hits <- gnames %in% gs
  Nh   <- sum(hits)
  Phit <- numeric(N); Pmiss <- numeric(N)
  NR   <- sum(abs(gl[hits])^exponent)
  if (NR > 0) Phit[hits] <- (abs(gl[hits])^exponent) / NR
  Pmiss[!hits] <- 1 / max(N - Nh, 1L)
  running_es <- cumsum(Phit) - cumsum(Pmiss)

  core_str   <- res[res$ID == nm, "core_enrichment"]
  core_genes <- if (length(core_str) == 1 && !is.na(core_str) && nzchar(core_str))
                  strsplit(core_str, "/")[[1]] else character(0)

  runsum_df <- data.frame(
    rank         = seq_len(N),                       # 1-based position in the ranked list
    gene         = gnames,                            # HGNC symbol at this rank
    stat         = as.numeric(gl),                    # signed moderated t (ranking metric)
    running_es   = running_es,                        # DOSE weighted running enrichment score
    hit          = hits,                              # TRUE = gene is a member of this set
    leading_edge = hits & (gnames %in% core_genes),   # TRUE = core / leading-edge gene
    gene_set     = nm,
    population   = pop_tag,
    contrast     = contrast,
    stringsAsFactors = FALSE
  )
  ia_path <- file.path(dirname(out_csv),
                       sprintf("runsum_interactive_%s_%s.csv", pop_tag, nm))
  write.csv(runsum_df, ia_path, row.names = FALSE)
  cat(sprintf("[fgsea_prerank] interactive runsum -> %s\n", ia_path))
}
