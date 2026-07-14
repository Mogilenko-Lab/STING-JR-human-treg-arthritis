#!/usr/bin/env Rscript
# percell_score.R — per-cell signature scoring with AUCell + UCell.
# =======================================================================
# The rigorous per-cell SECONDARY lens: rank-based, composition-robust scores
# that replace the old scanpy score_genes module score. Donor-pseudobulk fgsea
# NES stays the PRIMARY evidence (see fgsea_prerank.R) — this is corroborative.
#
# Called as a subprocess (table/file interop, like fgsea_prerank.R):
# python exports the expression matrix + gene sets to disk, invokes this, then
# reads the per-cell CSV back.
#
# Usage:
#   Rscript percell_score.R <expr.mtx> <genes.txt> <barcodes.txt> <out.csv> \
#       <n_cores> name1=set1.txt [name2=set2.txt ...]
#
#   expr.mtx     : MatrixMarket sparse matrix, GENES x CELLS (features in rows,
#                  cells in columns) — the AUCell/UCell convention. Values are
#                  log-normalized expression (log1p CP10k); rank-based scoring is
#                  invariant to the monotone transform, so raw counts also work.
#   genes.txt    : HGNC symbols, one per line, in matrix ROW order (rownames).
#   barcodes.txt : cell ids, one per line, in matrix COLUMN order (colnames).
#   out.csv      : per-cell scores; columns = cell, <set>_AUCell, <set>_UCell ...
#   n_cores      : parallel workers (AUCell block BPPARAM + UCell ncores).
#   nameK=fileK  : gene set K named `nameK`, symbols one per line in fileK.
#
# Output CSV schema (one row per cell, aligned to barcodes.txt order):
#   cell, <set1>_AUCell, <set1>_UCell, <set2>_AUCell, <set2>_UCell, ...

suppressPackageStartupMessages({
  library(Matrix)
  library(AUCell)
  library(UCell)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("need: <expr.mtx> <genes.txt> <barcodes.txt> <out.csv> <n_cores> name=genes.txt ...")

mtx_path  <- args[1]
genes_path <- args[2]
bc_path   <- args[3]
out_csv   <- args[4]
n_cores   <- max(1L, as.integer(args[5]))
set_specs <- args[6:length(args)]

# --- expression matrix: genes (rows) x cells (cols) ---
mat <- readMM(mtx_path)
mat <- as(mat, "CsparseMatrix")            # dgCMatrix; AUCell/UCell both accept sparse
genes <- readLines(genes_path, warn = FALSE)
barcodes <- readLines(bc_path, warn = FALSE)
if (nrow(mat) != length(genes))
  stop(sprintf("row/gene mismatch: %d rows vs %d genes", nrow(mat), length(genes)))
if (ncol(mat) != length(barcodes))
  stop(sprintf("col/barcode mismatch: %d cols vs %d barcodes", ncol(mat), length(barcodes)))
rownames(mat) <- genes
colnames(mat) <- barcodes

# --- gene sets: name=path/to/genes.txt ---
gene_sets <- list()
for (spec in set_specs) {
  kv <- strsplit(spec, "=", fixed = TRUE)[[1]]
  nm <- kv[1]; fp <- kv[2]
  g <- readLines(fp, warn = FALSE)
  g <- trimws(g); g <- g[nzchar(g)]
  gene_sets[[nm]] <- unique(g)
}
if (length(gene_sets) == 0L) stop("no gene sets supplied")

bpp <- if (n_cores > 1L) MulticoreParam(workers = n_cores) else SerialParam()

# --- AUCell: build rankings once, then AUC per set ---
# splitByBlocks keeps peak memory bounded for large (100k+ cell) matrices.
rankings <- AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE,
                                 splitByBlocks = TRUE, BPPARAM = bpp)
auc <- AUCell_calcAUC(gene_sets, rankings, verbose = FALSE)
auc_mat <- t(getAUC(auc))                  # cells x sets
auc_mat <- auc_mat[barcodes, , drop = FALSE]

# --- UCell: rank-based Mann-Whitney U scores; chunked + parallel ---
uc_mat <- ScoreSignatures_UCell(mat, features = gene_sets,
                                ncores = n_cores, chunk.size = 5000)
# UCell names columns "<set>_UCell"; strip the suffix to realign by set name.
colnames(uc_mat) <- sub("_UCell$", "", colnames(uc_mat))
uc_mat <- uc_mat[barcodes, , drop = FALSE]

# --- assemble per-cell CSV, preserving set order ---
out <- data.frame(cell = barcodes, stringsAsFactors = FALSE, check.names = FALSE)
for (nm in names(gene_sets)) {
  out[[paste0(nm, "_AUCell")]] <- if (nm %in% colnames(auc_mat)) auc_mat[, nm] else NA_real_
  out[[paste0(nm, "_UCell")]]  <- if (nm %in% colnames(uc_mat))  uc_mat[, nm]  else NA_real_
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("[percell_score] %d sets x %d cells (%d genes, %d cores) -> %s\n",
            length(gene_sets), length(barcodes), length(genes), n_cores, out_csv))
