#!/usr/bin/env Rscript
# =============================================================================
# 03b_pseudobulk_de.R  --  limma-voom DE on pseudobulk CSVs
# =============================================================================
# Reads raw integer pseudobulk counts emitted by Python (03a) and runs standard
# edgeR/limma-voom DE.
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
})

YAML_CONFIG <- yaml::read_yaml("02_analysis/config/analysis_config.yaml")


STAGE <- "03_pseudobulk"
tdir <- file.path("03_results", STAGE, "tables")

POP_TAG <- c(
    "Treg" = "treg", 
    "Tcon" = "tcon", 
    "CD8" = "cd8"
)

COND_KEY <- YAML_CONFIG$design$tissue_key
COND_NUM <- YAML_CONFIG$design$tissue_levels$synovial_fluid
COND_DEN <- YAML_CONFIG$design$tissue_levels$peripheral_blood
DONOR_KEY <- YAML_CONFIG$design$donor_key

counts_path <- file.path(tdir, "pseudobulk_counts.csv")
coldata_path <- file.path(tdir, "pseudobulk_coldata.csv")
genemap_path <- file.path(tdir, "gene_map.csv")

stopifnot(file.exists(counts_path) && file.exists(coldata_path))
# The counts matrix is keyed by Ensembl id; every consumer of ranked_*.tsv matches on
# HGNC symbol. Without this map the ranked lists would carry ENSG ids and silently
# intersect the reference gene sets at ~zero.
stopifnot(file.exists(genemap_path))

counts_df <- read.csv(counts_path, row.names=1, check.names=FALSE)
coldata <- read.csv(coldata_path, row.names=1, stringsAsFactors=FALSE)
gene_map <- read.csv(genemap_path, row.names=1, stringsAsFactors=FALSE)

common_strata <- intersect(rownames(coldata), rownames(counts_df))
coldata <- coldata[common_strata, , drop=FALSE]
counts_mat <- t(as.matrix(counts_df[common_strata, ]))  

summary_list <- list()

for (pop in names(POP_TAG)) {
  tag <- POP_TAG[[pop]]
  cd <- coldata[coldata$coarse_label == pop & coldata[[COND_KEY]] %in% c(COND_NUM, COND_DEN), , drop=FALSE]
  
  if (nrow(cd) == 0) next
  
  n_num <- sum(cd[[COND_KEY]] == COND_NUM)
  n_den <- sum(cd[[COND_KEY]] == COND_DEN)
  
  if (min(n_num, n_den) < 2) {
    cat(sprintf("[03b_pseudobulk_de] %s: underpowered (SF=%d, PB=%d) — skipping DE\n", pop, n_num, n_den))
    summary_list[[pop]] <- data.frame(population=pop, n_sf=n_num, n_pb=n_den, model="skipped", n_sig_de=0, n_ranked=0)
    next
  }
  
  c_mat <- counts_mat[, rownames(cd), drop=FALSE]
  
  cd[[COND_KEY]] <- factor(cd[[COND_KEY]], levels = c(COND_DEN, COND_NUM))
  cd[[DONOR_KEY]] <- factor(cd[[DONOR_KEY]])
  
  donors_per_arm <- table(cd[[DONOR_KEY]], cd[[COND_KEY]])
  shared_donors <- sum(rowSums(donors_per_arm > 0) == 2)
  
  if (shared_donors >= 2) {
    design_fmla <- as.formula(sprintf("~ %s + %s", DONOR_KEY, COND_KEY))
    design <- model.matrix(design_fmla, data = cd)
    model_str <- sprintf("~ %s + %s", DONOR_KEY, COND_KEY)
  } else {
    design_fmla <- as.formula(sprintf("~ %s", COND_KEY))
    design <- model.matrix(design_fmla, data = cd)
    model_str <- sprintf("~ %s", COND_KEY)
  }
  
  dge <- DGEList(counts = c_mat)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust=TRUE)
  
  coef_name <- sprintf("%s%s", COND_KEY, COND_NUM)
  res <- topTable(fit, coef=coef_name, number=Inf, sort.by="none")
  res$ensembl_id <- rownames(res)
  res$gene_symbol <- gene_map[res$ensembl_id, "gene_symbol"]

  # Canonical, engine-agnostic DE schema — this is the seam contract. Consumers
  # (03_pseudobulk_de_viz.py, 09_heat_hypoxia_viz.py) read these names, so swapping the
  # engine again must not ripple into them. limma's native columns are kept alongside.
  # avg_expr is limma's AveExpr (log2-CPM); it is deliberately NOT called baseMean,
  # which in DESeq2 is a normalised mean count on a different scale.
  res$log2FoldChange <- res$logFC
  res$stat           <- res$t
  res$pvalue         <- res$P.Value
  res$padj           <- res$adj.P.Val
  res$avg_expr       <- res$AveExpr
  res$model          <- model_str
  res$n_paired_donors <- shared_donors
  res$de_engine      <- "limma-voom"

  res <- res[order(res$pvalue), ]
  res <- res[, c("ensembl_id", "gene_symbol", "avg_expr", "log2FoldChange", "stat",
                 "pvalue", "padj", "B", "model", "n_paired_donors", "de_engine")]
  write.csv(res, file.path(tdir, sprintf("de_SFvsPB_%s.csv", tag)), row.names=FALSE)
  
  # Ranked list for fgsea by t-statistic
  res_ranked <- res[!is.na(res$t) & !is.na(res$gene_symbol) & res$gene_symbol != "nan", ]
  res_ranked <- res_ranked[order(abs(res_ranked$t), decreasing=TRUE), ]
  res_ranked <- res_ranked[!duplicated(res_ranked$gene_symbol), ]
  res_ranked <- res_ranked[order(res_ranked$t, decreasing=TRUE), c("gene_symbol", "t")]
  write.table(res_ranked, file.path(tdir, sprintf("ranked_%s.tsv", tag)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  n_sig <- sum(res$padj < YAML_CONFIG$thresholds$de_fdr & abs(res$log2FoldChange) >= 1.0, na.rm=TRUE)
  cat(sprintf("[03b_pseudobulk_de] %s: model %s; %d sig DE (padj<%.2f, |lfc|>=1); ranked %d genes\n", 
              pop, model_str, n_sig, YAML_CONFIG$thresholds$de_fdr, nrow(res_ranked)))
  
  summary_list[[pop]] <- data.frame(population=pop, n_sf=n_num, n_pb=n_den, model=model_str, n_sig_de=n_sig, n_ranked=nrow(res_ranked))
}

if (length(summary_list) > 0) {
  sum_df <- do.call(rbind, summary_list)
  write.csv(sum_df, file.path(tdir, "de_summary.csv"), row.names=FALSE)
}
cat("[03b_pseudobulk_de] done\n")
