#!/usr/bin/env bash
# Download recipe for GSE161426 (verbatim from STING-JR provenance crawl).
# Stage to /data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE161426_eTreg-signature/ then mount to 00_data/GSE161426_eTreg-signature/raw/.
set -euo pipefail

mkdir -p /data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE161426_eTreg-signature
curl -L --fail --max-time 600 -o /data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE161426_eTreg-signature/GSE161426_Gene_expression_table_log2.xlsx "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161426/suppl/GSE161426_Gene_expression_table_log2.xlsx"
# Verify size == 11259561 bytes:
# ls -l /data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE161426_eTreg-signature/GSE161426_Gene_expression_table_log2.xlsx
# LOAD (VERIFIED in scdock-r-dev:v0.5.5 via R readxl -- pandas needs openpyxl which is absent in the base venv):
#   Rscript -e 'library(readxl); df<-read_excel("GSE161426_Gene_expression_table_log2.xlsx"); samp<-grep("Teff|Treg",colnames(df),value=TRUE); print(length(samp))'  # 26 sample columns
# Python alt: pip install openpyxl; then pd.read_excel(..., index_col="gene")
# NOTE: This is the ONLY artifact. It contains log2 NORMALIZED expression (NOT raw counts) for 26 samples. Use it to DERIVE the eTreg signature (e.g., SF Treg vs PB Treg). No DE table or raw count matrix is hosted on GEO; do NOT attempt SRA/FASTQ reprocessing in this project.
