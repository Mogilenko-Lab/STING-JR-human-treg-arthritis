#!/usr/bin/env bash
# Download recipe for GSE160097 (verbatim from STING-JR provenance crawl).
# Stage to /data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE160097_JIA-SF-Treg/ then mount to 00_data/GSE160097_JIA-SF-Treg/raw/.
set -euo pipefail

## OPTION A (recommended, single shot): pull the _RAW.tar and untar into the staging dir
DEST=/data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE160097_JIA-SF-Treg
mkdir -p "$DEST"
curl -L --fail --max-time 1800 -o "$DEST/GSE160097_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160097/suppl/GSE160097_RAW.tar"
# verify size = 377845760 bytes, then:
tar -xvf "$DEST/GSE160097_RAW.tar" -C "$DEST"
# Result: 120 files = 40x *_filtered_*.h5 (raw counts) + 40x *_filtered_contig_annotations.csv.gz (TCR) + 40x *_filtered_contig.fasta.gz

## OPTION B (per-GSM FTP, if you want only H5 + TCR CSV and skip the FASTAs):
## NOTE (CORRECTED): the GSM directory token drops the LAST THREE digits of the accession, so GSM4859852 -> GSM4859nnn (NOT GSM485nnn). Earlier recipe used GSM485nnn which 404s.
DEST=/data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE160097_JIA-SF-Treg
mkdir -p "$DEST"
# Files are listed in filelist.txt with their exact <GSM>_<title>_ prefixes (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160097/suppl/filelist.txt).
for GSM in $(seq 4859835 4859874); do
  GSM="GSM$GSM"
  BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4859nnn/$GSM/suppl/"
  for f in $(curl -s --max-time 30 "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160097/suppl/filelist.txt" | grep -oE "$GSM[^[:space:]]*(_filtered_gene_bc_matrices_h5\.h5|_filtered_feature_bc_matrix\.h5|_filtered_contig_annotations\.csv\.gz)"); do
    curl -L --fail --max-time 600 -o "$DEST/$f" "$BASE$f"
  done
done
# (drop the contig_annotations pattern from the grep if you want H5 only)

## Authoritative reference for sizes/md5/exact filenames: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160097/suppl/filelist.txt
## VERIFIED single-file download used for test_verification:
## curl -s --max-time 600 -o /data2/users/JCRLab/STING-JR/_provenance/GSE160097/GSM4859852_PM2_Treg_SF_p5_filtered_gene_bc_matrices_h5.h5 "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4859nnn/GSM4859852/suppl/GSM4859852_PM2_Treg_SF_p5_filtered_gene_bc_matrices_h5.h5"
