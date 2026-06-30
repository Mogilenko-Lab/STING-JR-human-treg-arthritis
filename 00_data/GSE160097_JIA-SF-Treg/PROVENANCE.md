# Provenance — GSE160097  (human_treg_arthritis)

> Auto-generated from the STING-JR provenance crawl. Lives WITH the data in `00_data/GSE160097_JIA-SF-Treg/`.

## Summary
- **Hosting:** GEO
- **Downloadable now:** True  |  **Auth required:** False
  - Public GEO series (Public on Nov 27 2020). No authentication needed; standard anonymous HTTPS/FTP from ncbi.nlm.nih.gov. No dbGaP-controlled access required for the processed per-GSM matrices and TCR contigs (those are the count-level data; raw FASTQs live in SRA but are NOT needed).
- **Estimated total size:** ~360 MB (GSE160097_RAW.tar = 377,845,760 bytes). Per-component if pulled individually: ~330 MB of filtered 10x H5 (40 files) + ~8.5 MB contig_annotations CSVs + ~23 MB contig FASTAs.
- **Modalities:** GEX=True TCR=True BCR=False ADT=False

## Recommended ingestion (no reprocessing from raw)
- **Format:** `10x_filtered_h5`  |  **Reprocessing needed:** False
- **Rationale:** No author-provided integrated/processed Seurat or h5ad object exists in GEO suppl (only the _RAW.tar of per-GSM 10x outputs). Each of the 40 GSMs ships a CellRanger filtered 10x H5 carrying RAW integer UMI counts: 33 GSMs use *_filtered_gene_bc_matrices_h5.h5 (older CellRanger v2 layout) and 7 use *_filtered_feature_bc_matrix.h5 (v3 layout). Both load directly via sc.read_10x_h5(path) -> AnnData with raw integer counts (VERIFIED in scdock-r-dev:v0.5.5 on GSM4859852 -> 3426x32738, float32 storage but all-integer values, max 1339.0). TCR is paired per-GSM via *_filtered_contig_annotations.csv.gz (10x VDJ annotations) for scirpy/Dandelion ingestion. No normalization-only artifacts. reprocessing_needed=false: ingest per-GSM H5 + concatenate, attaching donor/tissue/population from the sample table.
- **Loader:** `import scanpy as sc; ad = sc.read_10x_h5('GSM4859852_PM2_Treg_SF_p5_filtered_gene_bc_matrices_h5.h5')  # or *_filtered_feature_bc_matrix.h5 for the v3 GSMs; TCR: pandas.read_csv('GSM4859852_..._filtered_contig_annotations.csv.gz')`

## Test verification (throwaway scdock-r-dev:v0.5.5 container)
- attempted=True loaded_ok=True tool=scanpy.read_10x_h5 in scdock-r-dev:v0.5.5
- file=`GSM4859852_PM2_Treg_SF_p5_filtered_gene_bc_matrices_h5.h5` format=10x_filtered_h5 shape=3426x32738 integer_counts=True
- notes: VERIFIED. Smallest representative 10x_filtered_h5 (CellRanger v2 layout, 507325 bytes) downloaded to /data2/users/JCRLab/STING-JR/_provenance/GSE160097/ and loaded with sc.read_10x_h5 in scdock-r-dev:v0.5.5. Result: AnnData (3426 cells x 32738 genes), X dtype float32 but all values integer (np.all(X.data==floor)=True), max 1339.0 -> raw UMI counts confirmed, matching the manifest's claimed 3426x32738. Recommended format 10x_filtered_h5 is CORRECT; reprocessing_needed=false stands. Minor: anndata emits a 'Variable names are not unique' UserWarning -> apply .var_names_make_unique() on ingest. PATH CORRECTION discovered during verification: the per-GSM FTP path uses GSM4859nnn (drop last 3 digits), not GSM485nnn -- the original recipe path 404'd; recipe and artifact url_or_recipe fields updated. File kept (under 200MB, 507 KB).

## Artifacts
| kind | format | direct? | n | size | where |
|---|---|---|---|---|---|
| counts | 10x_filtered_h5 | True | 40 | ~10 MB each; ~330 MB total across 40 GSMs | Per-GSM filtered 10x H5 (raw integer UMI counts). 33 GSMs: *_filtered_gene_bc_matrices_h5.h5 (CellRanger v2). 7 GSMs … |
| tcr | 10x_vdj_contig_annotations_csv | True | 40 | ~115 KB-330 KB each; ~8.5 MB total (one tiny outlier GSM4859852=2.6KB) | Per-GSM *_filtered_contig_annotations.csv.gz (10x Genomics VDJ filtered contig annotations, paired TCR alpha/beta). B… |
| tcr | fasta_gz | False | 40 | ~6.5 KB-1.1 MB each; ~23 MB total | Per-GSM *_filtered_contig.fasta.gz (raw filtered contig nucleotide sequences). Companion to the contig annotations CS… |
| metadata | geo_soft_text | True | 1 | ~220 KB | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160097&targ=gsm&form=text&view=brief -- source of truth for gsm… |
| counts | raw_tar_bundle | False | 1 | 377845760 bytes (~360 MB) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160097/suppl/GSE160097_RAW.tar -- single tarball bundling all 12… |

## Samples
40 GSM = 6 sorted T-cell populations x 7 JIA patients (paired SF + blood), minus 2 missing samples (patient 3 absent for PB CD4 Tcon and PB CD8). By tissue x population: Synovial Fluid / CD4+ Tcon (CD4+CD45RO+CD25-): 7 (p1-7); Synovial Fluid / CD4+ Treg (CD4+CD45RO+CD127loCD25+): 7 (p1-7); Synovial Fluid / CD8+CD45RO+: 7 (p1-7); Peripheral Blood / CD4+ Tcon: 6 (p1,2,4,5,6,7 -- no p3); Peripheral Blood / CD4+ Treg: 7 (p1-7); Peripheral Blood / CD8+CD45RO+: 6 (p1,2,4,5,6,7 -- no p3). All donors are JIA (juvenile idiopathic arthritis) patients; no healthy controls in this series. Condition = JIA throughout. Treg = 14 GSMs (7 SF + 7 blood); Tcon = 13; CD8 = 13.

## Open questions
- No author-provided integrated/processed object (Seurat .rds or .h5ad with cluster/cell-type labels) is present in GEO suppl -- only per-GSM CellRanger outputs. If the publication's annotated object is needed (e.g. Treg subset labels, UMAP), it would have to be requested from authors or reconstructed by reprocessing from the per-GSM H5 + TCR.
- Two samples are intentionally absent (patient 3 has no PB CD4 Tcon and no PB CD8 -- GSM numbering skips p3 in the GSM4859844 and GSM4859871 groups). Confirm this is a biological/sampling dropout, not a missing upload.
- Mixed CellRanger output versions: 33 GSMs use *_filtered_gene_bc_matrices_h5.h5 (v2) and 7 use *_filtered_feature_bc_matrix.h5 (v3). Both load via sc.read_10x_h5 but gene/feature ID conventions and matrix dims differ across versions -- the verified v2 file (GSM4859852) has 32738 genes; verify var_names harmonize on concatenation and run .var_names_make_unique (non-unique var_names warning observed on load).
- GSM4859852 (SF Treg p5) H5 is tiny (507325 bytes, 3426 cells) and its contig files are minimal (FASTA 6538 bytes / CSV 2582 bytes) -- low-input sample; flag for QC but counts load fine (VERIFIED).
- CORRECTION (PATH): the per-GSM FTP directory token is GSM4859nnn (last 3 digits dropped from GSM4859852), NOT GSM485nnn as written in the original recipe. The original GSM485nnn path returns HTTP 404. Recipe/artifact URLs updated accordingly.
- filelist.txt was regenerated 2024-02-01 (RAW.tar timestamp) while per-file timestamps are 2020 (the verified H5 Last-Modified is 2020-10-07) -- no content change implied, but confirm md5s in filelist.txt match downloaded files.

## Download
See `download.sh` (recipe captured verbatim from the crawl). Data bytes stage to `/data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE160097_JIA-SF-Treg/` and mount to `00_data/GSE160097_JIA-SF-Treg/raw/`.
