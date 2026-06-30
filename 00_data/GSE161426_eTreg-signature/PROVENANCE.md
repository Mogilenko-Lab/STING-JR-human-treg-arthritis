# Provenance — GSE161426  (human_treg_arthritis)

> Auto-generated from the STING-JR provenance crawl. Lives WITH the data in `00_data/GSE161426_eTreg-signature/`.

## Summary
- **Hosting:** GEO
- **Downloadable now:** True  |  **Auth required:** False
  - Public GEO series (Public on Feb 23 2021). Series-level suppl file is openly downloadable via FTP/HTTPS. No authentication needed. RAW reads (SRP292498 / PRJNA678290) are also public on SRA but are out of scope (reprocessing).
- **Estimated total size:** ~11 MB (single file: GSE161426_Gene_expression_table_log2.xlsx, exactly 11259561 bytes). No other GEO supplementary files exist (no _RAW.tar, no per-GSM files).
- **Modalities:** GEX=True TCR=False BCR=False ADT=False

## Recommended ingestion (no reprocessing from raw)
- **Format:** `expression_matrix`  |  **Reprocessing needed:** False
- **Rationale:** This dataset is a GENE-SIGNATURE source (Lutter/Mijnheer 2021, PMID 33976194), not a single-cell ingestion target. The ONLY processed artifact in GEO suppl is GSE161426_Gene_expression_table_log2.xlsx (11.3 MB): a log2-transformed normalized bulk RNA-seq table. VERIFIED structure: 32584 gene rows x 34 columns = 8 leading annotation columns + 26 sorted-population sample columns (sample colnames match the GSM titles exactly: adultPBTeff1-4, adultPBTreg1-4, aJIAPBTreg*, childPBTreg*, JIASFTeff1-4, JIASFTreg1-4, rJIAPBTreg*). The 26 expression columns range from -3.32 to 12.27 with ~48% negative, non-integer values, confirming log2 NORMALIZED (NOT raw counts). It is directly usable for deriving the human effector-Treg (eTreg) core signature by between-group comparison/ranking (SF Treg vs PB Treg, etc.), which is the project's intended use. There is NO raw count matrix and NO precomputed DE table in the GEO supplement, and all 26 GSMs have supplementary_file=NONE. Raw integer counts exist only as FASTQ in SRA (SRP292498/PRJNA678290) and would require reprocessing (explicitly NOT recommended). So reprocessing_needed=false for the intended signature-derivation use of the log2 table; raw-count AnnData/Seurat ingestion is not achievable from GEO suppl. Loader caveat: the original pandas.read_excel one-liner does NOT run in scdock-r-dev:v0.5.5 because openpyxl is absent in the base venv -- use R readxl (verified working) or install openpyxl.
- **Loader:** `# Sheet has 8 annotation cols (gene, RS, NA, chr, start, end, strand, transcript_length) + 26 sample cols. NOTE: scdock-r-dev:v0.5.5 base python lacks openpyxl, so pd.read_excel FAILS there. Use R readxl (verified) OR install openpyxl first. R: library(readxl); df <- read_excel('GSE161426_Gene_expression_table_log2.xlsx'); samp <- grep('Teff|Treg', colnames(df), value=TRUE)  # 26 sample columns; values are log2 NORMALIZED (NOT raw counts). Python alt: pip install openpyxl; df = pd.read_excel('GSE161426_Gene_expression_table_log2.xlsx'); set index to 'gene', subset the 26 Teff/Treg columns.`

## Test verification (throwaway scdock-r-dev:v0.5.5 container)
- attempted=True loaded_ok=True tool=scdock-r-dev:v0.5.5 R readxl::read_excel (docker run, timeout 240)
- file=`/data2/users/JCRLab/STING-JR/_provenance/GSE161426/GSE161426_Gene_expression_table_log2.xlsx` format=expression_matrix shape=32584 genes x 34 cols (8 annotation cols + 26 sample cols) integer_counts=False
- notes: Download OK: exactly 11259561 bytes, matches manifest. CORRECTION TO LOADER: the manifest's pandas.read_excel one-liner FAILS in scdock-r-dev:v0.5.5 -- ModuleNotFoundError: No module named 'openpyxl' (not installed in base venv, nor any /opt/venvs/* venv). Successfully loaded instead with R readxl (requireNamespace('readxl')==TRUE in this image). VERIFIED CONTENT: 32584 gene rows; 34 columns = 8 annotation columns (gene, RS, NA, chr, start, end, strand, transcript_length) + exactly 26 sample expression columns whose names match the GSM titles (adultPBTeff1-4, adultPBTreg1-4, aJIAPBTreg1/2/p, childPBTreg1/2/p, JIASFTeff1-4, JIASFTreg1-4, rJIAPBTreg1-3/p). The 26 sample expression columns range -3.32 to 12.27, ~48% negative, non-integer => confirmed log2 NORMALIZED values, NOT raw integer counts (integer_counts=false). recommended_ingestion.format=expression_matrix and reprocessing_needed=false are CORRECT; only the loader engine was corrected (use R readxl or pip install openpyxl). File is ~11 MB (<200 MB) so retained in _provenance, not deleted.

## Artifacts
| kind | format | direct? | n | size | where |
|---|---|---|---|---|---|
| processed_object | xlsx (log2 normalized bulk expression matrix; 32584 genes x [8 annotation cols + 26 sample cols]) | True | 1 | 11259561 bytes (~11 MB) | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161426/suppl/GSE161426_Gene_expression_table_log2.xlsx |
| signature | n/a (no precomputed signature/DE table present in GEO suppl) | False | 0 | 0 | NOT AVAILABLE in GEO supplement. The eTreg core signature must be DERIVED from GSE161426_Gene_expression_table_log2.x… |
| counts | n/a (no raw count matrix in GEO; FASTQ only) | False | 0 | see SRA | RAW integer counts NOT in GEO suppl. Only FASTQ via SRA SRP292498 / BioProject PRJNA678290 (https://www.ncbi.nlm.nih.… |

## Samples
26 bulk RNA-seq samples (sorted CD4 populations), Illumina NextSeq 500 (GPL18573), Homo sapiens. By tissue x population x group: Peripheral blood (PB): adult Teff = 4 (GSM4906863-66); adult Treg = 4 (GSM4906867-70); active-JIA (aJIA) child Treg = 3 (GSM4906871-73); healthy child Treg = 3 (GSM4906874-76); remission-JIA (rJIA) child Treg = 4 (GSM4906885-88). Synovial fluid (SF, JIA child): SF Teff = 4 (GSM4906877-80); SF Treg = 4 (GSM4906881-84). Totals by population: Treg = 18, Teff = 8. Totals by tissue: PB = 18, SF = 8. The key contrast for the eTreg signature is SF Treg vs PB Treg (inflamed-joint effector Treg programming). VERIFIED: all 26 sample column names present in the xlsx match these GSM titles.

## Open questions
- The xlsx contains log2 NORMALIZED expression only (verified: range -3.32 to 12.27, ~48% negative, non-integer) -- there is NO precomputed eTreg signature gene list and NO DE/statistics table in the GEO supplement. The signature must be derived from this matrix, or pulled from the paper's supplementary tables (Mijnheer/Lutter et al., Nat Commun 2021, PMID 33976194, DOI 10.1038/s41467-021-22975-7). Confirm whether the project wants to (a) recompute the signature from the log2 matrix or (b) use the published gene list from the paper supplement.
- LOADER/ENVIRONMENT: scdock-r-dev:v0.5.5 base python venv does NOT include openpyxl, so the original pandas.read_excel loader fails in-container. Verified the file loads via R readxl. For a python pipeline either install openpyxl in the target env or read via R; flagging so the ingestion recipe uses the correct engine.
- STRUCTURE: the table is 32584 genes x 34 columns = 8 annotation columns (gene, RS, NA, chr, start, end, strand, transcript_length) followed by the 26 sample columns -- the original manifest description ('genes x 26 samples') omitted the leading annotation block. When building the expression matrix, subset to the 26 Teff/Treg-named columns and use the 'gene' column as the index.
- Normalization/units of the log2 table are not documented in GEO metadata (likely log2 of normalized counts/CPM/TPM). Confirm exact normalization if precise units matter for signature scoring; presence of negative log2 values indicates normalized (not raw) data.
- Patient/donor IDs in the sample table are inferred from the GEO 'patient:' characteristic (integers 1-18). For SF, Teff and Treg share patient numbers 11-14 (paired sort from same donor); PB groups use distinct patient numbers. Verify pairing assumptions against the paper if donor-level pairing is needed.
- Raw integer counts are only available as FASTQ on SRA (SRP292498 / PRJNA678290); reprocessing is explicitly out of scope, so no raw-count AnnData/Seurat object can be built for this dataset.

## Download
See `download.sh` (recipe captured verbatim from the crawl). Data bytes stage to `/data2/users/JCRLab/STING-JR/human_treg_arthritis/GSE161426_eTreg-signature/` and mount to `00_data/GSE161426_eTreg-signature/raw/`.
