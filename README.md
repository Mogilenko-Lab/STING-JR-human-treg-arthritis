# JIA sorted synovial and blood T cells — GSE160097

Single-cell RNA-seq of FACS-sorted Treg, Tcon and CD8 T cells from the synovial fluid and the paired
peripheral blood of the same juvenile idiopathic arthritis patients.

**The question.** Does a mouse 39 °C-derived signature separate the inflamed synovial niche from
paired blood within a frozen cell state, and does any part of that separation hold up against the
other stresses the same joint imposes?

This is the Treg-relevant compartment of a project testing whether human inflammatory and
autoinflammatory disease states carry programs consistent with a temperature-stress axis in T cells
and with cGAS–STING-related biology. The signatures scored here come from a mouse 2×2 temperature ×
cGAS-genotype experiment (GSE329522), and every gene set in the tree originates outside the JIA data.

**Abbreviations.** SF = synovial fluid (inflamed joint). PB = peripheral blood. Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory. Tcon = CD4⁺CD25⁻ conventional. CD8 = CD8⁺CD45RO⁺ memory. NES =
normalised enrichment score.

## Design

[GSE160097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160097) — forty libraries, seven
donors, six paired in each sorted population after QC. That donor pairing is what makes this a niche
contrast: synovial fluid against blood within the same patient and the same sorted population. A
second accession, [GSE161426](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161426), supplies
an effector Treg signature derived here from its deposited matrix.

The analysis builds and QCs the object, freezes the sort labels, aggregates to donor-level pseudobulk,
fits the paired contrast, and scores the mouse arms against those rankings. The rest of the tree
presses on the answer: purging hypoxia gene content, decomposing the arm against curated lenses,
running the unbiased sweep for calibration, and localising each score on two embeddings.

## Two tracks, ranked

**The confirmatory spine carries every claim.** Donor-level pseudobulk on raw counts within frozen
sort labels, limma-voom, then pre-ranked enrichment on the moderated *t*. Each donor casts one vote,
and effect sizes accumulate in `03_results/master/`.

**Per-cell readouts are annotation.** A per-cell score localises a program on a map, pooling thousands
of cells from donors of unequal yield. Read those figures for *where* a score sits, and the
confirmatory ones for whether it *separates*.

## Findings

- **The mouse 39 °C up arm separates the niche in every sorted population.** NES 2.59 in Treg on 120
  of 202 arm genes, 2.68 in Tcon on 131, 2.06 in CD8 on 114, every pooled FDR below 1e-4. The Treg
  score sits between the Tcon and CD8 scores, so the separation is pan-T, with Treg one of the three
  populations carrying it.
- **Both arms carry information.** The down arm reaches NES 1.432 at FDR 0.0354 on 64 genes in Tcon,
  at the same sign as the up arm — a shared non-directional shift.
- **By composition the arm is largely inflammatory.** Nine curated anchor-independent lenses claim 67
  of the 202 `WT_heat_up` genes: 35 TNFα/NF-κB, 21 inflammatory-response, 2 in the curated heat-shock
  core, and 2 of the 21 published interferon-independent STING genes. The remaining 135 stay
  unclaimed, the largest single part.
- **The enrichment survives deleting its hypoxia gene content.** Removing the 18 HALLMARK_HYPOXIA
  overlap genes costs 0.126 to 0.164 NES and leaves all three populations significant, which bounds
  the result to gene content.
- **The niche moves an enormous amount at once.** 1,443 of 11,514 tests reach pooled FDR 0.05 in Treg,
  2,043 of 11,752 in Tcon and 939 of 11,532 in CD8, and the largest single effects run toward blood.
  Among the programs moving toward the joint in sorted CD4 cells the arm sits near the top of the
  distribution, and that calibration is what makes the enrichment readable.
- **Temperature and hypoxia are jointly imposed by the inflamed joint** and stay entangled in
  cross-sectional human data.

Each number is reproduced in a table under `03_results/`, named in that stage's `README.md`.

## Two implementation details

**Pseudobulk differential expression runs in R, across a file seam.** Python aggregates raw integer UMI
counts to the donor × condition × label matrix and writes plain CSVs; R reads those, fits
edgeR/limma-voom, and writes the DE tables and signed ranked lists. The seam carries a gene map on
purpose: counts are keyed by Ensembl id while reference sets match on HGNC symbol, so the R side joins
the map before ranking. Skipping that join yields ranked lists that intersect the references at
approximately zero, which enrichment tools report as empty results.

**A gene symbol that fails to match is a vocabulary result until proven a biological one.** The count
matrix is frozen to the symbol vocabulary of the reference it was quantified against, while reference
sets ship current symbols — here cGAS is `MB21D1` and STING is `TMEM173`. Aliases resolve where
reference symbols enter, and matching is reported as a three-way ledger: matched, matched-via-alias,
genuinely absent.

## Layout

| Path | Contents |
|---|---|
| `00_data/` | Per accession: a provenance record, a sample table, a manifest and the download script. Raw bytes stage off-repo and mount read-only, so `raw/` carries a placeholder here. |
| `02_analysis/scripts/` | 54 numbered scripts, 34 in Python and 20 in R. Compute scripts write tables; `_viz` scripts draw from them. |
| `02_analysis/config/` | `analysis_config.yaml` — every path and threshold the scripts read. |
| `02_analysis/helpers/` | Pseudobulk aggregation, gene-set handling, figure styling, symbol-alias resolution. |
| `02_analysis/notebooks/` | Working notebooks for QC, annotation and score review. |
| `03_results/` | One directory per stage, each with `tables/`, `figures/` and a `README.md` captioning every file. |

`03_results/README.md` is the entry point: reading order stage by stage, each finding with the table it
comes from, and the accession, derivation and citation for every gene set scored.

## Reproducing

Scripts run in numeric order inside the container under `.devcontainer/`, reading paths and thresholds
from `02_analysis/config/analysis_config.yaml`. `00_data/<accession>/download.sh` fetches the data from
GEO, with the manifest and sample table recording what should arrive. Tables and stage READMEs are the
tracked record — figures regenerate from the script each caption names, since `.png` and `.pdf` are
ignored repo-wide. Clone with `--recursive` to populate `01_modules/`, which GitHub source tarballs
carry as empty directories.

## Environment

The analysis runs in `scdock-r-dev:v0.5.10`, pinned on the `dev-core` service in
`.devcontainer/docker-compose.yml`. That image is defined by
[scbio-docker](https://github.com/tony-zhelonkin/scbio-docker) at commit
[`5885cd3`](https://github.com/tony-zhelonkin/scbio-docker/commit/5885cd306ea908cb1949e7238b9186074b938953).

[RNAseq-toolkit](https://github.com/tony-zhelonkin/RNAseq-toolkit) supplies the canonical volcano and
GSEA running-sum plotters, plus the MitoCarta 3.0 human gene-set build. This release records it at
commit
[`752481f`](https://github.com/tony-zhelonkin/RNAseq-toolkit/commit/752481fd13542ccb81d2b9b92ba57305cf13d6fc)
(`v0.2.0-9-g752481f`, on `dev`) under `01_modules/`.

## Data availability

Both accessions are public: [GSE160097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160097)
(Mijnheer / Lutter et al. 2021, *Nature Communications*, doi 10.1038/s41467-021-22975-7) and
[GSE161426](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161426). Every external gene set,
lens and network is listed with its version, derivation and citation in `03_results/README.md`.

## License

**MIT** for code — `02_analysis/`, `.devcontainer/`, `config.py`, and scripts anywhere in the tree
([LICENSE](LICENSE)). **CC BY 4.0** for results and prose — the tables, figures and README text under
`03_results/`, and the documentation ([LICENSE-CC-BY-4.0.txt](LICENSE-CC-BY-4.0.txt)). Both require
attribution. Sequence data carries the terms of its GEO deposition.
