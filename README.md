# JIA sorted synovial and blood T cells — GSE160097

Single-cell RNA-seq analysis of FACS-sorted Treg, Tcon and CD8 T cells from the synovial fluid and
the paired peripheral blood of the same juvenile idiopathic arthritis patients. This repository
asks whether a mouse temperature-derived signature separates the inflamed joint from blood within
a frozen cell state, and whether any part of that separation resists being reduced to the other
stresses the same joint imposes.

It is the most directly Treg-relevant compartment of a larger project testing whether human
inflammatory and autoinflammatory disease states contain programs consistent with a mouse 39 °C
T-cell stress axis and with cGAS–STING-related biology. The signatures scored here are derived
elsewhere, from a mouse 2×2 temperature × cGAS-genotype experiment (GSE329522), and every gene set
in the tree comes from outside the JIA data.

**Abbreviations.** SF = synovial fluid (inflamed joint). PB = peripheral blood. Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory. Tcon = CD4⁺CD25⁻ conventional. CD8 = CD8⁺CD45RO⁺ memory.
NES = normalised enrichment score.

## The data

[GSE160097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160097) — forty libraries,
seven donors, six paired in each sorted population after QC. That donor pairing is what lets the
compartment ask a niche question rather than a cohort question. A second accession,
[GSE161426](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161426), supplies an effector
Treg signature derived here from its deposited matrix.

## Two tracks, and how they rank

This is the load-bearing distinction in the repository, and it is worth reading before any figure.

**The confirmatory spine is the only track that supports a claim.** Donor-level pseudobulk on raw
counts within frozen sort labels, limma-voom, then pre-ranked enrichment on the moderated *t*. Each
donor casts one vote. Its effect sizes accumulate in `03_results/master/`.

**Everything per-cell is annotation.** A per-cell score localises a program on a map. It pools
thousands of cells from donors of unequal yield, so a tissue difference read off a colour or a
violin is pseudoreplicated. Read those figures for *where* a score sits; read the confirmatory ones
for whether it *separates*.

## What the analysis found

**The mouse 39 °C up arm separates the niche, and it does so in every sorted population.** NES 2.59
in Treg on 120 of 202 arm genes, 2.68 in Tcon on 131, 2.06 in CD8 on 114, every pooled FDR below
1e-4. The Treg score sits between the Tcon and CD8 scores, so the separation is pan-T rather than
Treg-selective — Treg is one of three populations carrying it.

**Its up and down arms move together.** The down arm reaches NES 1.432 at FDR 0.0354 on 64 genes in
Tcon, at the same sign as the up arm, and carries no direction in Treg or CD8. Both arms carry
information, and the pattern is a shared non-directional shift rather than a recapitulation of the
mouse contrast.

**By composition the arm is largely inflammatory.** Nine curated anchor-independent lenses contain
67 of the 202 `WT_heat_up` genes, leaving 135 unclaimed as the largest single part. What the lenses
do claim runs 35 TNFα/NF-κB genes and 21 inflammatory-response genes against 2 in the curated
heat-shock core and 2 of the 21 published interferon-independent STING genes.

**Its enrichment survives deleting its hypoxia gene content.** Removing the 18 HALLMARK_HYPOXIA
overlap genes costs 0.126 to 0.164 NES and leaves all three populations significant. That is a
statement about gene content and nothing more.

**The niche moves an enormous amount at once.** 1,443 of 11,514 tests reach pooled FDR 0.05 in Treg,
2,043 of 11,752 in Tcon and 939 of 11,532 in CD8, and the largest single effects run toward blood.
Among the programs moving toward the joint in sorted CD4 cells the arm sits near the top of the
distribution — that calibration is what makes the enrichment interpretable.

**Temperature and hypoxia are jointly imposed by the inflamed joint** and stay entangled in
cross-sectional human data. Nothing here separates them, and no artifact in this tree asserts that
one is a confound of the other.

Each number is reproduced in a table under `03_results/`, named in that stage's `README.md`.

## Layout

| Path | Contents |
|---|---|
| `00_data/` | Read-only inputs. Per accession: a provenance record, a sample table, a manifest and the download script. Raw bytes are staged off-repo and mounted read-only, so `raw/` carries only a placeholder here. |
| `02_analysis/scripts/` | 54 numbered scripts, 34 in Python and 20 in R. Every artifact under `03_results/` is reproducible from one of them; a `_viz` script never computes and a compute script never plots. |
| `02_analysis/config/` | `analysis_config.yaml` holds every path and threshold. Nothing is hardcoded in a script. |
| `02_analysis/helpers/` | Shared functions: pseudobulk aggregation, gene-set handling, figure styling, symbol-alias resolution. |
| `02_analysis/notebooks/` | Working notebooks for QC, annotation and score review. |
| `03_results/` | One directory per analysis stage, each with `tables/`, `figures/` and a `README.md` captioning every file. `03_results/README.md` is the reading order. |
| `.devcontainer/` | The container the analysis runs in. |

`03_results/README.md` is the best entry point: reading order stage by stage, the findings with the
table each comes from, and the accession, derivation and citation for every gene set scored.

## Two implementation details that matter

**Pseudobulk differential expression runs in R, across an explicit file seam.** A Python script
aggregates raw integer UMI counts to the donor × condition × label matrix and writes plain CSVs,
computing no statistics; an R script reads those, runs edgeR/limma-voom, and writes the DE tables
and signed ranked lists. The seam carries a gene map on purpose. Count matrices are keyed by
Ensembl id while every reference gene set matches on HGNC symbol, so ranking without joining the
map produces lists that intersect the references at approximately zero — and enrichment tools
report that as empty rather than as an error, so it fails silently and looks like a biological
null.

**A gene symbol that fails to match is a vocabulary result until proven a biological one.** The
count matrix is frozen to the symbol vocabulary of the reference it was quantified against, and
reference gene sets ship current symbols. Here cGAS is `MB21D1` and STING is `TMEM173`, so exact
string matching silently drops genes that are present. Aliases are resolved where reference symbols
enter, and matching is reported as a three-way ledger — matched, matched-via-alias, genuinely
absent — never as a pass/fail count.

## Reproducing

The scripts run in numeric order inside the provided container, reading paths and thresholds from
`02_analysis/config/analysis_config.yaml`. Three things are worth knowing before starting:

- **The raw data is not redistributed here.** `00_data/<accession>/download.sh` fetches it from GEO;
  the manifest and sample table record what should arrive.
- **Figures are not tracked.** `.png` and `.pdf` are ignored repo-wide. Each stage's tables and its
  `README.md` are the tracked record, and every figure regenerates from the script its caption
  names.
- **`01_modules/` is empty in a source archive.** The two toolkits under it are git submodules, and
  GitHub's release tarballs do not include submodule contents. Clone with `--recursive` to get them.

## Data availability

Both accessions are public: [GSE160097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160097)
(Mijnheer / Lutter et al. 2021, *Nature Communications*, doi 10.1038/s41467-021-22975-7) and
[GSE161426](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161426). Every external gene set,
lens and network is listed with its version, derivation and citation in `03_results/README.md`, so
each is traceable to its source rather than to this repository.

## License

Two licenses, split by what the file is:

- **Code** — everything under `02_analysis/`, `.devcontainer/`, `config.py` and the scripts anywhere
  in the tree: [MIT](LICENSE).
- **Results and prose** — the tables, figures and README text under `03_results/`, and the
  documentation: [CC BY 4.0](LICENSE-CC-BY-4.0.txt).

Reuse of either requires attribution. The underlying sequence data carries the terms of its GEO
deposition, not these.
