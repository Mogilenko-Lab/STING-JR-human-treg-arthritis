# 12_treg_localisation — Per-cell score distributions inside the Treg gate

Five signatures scored per cell with AUCell, summarised across the synovial and blood arms inside
the sorted CD4⁺ Treg gate. This stage is a compute resource and publishes no figure. Its two
tables carry the per-cell scores and their per-signature summaries.

Per-cell distributions are secondary and corroborative. The statistical claims about the same
contrast rest on donor-level pseudobulk differential expression in
[`../03_pseudobulk/`](../03_pseudobulk/) and the enrichment computed on it.

## Signature provenance

The five scored signatures are `WT_heat_up`, `WT_heat_up_purged_hypoxia`, `Interaction_up`,
`Interaction_fdrOnly_up` and `HALLMARK_HYPOXIA`.

| Signature | Origin | How the list was derived |
|---|---|---|
| `WT_heat_up`, `Interaction_up`, `Interaction_fdrOnly_up` | **GSE329522**, this project's own mouse anchor. No paper reference recorded. | Bulk RNA-seq of induced regulatory T cells from primary murine splenic CD4⁺ T cells, in a 2×2 design of genotype (WT, cGAS-KO) × temperature (37 °C, 39 °C), 5 biological replicates per group over 20 libraries. Each is one thresholded model contrast — the WT 39 °C-against-37 °C up arm and the heat-by-genotype interaction at two gates — projected to human orthologs with pinned offline babelgene. |
| `WT_heat_up_purged_hypoxia` | The same mouse arm, with its `HALLMARK_HYPOXIA` members deleted. | The 18 overlap genes removed as a plain set difference, leaving 184. |
| `HALLMARK_HYPOXIA` (200) | **MSigDB Hallmark collection H**, *Homo sapiens*, **v2026.1.Hs**, retrieved offline through **msigdbr 26.1.0**. | Frozen one symbol per line under `00_data/references/msigdb_hallmark/` with a validated expected size, and used whole. |

---

## Tables

### `tables/treg_localisation_summary.csv`

One row per signature × tissue arm (`synovial_fluid`, `peripheral_blood`) within
`coarse_label == "Treg"`. Carries the cell and donor counts, the mean, median and interquartile
range of the per-cell score, and both set sizes: the nominal gene-list count and the effective
in-dataset match.

`power_band` classifies the effective size on the project's standard thresholds — `testable` at
15 genes or more, `underpowered_reported` at 5 to 14, `untestable` below 5 — so a thin set stays
in the table, reported with its size.

### `tables/treg_per_cell_scores.csv`

One row per cell barcode: the rank-based AUCell score for each of the five evaluated signatures,
alongside donor, tissue, coarse label and the embedding coordinates. AUCell is unsigned and
bounded in [0, 1], and its scale depends on set size, so values compare across tissue within a
signature.
