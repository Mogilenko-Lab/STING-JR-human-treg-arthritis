# 12_treg_localisation — Per-cell score distributions inside the Treg gate

Five signatures scored per cell with AUCell, summarised across the synovial and blood arms inside
the sorted CD4⁺ Treg gate. This stage is a compute resource and publishes no figure; its two
tables carry the per-cell scores and their per-signature summaries.

Per-cell distributions are secondary and corroborative. The statistical claims about the same
contrast rest on donor-level pseudobulk differential expression in
[`../03_pseudobulk/`](../03_pseudobulk/) and the enrichment computed on it.

---

## Tables

### `tables/treg_localisation_summary.csv`

One row per signature × tissue arm (`synovial_fluid`, `peripheral_blood`) within
`coarse_label == "Treg"`. Carries the cell and donor counts, the mean, median and interquartile
range of the per-cell score, and both set sizes: the nominal gene-list count and the effective
in-dataset match.

`power_band` classifies the effective size on the project's standard thresholds — `testable` at
15 genes or more, `underpowered_reported` at 5 to 14, `untestable` below 5 — so a thin set is
reported with its size rather than dropped.

### `tables/treg_per_cell_scores.csv`

One row per cell barcode: the rank-based AUCell score for each of the five evaluated signatures,
alongside donor, tissue, coarse label and the embedding coordinates. AUCell is unsigned and
bounded in [0, 1], and its scale depends on set size, so values compare across tissue within a
signature.
