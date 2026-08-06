# 09_heat_hypoxia — Is the enrichment reducible to its hypoxia gene content?

The mouse 39 °C up arm enriches toward synovial fluid in every sorted population. This stage asks
one bounded question of that result: **is it reducible to the set's own HALLMARK_HYPOXIA-overlap
gene content?** That is a membership question, and it is answered by deleting those genes and
re-running the same donor-pseudobulk enrichment.

**It survives.** Removing the 18 overlap genes takes 12 to 15 testable genes out of the arm and
costs 0.126 to 0.164 NES — 2.5914 → 2.4271 in Treg, 2.6826 → 2.5565 in Tcon, 2.0614 → 1.9181 in
CD8 — leaving all three significant at FDR ≤ 4.1e-5
(`tables/gene_purge_nes_comparison.csv`).

**What that licenses is narrow, and stating the bound is the point of the stage.** A positive
purged score means the enrichment holds without those genes. Temperature is untested here.
Temperature and hypoxia are jointly imposed by the inflamed joint and stay entangled in
cross-sectional human data, so the causal structure between them stays undetermined and this
stage asserts none. Hypoxia is a transcriptional readout throughout.

One figure carries the question and stops there. The reader sequence then moves to the curated
whole-arm coverage panel in [`../11_heat_decomposition/`](../11_heat_decomposition/), which alone
answers what the set contains.

**Two tables remain as compute resources with no panel of their own.**
`heat_hypoxia_colocalization.csv` carries the per-cell agreement between the heat and hypoxia
scores, read by the reactive review notebook and by the cross-dataset layer.
`leadingedge_composition.csv` carries a model-assigned taxonomy of the genes at the synovial end
of each ranking; its visualisation is withdrawn, because a fraction over genes selected because
they enriched describes that leading edge rather than the arm.

---

## Figures

### `figures/_overview/heat_purge_nes_paired.png`

**The cost of deleting the hypoxia-overlap genes, arm by arm and population by population.**
Six rows, one per population × mouse arm. x, normalised enrichment score, positive toward
synovial fluid. Each row pairs the full set (large diamond) with its purged form (small circle),
and the connecting bar is the NES cost. Warm brown gives the up arm and cool blue the down arm; a
dark outline marks FDR below 0.05. Right-hand text reports effective and nominal set sizes, the
NES cost, and the purged FDR.

**Two gene counts differ and both are given:** 18 genes come out of the frozen set, of which 12 to
15 were present in a ranked list. The Tcon down arm stays significant at the up arm's sign. This
is the confirmatory tier — donor-level pseudobulk within frozen sort labels, limma-voom then
pre-ranked enrichment — and it licenses a membership statement.
*Source* `tables/gsea_{full,purged}_{treg,tcon,cd8}.csv` ·
`02_analysis/scripts/09_heat_hypoxia_viz.py`.

---

## Tables

### `tables/gene_purge_nes_comparison.csv` — the paired answer

One row per sorted population. `NES_full` is the original `WT_heat_up` score and `NES_purged` the
same engine after removing the hypoxia-overlap genes; positive means enriched toward the synovial
end. `genes_removed` echoes the deleted list.

A positive, significant purged score is evidence that the enrichment holds without its
HALLMARK_HYPOXIA-overlap gene content. It says nothing about temperature.

### `tables/gsea_full_{treg,tcon,cd8}.csv` — the unpurged reference

One file per population, two rows each. `set_size` counts signature genes surviving intersection
with the ranked list — 119 / 130 / 113 of 199 up, 56 / 61 / 57 of 94 down — so roughly half to
two thirds of the projected signature is testable here. `core_enrichment` is the slash-separated
leading edge.

`WT_heat_up` reaches NES 2.59 (Treg), 2.68 (Tcon) and 2.07 (CD8), all at FDR ≤ 3.6e-7.
`WT_heat_down` leans positive and reaches significance in Tcon alone (1.47, FDR 0.026) against
Treg (0.97, 0.51) and CD8 (1.09, 0.31). **The positive down arm is the caveat**: both arms move
the same way, so the synovial shift is a shared shift rather than a clean bidirectional
recapitulation of the mouse contrast.

### `tables/gsea_purged_{treg,tcon,cd8}.csv`

The same runs after the purge. `WT_heat_up` reaches NES 2.43 (Treg), 2.55 (Tcon), 1.93 (CD8) at
FDR ≤ 4.1e-5, and effective size falls from 119 / 130 / 113 to 107 / 115 / 101. `WT_heat_down` is
untouched at 0.97 / 1.47 / 1.09, because none of its 94 genes overlaps hypoxia. The contrast tag
distinguishes purged rows, and these feed `gene_purge_nes_comparison.csv`, which owns the paired
comparison.

### `tables/_signatures_full/WT_heat_{up,down}.txt`

The frozen mouse-anchor human-ortholog sets exactly as handed to the enrichment engine — 199 up
genes and 94 down. Plain newline-delimited HGNC symbols, alphabetically ordered, with the
direction carried by the filename. These are inputs rather than results: their value is
provenance, and they are regenerated verbatim from the frozen contract on every run. Diff them
against `_signatures_purged/` to see exactly what the purge removed.

### `tables/_signatures_purged/WT_heat_{up,down}.txt`

181 up genes after dropping the 18 HALLMARK_HYPOXIA members — ADM, ADORA2B, AK4, ANXA2, ATF3,
CCN1, CDKN1A, EGFR, F3, FOSL2, HK2, IER3, P4HA2, PDGFB, PLAUR, SDC4, SERPINE1, TGM2 — which is 9.0%
of the up set. The 94-gene down list is identical to the full one. The purge is a plain set
difference against the 200-gene HALLMARK_HYPOXIA reference applied to both arms, and that it
changes only the up list is itself informative: the hypoxia overlap sits entirely on the
synovial-high side.

### The two compute resources

**`tables/heat_hypoxia_colocalization.csv`** — rows stratified by population, level and
correlation method. `level = cell` uses synovial cells directly; `level = donor_sf_mean`
correlates per-donor synovial means. Positive `r` means a higher `WT_heat_up` score tends to sit
with a higher hypoxia score.

The cell-level correlation is weak, Spearman 0.08 to 0.20, so `WT_heat_up`-high and hypoxia-high
cells are largely distinct cells by this measure. The donor-level correlation rests on six to
seven donors and is effectively unpowered, so its sign supports no reading. Secondary per-cell
tier.

**`tables/leadingedge_composition.csv`** — one row per population. `n_leading_edge` counts the
core-enrichment genes from the full run, 49 to 71 per population, and each is assigned to one
program by the frozen model-assigned taxonomy named in `taxonomy_source`
(`agy_gemini_3.1_pro_2026-07-14`); `n_unclassified` counts genes assigned to nothing.

The taxonomy puts 55% to 61% of those genes in immediate-early or effector-activation categories.
**Those are fractions of leading-edge genes**, and they sit on a different denominator from the
curated whole-arm counts of 5.6% in mouse and 6.0% after human projection. This table was never
used to re-test a leading-edge subset. **Read whole-arm composition from
[`../11_heat_decomposition/`](../11_heat_decomposition/)**, where curated versioned lenses leave
137 of 199 genes unassigned and contain 2 curated-HSR and 1 type-I-interferon gene.

### `tables/runsum_interactive_gsea_{full,purged}_{treg,tcon,cd8}_WT_heat_{up,down}.csv`

Twelve files, the gene-by-gene walk behind every curve, in the shared running-sum schema: `rank`,
`gene`, `stat`, `running_es`, `hit`, `leading_edge`, plus the set, population and contrast keys.

### `tables/_overview/heat_purge_nes_paired.csv` · `tables/source_hash_manifest.csv`

The figure's same-stem source pairs the full and purged score side by side so the cost reads as
one subtraction — six rows, six markers, nothing aggregated. The manifest pins every
cross-compartment source this stage reads by SHA-256; the compute script stops if a pinned file
drifts.
