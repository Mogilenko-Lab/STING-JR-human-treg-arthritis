# 16_narrative_scoring — The mouse arms and the curated lenses on one map

One per-cell substrate for the sorted JIA T-cell compartment (GSE160097), built so that the same
99,915-cell embedding can be coloured by an empirical mouse-derived up arm and by a curated
program lens, and the two colourings read side by side.

**The question.** Where on this map does each program sit, and does the mouse 39 °C-derived arm
mark territory that a curated lens leaves unmarked?

**What was computed.** Fourteen gene sets scored per cell with AUCell on the frozen annotation,
summarised per (set × sort gate × tissue), with a manifest recording how many symbols each score
was actually computed over.

**What was drawn.** Three six-panel strips over one sampled frame — the reference layout, the
mouse arms beside three lenses, the cGAS-STING and interferon family beside the inflammation
programs — plus the two violin figures carrying the distributions those colourings compress, and
one stacked patchwork.

**Tier.** Annotation only. A per-cell score localises a program on a map and tests it in no
case. Nothing here pools with the donor-level pseudobulk spine, and this stage writes no
effect-size row.

## The three kinds of set, kept apart

- Four **mouse-derived, human-projected up arms** — `WT_heat_up`, `KO_heat_up`,
  `Interaction_up`, `Interaction_fdrOnly_up` — which depend on the mouse experiment.
- Nine **curated, anchor-independent lenses** — six frozen MSigDB Hallmark programs, the frozen
  `HSR_core` proteostasis lens, the published 21-gene interferon-independent STING signature (de
  Cevins et al. 2023, Table S6), and a 200-gene generic type-I interferon axis.
- One **project-derived lens**, `eTreg_up`, anchor-independent with a weaker provenance than the
  curated nine: derived here from the GSE161426 supplementary log2 matrix by mean difference and
  Welch t over 4 synovial against 14 blood donors, gated at |log2FC| ≥ 1 and p < 0.05 and capped
  at 200 genes. It reads `kind = project_derived_lens` in the manifest and sits behind its own
  rule in every figure that draws it.

**Down arms are out of scope here.** A per-cell colouring answers "where on this map is this
program high", and a continuous colour scale renders an inverted arm as an absence.

## What the map carries, and what correction: none

The embedding coordinates come verbatim from
[`../interactive/08_harvest_readout.parquet`](../interactive/), which inherits the QC-filter
recipe — `highly_variable_genes` → subset → `scale(max_value=10)` → `pca` → `neighbors` → `umap`
— run with no batch key. Harmony, scVI and every other correction stay out of it.

The donor structure that leaves is measured. On the Treg cells of this map, 42.0% of a cell's 30
nearest neighbours share its donor against 14.6% expected
(`../17_treg_reembedding/tables/treg_reembedding_mixing.csv`, rows `full_object_restricted`).
Donor is crossed with tissue by design — all seven patients contribute to both arms — so that
structure sits inside each tissue. The Treg-only map in
[`../17_treg_reembedding/`](../17_treg_reembedding/) is the corrected one.

## Three conventions that make the strips comparable

**One metric.** This stage ships AUCell and nothing else, so two colourings of the same map
differ by gene set alone. Scores are computed on log-normalised expression from the frozen
annotation checkpoint.

**One geometry.** Every figure is a single row of six panels drawn through
`02_analysis/helpers/umap_grid.py`, which fixes panel size, aspect, margins and colourbar slots.
[`../17_treg_reembedding/`](../17_treg_reembedding/) draws its strips through the same module, so
a strip from either stage stacks against a strip from the other with a column reading top to
bottom. Each map's bounding box is padded to the shared aspect, so a UMAP keeps equal aspect
wherever it appears.

**Panel titles give the size that was scored.** The second title line is
`n_genes_found_in_object` — the symbols the score was computed over after alias resolution. That
sits below the nominal set size, so `WT_heat_up` reads 177 where its carried human set holds 205
symbols, and the generic interferon axis reads 119 of 203. Read `frac_found` in the manifest
before reading any colouring.

**Which column carries the mouse arm.** Use `WT_heat_up_AUCell`. The substrate also carries
`published_WT_heat_up`, a stale mean-centred scanpy `score_genes` module score inherited from the
earlier readout, retained so the discrepancy stays visible. Against a genuine AUCell reference
the column computed here reproduces at r = 1.000000 over all 99,915 cells, where
`published_WT_heat_up` reaches r = 0.755, because it is a different metric.

---

## Figures

### `figures/_overview/umap_full_reference.png`

**The layout every score colouring is read against.**
Six panels over one sampled 60,000-cell frame at identical coordinates, so a cell sits in the
same place in all six. Panels 1 and 2 are the annotation this compartment is built on: tissue of
origin, synovial fluid in vermillion and paired blood in blue, then the frozen FACS sort gate.
Panels 3 to 5 are log-normalised FOXP3, IL2RA and CTLA4, joined on barcode, pooled onto one clip
at the 2nd and 98th percentile (0.00 to 2.88) with the highest-expressing cells drawn last, so
one bar in real units serves all three. Panel 6 is donor.

All seven donors contribute cells to both the synovial-fluid and the paired-blood side, the
three sort gates occupy largely distinct territory within each tissue, and the identity genes
corroborate the Treg gate from expression — FOXP3 mean 1.693 in the Treg gate against 0.059 in
Tcon and 0.018 in CD8. These coordinates carry no batch correction, so the donor panel shows
real donor structure at 0.420 same-donor neighbours against 0.146 expected.
*Source* `tables/_overview/umap_full_reference.csv` ·
`02_analysis/scripts/16_narrative_embedding_viz.py`.

### `figures/_overview/umap_full_arms.png`

**Three mouse up arms beside three anchor-independent lenses.**
Six panels on the reference strip's frame and bounding box. The vertical rules carry the
provenance. Left three, anchor-dependent: `WT_heat_up`, `KO_heat_up`, and `Interaction_up` at 7
genes, small enough that one gene moves the score, so read it for location. Middle two, curated
and versioned: `HALLMARK_HYPOXIA` and the activation-free `HSR_core`. Right one: `eTreg_up`,
this compartment's own GSE161426 effector-Treg contrast, ruled off for its provenance.

One bar serves the row and reads 0 to 1: each panel is clipped to its own 2nd and 98th
percentile and rescaled across that clip, so brightness compares tissue within a panel and
carries no meaning between panels. `WT_heat_up` and `KO_heat_up` share 182 genes and one clip
(0.0023 to 0.0328), so those two compare pixel for pixel.

All three mouse arms colour the synovial-fluid side brighter than paired blood in every sort
gate — and so do all three lenses. Over Treg cells the arm runs 0.0111 to 0.0188 in AUCell mean
while curated hypoxia runs 0.0718 to 0.0950 and the curated proteostasis core 0.0989 to 0.1075.
So the arm's tissue colouring is shared with the lenses beside it. Cells pool across donors, so
a tissue difference here is pseudoreplicated. Temperature and hypoxia are both imposed by the
inflamed joint and stay entangled in cross-sectional human data, so the hypoxia panel carries no
HIF claim.
*Source* `tables/_overview/umap_full_arms.csv` ·
`02_analysis/scripts/16_narrative_embedding_viz.py`.

### `figures/_overview/umap_full_programs.png`

**The cGAS-STING and interferon family, ruled off from inflammation.**
Six panels on the arm strip's frame and colormap. Left of the rule: the published 21-gene
interferon-independent STING signature, the 200-gene generic type-I interferon axis carrying the
thinnest intersection in the strip, and `HALLMARK_INTERFERON_ALPHA_RESPONSE`. Right of it:
`HALLMARK_TNFA_SIGNALING_VIA_NFKB`, `HALLMARK_INFLAMMATORY_RESPONSE` and
`HALLMARK_IL2_STAT5_SIGNALING`, the programs the first family has to be distinguished from. One
bar serves the row on the same rescaled 0-to-1 convention.

The reading rule is the rule itself: a synovial-high colouring shared by both families is
generic inflammation, and only a pattern the left family carries alone would be specific to
STING or interferon. The generic interferon axis and all three inflammatory programs brighten
the synovial side in every gate. The published STING panel scores far lower in the Treg gate
than in Tcon or CD8 in both tissues — per-cell AUCell mean 0.0220 in Treg blood against 0.0547
in Tcon and 0.0672 in CD8 — so that panel reports a sort-gate difference alongside a tissue one.

The published STING set is 21 genes and its own IFN-β validation is underpowered at three
donors, so a bright or dim panel there is consistent with STING pathway activity and stops short
of proving it.
*Source* `tables/_overview/umap_full_programs.csv` ·
`02_analysis/scripts/16_narrative_embedding_viz.py`.

### `figures/_overview/arm_score_violins.png`

**The distributions behind the arm strip, panel for panel.**
One panel per gene set, in the panel order of `umap_full_arms.png`, so the two figures lay side
by side. Inside a panel the x axis is the three frozen sort labels and the two violins of a
label are the two tissues: warm synovial fluid, cool paired blood, black line at the median.
Each panel keeps its own y axis, because the sets range from 7 to 195 genes scored.

The grey row under each panel answers what the violins are worst at. It gives Cliff's δ of
synovial fluid against paired blood inside one sort label — the probability that a random
synovial cell outscores a random blood cell of the same label, ties counted as half, rescaled
onto −1 to +1. It is unit-free, so one ±1 axis serves every panel and a label-selective tissue
effect reads as a δ standing away from the other two in its own panel while staying level in the
panels beside it. Behind each pooled δ sit the per-donor δ values.

All six sets score higher in synovial fluid in all three sort labels. The mouse arm separates at
δ 0.623 in Treg, 0.508 in Tcon and 0.329 in CD8 — largest in Treg — and the three lenses beside
it separate the same tissues at least as far, with `eTreg_up` reaching 0.958 in Treg. Every
donor carrying both tissues agrees in sign with the pooled value on that arm. Ranking the sort
labels belongs to the donor-level panel in
[`../14_unbiased_enrichment/`](../14_unbiased_enrichment/), where the same arm reaches NES 2.68
in Tcon against 2.59 in Treg. The 7-gene interaction arm leaves 33% of cells at exactly zero.
*Source* `tables/_overview/arm_score_violins.csv` ·
`02_analysis/scripts/16_narrative_scoring_arms_viz.py`.

### `figures/_overview/program_score_violins.png`

**The distributions behind the program strip, panel for panel.**
Same geometry, same δ row, in the panel order of `umap_full_programs.png`.

All six lenses score higher in synovial fluid in all three sort labels, so a synovial-side
colouring is shared by the cGAS-STING family and by the inflammation programs it has to be told
apart from. Separations span δ 0.160 (`HALLMARK_TNFA_SIGNALING_VIA_NFKB` in CD8) to 0.739
(`HALLMARK_IL2_STAT5_SIGNALING` in Treg). The widest Treg lead over both other labels belongs to
`HALLMARK_IL2_STAT5_SIGNALING` (0.739 against 0.466 and 0.208), a gap of 0.273, which is the
scale of Treg-selectivity this per-cell channel carries at all.

This figure supports one reading its map does not: the zero fraction. The 18-gene published
STING panel sits lowest of the six in the Treg gate (per-cell mean 0.0424 against 0.0809 in Tcon
and 0.0816 in CD8) and leaves 40% of Treg blood cells at exactly zero, drawing as a body pinned
to the axis. A tissue difference resting on that baseline tracks how many cells score anything
at all, and the zero fraction is given per violin in the source table.
*Source* `tables/_overview/program_score_violins.csv` ·
`02_analysis/scripts/16_narrative_scoring_programs_viz.py`.

### `figures/_overview/umap_full_patchwork.png`

**The reference layout above its arm colouring, on one canvas.**
The two strips this stage ships separately, stacked at identical panel size so a column reads
top to bottom. Nothing new is drawn, and both rows hold the identical frame of cells at
identical coordinates. The units differ: the top row is categorical annotation plus
log-normalised expression, the bottom row per-cell AUCell rescaled per panel onto one 0-to-1
bar. Each row's own caption carries its full reading.

The arms brighten the synovial-fluid side of every gate, and the gate FOXP3 marks is not the
gate where they brighten most.
*Source* `tables/_overview/umap_full_patchwork.csv` ·
`02_analysis/scripts/16_narrative_embedding_viz.py`.

---

## Tables

### `tables/narrative_scoring_manifest.csv` — read this before any colouring

One row per scored gene set. `kind` separates the three provenances. `n_genes_in_set` is the
nominal size on disk, `n_genes_found_in_object` the size the score was actually computed over
after symbol matching, and `frac_found` their ratio. `gate` is a power band on the found size —
`testable` at 15 or more, `underpowered_reported` at 5 to 14, `untestable` below 5 — so a thin
set is reported with its size. `source_path` is repo-relative.

Every set intersects the object thickly enough to be read except two: the 7-gene
`Interaction_up` arm falls below the testable floor, and the generic interferon axis matches
only 119 symbols (59%), the thinnest intersection in the panel.

### `tables/narrative_score_summary.csv` — the substrate as numbers

One row per (`set_name` × `coarse_label` × `tissue`): mean, median and standard deviation of the
per-cell AUCell score, with the cell and donor counts behind it. `coarse_label` is the FACS sort
gate. Scores are bounded in [0, 1] and comparable across tissues **within** a `set_name`, and
AUCell's scale depends on set size, so a larger mean across sets means no stronger program.

Cells pool across donors, so the unit of replication is the cell and every tissue difference
here is descriptive. `evidence_tier` reads `secondary_percell` throughout.

Read across the panel, means sit higher in synovial fluid than in paired blood for the mouse arm
and for the curated hypoxia lens alike, and the mouse arm's shift appears in Tcon and CD8 as
well as Treg.

### `tables/_overview/<figure stem>.csv`

Six same-stem sources, one per figure. `umap_full_reference.csv` and `umap_full_patchwork.csv`
carry one row per (`coarse_label` × `tissue` × `donor`) stratum with the drawn and full-object
cell counts and one `mean_<channel>` column per panel. `umap_full_arms.csv` and
`umap_full_programs.csv` restrict the stage summary to the six sets their figure draws. The two
violin tables carry the cell count, mean, median, quartiles, range and zero fraction of all 36
violins together with each pair's δ, its median shift and its donor range.

Read `metric` before comparing columns: log-normalised expression and AUCell are unrelated
scales.

---

The per-cell feed for the reactive review notebook is
[`../interactive/16_narrative_embedding.parquet`](../interactive/), regenerable and untracked.
