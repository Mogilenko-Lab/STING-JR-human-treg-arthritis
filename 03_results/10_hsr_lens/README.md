# 10_hsr_lens — What an anchor-independent proteostasis lens returns

The mouse arm enriches. That arm was derived from a 39 °C contrast, so a natural next question is
whether an **independently curated** heat-shock lens returns anything in the same rankings. This
stage asks it.

**The lens is anchor-independent and MSigDB-derived.** `HSR_sensitivity` is the union of three
human MSigDB v2026.1.Hs sets pulled from the offline `msigdbr` 26.1.0 package by exact name with
validated sizes — `REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS` (101),
`REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` (82), `GOBP_RESPONSE_TO_HEAT` (104) —
totalling 176 genes. A per-gene functional taxonomy splits those 176 into `hsf1_core_hsr` (45),
`co_chaperone` (11), `generic_stress` (72), `npc_transport` (30), `thermosensory` (10) and
`upr_er` (8). **`HSR_core` (56) is `hsf1_core_hsr` plus `co_chaperone`**, and it is the term the
figure uses. It is named for the curated categories it was built from.

| Asset | Source | n |
|---|---|---|
| `HSR_sensitivity` | Union of three MSigDB v2026.1.Hs sets via `msigdbr` 26.1.0 | 176 |
| `HSR_core` | Taxonomy categories `hsf1_core_hsr` + `co_chaperone` | 56 |
| `WT_heat_up` | Mouse-anchor 39 °C up arm, projected to human orthologs | 199 |

Three candidate sets were kept out of the union deliberately.
`GOBP_DETECTION_OF_TEMPERATURE_STIMULUS` and its thermoception sibling are thermosensory-neuron
programs with no bearing on T cells, and `HP_FEVER` is a mutation-etiology panel.

**The independence is measurable.** `HSR_core` shares two genes with the 199-gene `WT_heat_up`
(HSPA1A, HSPH1, Jaccard 0.008, tallied in `tables/hsr_wtheatup_overlap.csv`). The lens is a
separate probe, and that separation is the reason for carrying it.

Every run here reads a frozen copy. `freeze_hsr_lens.R` takes the byte-identical lens from the
mouse anchor's `temp_hsr_human_lens.rds`, so the JIA lists and the anchor lists are the same
genes.

## What the lens returns

**A sign flip at trend level.** The curated core points toward synovial fluid in Treg and away
from it in Tcon and CD8: NES +1.4852 at FDR 0.0651, −1.3284 at 0.1688, −1.1296 at 0.4279. Every
population sits above FDR 0.05, so this supplies directional context and a Treg-selective effect
stays untested here.

**The honest ceiling.** Even a clean heat-shock core is proteotoxic-stress-general. The mouse
anchor is the only setting in this project with an experimental 37/39 °C contrast, so this lens
is read correlatively and stops short of naming a temperature driver. Its scores are annotation
tier, firewalled from the confirmatory `WT_heat` spine, and stay out of
`effect_sizes_treg_arthritis.csv`.

---

## Figures

### `figures/_overview/hsr_core_running_sum.png`

**The curated heat-shock core along all three rankings, on one axis.**
Three stacked panels sharing an x axis of **fractional rank** — each gene's position in its own
population's ranked list as a fraction of that list's length, most synovial-up at 0 and most
blood-up at 1 — because the three rankings differ in length.

Top panel: the weighted running enrichment score as each list is walked left to right, so a
positive left-shifted excursion is synovial enrichment. Its y range is pinned to [−1.0, 1.0], the
one range every running sum in this project uses, so a curve's height means the same thing here
as anywhere else.

Middle panel: where each population's core genes sit in its ranking, one labelled row per
population in matching colour. Bottom panel: the ranked moderated t each curve was computed on,
which shows how much signal each rank carries and where the three rankings cross zero — the
assumption the shared fractional axis rests on.

Legend labels carry each NES and FDR, the effective size against the 56-gene nominal set, and its
testability band. 44 of 56 genes are testable in every ranking. The Treg trace is a trend at FDR
0.065.
*Source* `tables/_overview/hsr_core_running_sum.csv` and the `runsum_interactive_*` traces ·
`02_analysis/scripts/10_hsr_lens_viz.py`.

---

## Tables

### `tables/hsr_lens_nes.csv` — the summary

One row per population and HSR term. Positive `nes` means the term is enriched toward the
synovial-up end. `padj` is the FDR, and `leading_edge` is semicolon-delimited. `evidence_tier`
reads `secondary_annotation` throughout, so the tier statement travels in the data itself.

This is the table that carries the sign flip: `HSR_core` positive in Treg, negative in Tcon and
CD8, with no population below FDR 0.05.

### `tables/hsr_gsea_{treg,tcon,cd8}.csv` · `.rds`

The per-population enrichment output. Positive NES means enrichment toward synovial-up genes,
`padj` is the FDR, `core_enrichment` lists the leading edge. The `.rds` files preserve the
`gseaResult` objects, so a running sum can be reconstructed exactly if a later display needs it.
They are compute substrate, and their numbers read through the sibling CSV and the summary above.

Each population answers its own question: Treg asks whether a clean proteostasis signal survives
where `WT_heat_up` enriched, and Tcon and CD8 ask whether any such signal is Treg-specific or
pan-T.

### `tables/hsr_wtheatup_overlap.csv`

The independence check. `n_a` and `n_b` are the sizes of `WT_heat_up` and the HSR term,
`n_intersect` and `jaccard` quantify the direct overlap, and `genes_intersect` lists the shared
symbols. A gene-list annotation check, carrying no enrichment statistic.

Two genes shared with `HSR_core` at Jaccard 0.008 is what makes the lens a separate probe.

### `tables/hsr_colocalization.csv`

Rows stratified by population, HSR term, correlation level and method. `level = cell` uses
individual synovial cells, and `level = donor_sf_mean` correlates donor-level synovial means.
Positive `r` means higher `WT_heat_up` AUCell tends to coincide with higher HSR AUCell.

A low cell-level `r` means the empirical arm and the curated lens label different cells, which is
consistent with activation carrying much of `WT_heat_up`. Secondary per-cell tier.

### `tables/per_cell_hsr_scores.csv`

One row per cell barcode, with donor, tissue and frozen coarse label alongside four score
columns: `HSR_core_AUCell`, `HSR_core_UCell`, `HSR_sensitivity_AUCell`,
`HSR_sensitivity_UCell`. Both scorers are unsigned and rank-based, so a higher value means the
set is more represented in that cell's expression ranking. This is the substrate the
colocalisation table reads, and it allows the curated lens to be compared with `WT_heat_up`
without selecting cells on either score.

### The remaining files

| File | What it holds |
|---|---|
| `tables/_signatures_hsr/HSR_{core,sensitivity}.txt` | The stage-local copies handed to the enrichment engine — one sorted HGNC symbol per line. Inputs, kept so the command is self-contained. |
| `tables/runsum_interactive_hsr_gsea_{treg,tcon,cd8}_HSR_{core,sensitivity}.csv` | Six running-sum substrates in the shared schema: one row per ranked gene with `running_es`, `hit` and `leading_edge`. |
| `tables/_overview/hsr_core_running_sum.csv` | The figure's same-stem source: per population, the NES, FDR, testable set size and ranked-list length its trace walks, which is why the three traces end at slightly different x. |
| `tables/source_hash_manifest.csv` | The SHA-256 pin on the mouse-anchor projection this stage reads for the overlap table. A changed anchor stops the stage. |
