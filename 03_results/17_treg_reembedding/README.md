# 17_treg_reembedding: artifact captions

A second map of the sorted JIA synovial-fluid/blood T-cell compartment (GSE160097), computed on
the Treg gate alone. The published embedding lays out all three sorted lineages in one geometry,
so the variance separating Treg from Tcon from CD8 sets the axes and Treg substructure occupies
what room is left; here the 27,175 cells of the `coarse_label == "Treg"` gate get the whole
canvas. Per-cell scores are joined in from `03_results/interactive/16_narrative_embedding.parquet`
on barcode and are not recomputed, so every map can be coloured by identical values and a
difference between them is a difference of layout.

**Recipe.** The six calls of the full-object embedding, in the same order, on the subset:
`highly_variable_genes(n_top_genes=3000)` -> subset to those genes -> `scale(max_value=10)` ->
`pca(n_comps=30)` -> `neighbors(n_pcs=30)` -> `umap`, with `random_state = 0`. `hvg_n_top` and
`n_pcs` are read from the same `thresholds:` block the full-object run read, and `n_neighbors` is
left at the scanpy default because the full-object call left it there, so the maps differ by
which cells went in and by the batch correction alone. Expression is the log-normalised `X` of
the frozen annotation checkpoint; per-cell normalisation does not depend on which other cells are
present, so subsetting does not require re-deriving it.

**Two coordinate pairs, one PCA.** The uncorrected recipe applies no batch correction, and on this
subset it resolves donor-by-tissue sample of origin (66.1% same-donor neighbours against 14.6%
expected). So the neighbours-and-UMAP tail runs twice off the same PCA: once on the raw PCA, once
on the PCA after Harmony over `donor` at harmonypy standard settings, seeded, converged in 5
iterations. `x`/`y` are the Harmony-corrected coordinates and are the pair to draw.
`x_uncorrected`/`y_uncorrected` are the uncorrected pair, carried for comparison on the same
cells. Harmony is a visualisation aid under the umbrella embeddings guardrail and supports no
claim.

**Verdict on the correction, at k = 30.** Same-donor neighbours fall from 0.661 uncorrected to
0.201 corrected, against a 0.146 chance level, so the corrected map sits at 1.38x chance
(excess index 0.064) and below the full-object map's 0.420. Every donor drops, and the residual
concentrates in one donor (JIA_patient_5, 0.271). Same-tissue neighbours go 0.975 to 0.923
against 0.500 expected, so removing almost all donor structure leaves the synovial-fluid/blood
split at 1.85x chance: donor and tissue are separable in this representation, even though the two
are nested at the sample level. The corrected coordinates are usable for reading Treg
substructure. A small tight cluster still deserves a donor-composition check, and Harmony
reshapes the space it corrects, so this map remains a map.

**Tier.** Annotation / visualisation only. An embedding places cells; it tests nothing. Nothing
here is pooled with the donor-level pseudobulk spine and this stage writes no effect-size row.
Confirmatory claims are carried by donor-level pseudobulk differential expression within the
frozen labels.

## tables/treg_reembedding_mixing.csv

Harmony over `donor` brings the Treg-only map from 0.661 same-donor neighbours to 0.201 against a
0.146 chance level, which is nearer chance than the full-object map restricted to the same 27,175
Tregs (0.420) and holds for all seven donors individually, while same-tissue neighbours fall only
0.975 to 0.923 against 0.500 expected, so the synovial-fluid/blood separation survives the
removal of nearly all donor structure, and the uncorrected map's donor dominance is a property of
the recipe that omits batch correction.

**How to read:** One row per (`embedding` x `space` x `grouping_key` x `group`).
`embedding` names which map was measured: `treg_only_harmony` is this stage's corrected map
(`x`/`y` in the parquet), `treg_only` is the same cells and same PCA without correction
(`x_uncorrected`/`y_uncorrected`), `full_object_restricted` is the published
`X_umap_unsupervised` with the rows restricted to the same Treg cells. The two `_latent` rows
measure the 30-dimensional representations the two subset layouts came from, `pca_30` and
`harmony_30`. `space` says which of those was searched; the three `umap_2d` rows are the
head-to-head, because the published embedding has no stored latent space in the annotation
checkpoint and 2D is also the space a widget draws.
`grouping_key` is the factor whose neighbourhood composition is measured; `group` is one factor
value, with `_all_` the cell-weighted overall row.
`observed_same_frac` is the mean over cells of the fraction of that cell's `k` = 30 nearest
neighbours carrying the same factor value. `expected_same_frac` is what a random relabelling
would give, `(n_group - 1) / (N - 1)` for a cell of that group, so unequal group sizes are
already absorbed and the two columns are directly comparable. `excess_over_chance` is
`(observed - expected) / (1 - expected)`: 0 means neighbourhoods look like the dataset
composition, 1 means every neighbour shares the factor value. For `donor`, lower is better,
since donor structure is what the correction is meant to remove; for `tissue`, a value that
stays high after donor correction says the two factors are separable.
Donor and tissue are nested in this design: each GSM is one donor x one tissue x one sorted
population, so a same-sample neighbour is same-donor and same-tissue at once, which is why the
`tissue` rows are read for whether they collapse and not for their absolute size.
`evidence_tier` reads `annotation_embedding` throughout; this table is a property of a layout and
carries no biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `mixing_table` | `mixing_k = 30`, `mixing_keys = [donor, tissue]`, `harmony_batch_key = donor`, `reference_embedding_key = X_umap_unsupervised`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/objects/17_treg_reembedding.parquet` |

## tables/treg_reembedding_manifest.csv

Every parameter of the Treg-only maps is the parameter the full-object embedding used, sourced
from the same config keys the full-object run read (`hvg_n_top = 3000` and `n_pcs = 30` from
`thresholds:`, `scale(max_value=10)` and the scanpy default `n_neighbors` from the recipe), the
correction adds only `harmony_batch_key = donor` at the harmonypy default `max_iter_harmony = 10`,
the subset holds exactly the 27,175 cells the frozen `coarse_label == "Treg"` gate declares, and
all 16 score columns of the narrative substrate are carried across the barcode join.

**How to read:** One row per recorded parameter: `parameter`, the `value` used in the run, and
`source`, which says where that value came from. `thresholds.*` means the shared config block
the full-object embedding also read, so those two rows come from one place for every map;
`treg_reembedding.*` means this stage's own config block; `recipe` means a value fixed by the
embedding recipe being reproduced; `derived` means the run measured it.
`corrected_coord_columns` and `uncorrected_coord_columns` say which parquet columns hold which
variant, so a consumer does not have to infer it. `harmony_max_iter` is the harmonypy default,
stated here because the correction was run once at standard settings and reported as it came.
`coordinate_state` reads `recomputed` when the embedding ran and `checkpoint_reused` when the
coordinate cache under `03_results/objects/` was accepted, which requires the cache to cover
exactly the current subset's barcodes and to carry both representations and both coordinate
pairs. `n_cells_embedded` against `n_cells_expected` is a reproduction check, and a disagreement
aborts the run. No biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_manifest` | `subset_key = coarse_label`, `subset_value = Treg`, `expected_n_cells = 27175`, `random_seed = 0`, `harmony_batch_key = donor`, `harmony_max_iter = 10`, `hvg_n_top = 3000`, `n_pcs = 30` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |

## ../interactive/17_treg_reembedding.parquet

The per-cell substrate for a Treg-only companion widget: 27,175 rows, one per sorted Treg cell,
carrying both coordinate pairs alongside the identical 16 score columns the full-object map is
coloured by, so a set can be viewed on any of the maps without rescoring. Written alongside the
other per-cell interactive substrates in `03_results/interactive/`, and regenerable (untracked)
as those are.

**How to read:** `barcode` is the annotation object's cell name (10x barcode plus GSM suffix) and
is the join key. `x`/`y` are the Harmony-corrected Treg-only UMAP coordinates and are the pair to
draw; `x_uncorrected`/`y_uncorrected` are the same cells and same PCA without the correction, for
side-by-side comparison. These four are the only columns this stage computes. Each pair is on its
own arbitrary scale, so the two pairs are not comparable coordinate for coordinate with each other
or with the `x`/`y` of `16_narrative_embedding.parquet`. `coarse_label` (constant `Treg` here),
`tissue`, `donor` and `pct_counts_mt` are carried verbatim. The 13 `*_AUCell` columns are
per-cell AUCell of the four mouse-derived human-projected up arms and the nine curated
anchor-independent lenses, bounded in [0, 1] and comparable across tissue **within** a column but
not across columns, since AUCell's scale depends on set size. The three `published_*` columns are
carried from the earlier readout; `published_WT_heat_up` holds a scanpy `score_genes` module
score, so colour the mouse up arm with `WT_heat_up_AUCell`. On `x_uncorrected`/`y_uncorrected`
two thirds of a cell's neighbours share its donor, so a coloured patch there is a statement about
a sample; on `x`/`y` that falls to 1.38x chance. Annotation tier throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_treg_reembedding.py` | `build_substrate` | `subset_key = coarse_label`, `subset_value = Treg`, `harmony_batch_key = donor`, `source_parquet = 16_narrative_embedding.parquet`, `output_parquet = 17_treg_reembedding.parquet` | `03_results/objects/02_annotation.h5ad`, `03_results/interactive/16_narrative_embedding.parquet` |
