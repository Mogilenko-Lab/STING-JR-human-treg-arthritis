# 16_narrative_scoring — artifact captions

One per-cell substrate for the sorted JIA synovial-fluid/blood T-cell compartment
(GSE160097), built so the same embedding can be coloured by an empirical mouse-derived up
arm and by a curated program lens and the two colourings read side by side. Thirteen gene
sets are scored — four mouse-derived, human-projected **up arms** (`WT_heat_up`,
`KO_heat_up`, `Interaction_up`, `Interaction_fdrOnly_up`) and nine curated,
anchor-independent lenses (six frozen MSigDB Hallmark programs, the frozen `HSR_core`
proteostasis lens, the 21 published IFN-independent STING genes, and the 200-gene generic
type-I interferon axis).

**One metric: AUCell.** This stage ships AUCell and nothing else, so that two colourings of
the same map differ by gene set and by nothing else. UCell is a by-product of the shared
rank-based scorer and is dropped before the substrate is written. Scores are computed on
log-normalised expression from the frozen annotation checkpoint.

**No down arms.** Down arms are deliberately out of scope here: not scored, not carried,
not named. A per-cell colouring answers "where on this map is this program high", and a
down arm inverted onto that question reads as an absence, which a continuous colour scale
cannot render honestly.

**Tier.** Secondary / annotation only. A per-cell score localises a program on a map; it
does not test it. Nothing here is pooled with the donor-level pseudobulk spine, and this
stage writes no effect-size row. Confirmatory claims are carried by donor-level pseudobulk
differential expression.

**Which column to colour the mouse arm with.** Use `WT_heat_up_AUCell`. The substrate also
carries `published_WT_heat_up`, and that column is **not** AUCell — it is a stale
mean-centred scanpy `score_genes` module score inherited by the published per-cell readout,
and it is retained only so the discrepancy stays visible rather than being quietly dropped.
The seam-check panel below is where that is established.

## figures/_overview/narrative_score_seam_check.png

Re-deriving the two genuine AUCell columns of the published per-cell
readout reproduces them exactly (Pearson and Spearman r = 1.000000
over all 99,915 cells), so this substrate sits on the published AUCell
scale; the third comparison reaches only r = 0.755 because the
published `WT_heat_up` column is a stale mean-centred scanpy
score_genes module score rather than AUCell, and the same mouse up arm
reproduces at r = 1.000000 against the canonical AUCell column for
that arm.

**How to read:** One grey point per cell: the published score on x, the score newly
computed in this stage on y. Green statistics and a dashed blue
identity line mark a same-metric comparison (AUCell vs AUCell), where
landing on the line means the scorer reproduced the published column;
vermillion statistics and a vermillion title mark a comparison between
two DIFFERENT metrics, where no identity line applies and the scatter
is expected to spread. Correlations are Pearson and Spearman over the
shared barcodes, with the r >= 0.98 pass floor stated as a verdict in
each box. This is a provenance panel and carries no biological claim;
confirmatory claims are made by donor-level pseudobulk differential
expression.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring_viz.py` | `build_figure` | `percell_score_ncores = 8` | `03_results/interactive/16_narrative_embedding.parquet` |

## tables/_overview/narrative_score_seam_check.csv

The AUCell scorer used here reproduces both genuine AUCell columns of the published per-cell
readout at Pearson and Spearman r = 1.000000 over all 99,915 cells, so this substrate sits on
the published scale; the only comparison that misses the floor is the one whose reference
column is a different metric altogether, and the mouse up arm it concerns reproduces at
r = 1.000000 once compared against a genuine AUCell reference.

**How to read:** One row per re-derivation check: `new_column` computed here,
`reference_column` the published column it is checked against, `reference_source` where that
reference lives. `comparison_kind` carries the row's interpretation. `same_metric` rows compare
AUCell against AUCell and are the only rows that test drift; they must clear `r_floor` (0.98),
and a failure there invalidates every colouring built on the substrate. The single
`cross_metric` row compares AUCell against a scanpy `score_genes` score, so its low r measures a
metric difference, not drift — it is retained to put on record why the mouse arm must be
coloured with `WT_heat_up_AUCell`. `passes_floor` is `pearson_r >= r_floor`, so `False` there is
expected by construction. No biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring.py` | `seam_check` | `percell_score_ncores = 8` | `03_results/interactive/08_harvest_readout.parquet`, `03_results/05_scoring/tables/per_cell_scores.csv` |

## tables/narrative_scoring_manifest.csv

Every scored set intersects the annotation object thickly enough to be read except two: the
7-gene `Interaction_up` arm falls below the 15-gene testable floor, and the 200-gene generic
interferon axis matches only 116 symbols (58%) — the thinnest intersection in the panel and the
one whose colouring carries the least support per gene.

**How to read:** One row per scored gene set. `kind` separates the anchor-dependent
`mouse_derived_arm` sets from the versioned, anchor-independent `curated_lens` sets.
`n_genes_in_set` is the nominal size on disk; `n_genes_found_in_object` is the size the score was
**actually** computed over after matching HGNC symbols against the object; `frac_found` is their
ratio. Check `frac_found` before reading any colouring — a set that barely intersects the object
cannot support one, whatever its nominal size claims. `gate` is a power band on
`n_genes_found_in_object`: `testable` (>= 15), `underpowered_reported` (5-14), `untestable`
(< 5); a set below the floor is reported *with* its size rather than dropped. `source_path` is
repo-relative, so provenance is checkable from the row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring.py` | `build_manifest` | `percell_score_ncores = 8` | `03_results/objects/02_annotation.h5ad` |

## tables/narrative_score_summary.csv

Read across the panel, per-cell AUCell means sit higher in synovial fluid than in paired blood
for the mouse `WT_heat_up` arm and for the curated hypoxia lens alike, and the mouse arm's shift
is not confined to the Treg gate — it appears in Tcon and CD8 too; whether any of that survives
donor-level testing is a question for the pseudobulk spine, not for this table.

**How to read:** One row per (`set_name` x `coarse_label` x `tissue`): mean, median and standard
deviation of that set's per-cell AUCell score, with the cell and donor counts behind it, so the
substrate is legible without loading the parquet. `coarse_label` is the FACS sort gate
(`Treg`/`Tcon`/`CD8`), not a score-derived selection. Scores are bounded in [0, 1] and comparable
across tissues **within** a `set_name`, but **not** across `set_name` — AUCell's scale depends on
set size, so a larger mean is not a stronger program. Cells are pooled across donors, so the unit
of replication is the cell: every SF-versus-PB difference here is pseudoreplicated and purely
descriptive. `evidence_tier` reads `secondary_percell` throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring.py` | `summarise` | `percell_score_ncores = 8` | `03_results/interactive/16_narrative_embedding.parquet` |
