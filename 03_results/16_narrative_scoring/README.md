# 16_narrative_scoring: artifact captions

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
Against a genuine AUCell reference for the same mouse arm, the column computed here
reproduces at Pearson and Spearman r = 1.000000 over all 99,915 cells, while
`published_WT_heat_up` reaches r = 0.755 because it is a different metric.

**Scale guard.** Three of the sets scored here already had a per-cell column, so re-deriving
them checks the scale of the whole substrate. `assert_seam_reproduces` in
`02_analysis/scripts/16_narrative_scoring.py` holds every same-metric comparison to
r = 1 within 1e-06 and halts the run otherwise, so a substrate that stopped matching its
references could not reach this directory. It is a run-time assertion and produces no
artifact here: the comparison carries no biological reading, and nothing in this stage is
read from it.

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

## figures/_overview/umap_full_reference.png

One sampled frame of the frozen 99,915-cell sorted JIA T-cell map
coloured three ways, by tissue of origin, by the frozen FACS sort gate
and by donor, establishing the layout that every score colouring of
this substrate is read against: all seven JIA donors contribute cells
to both the synovial-fluid and the paired-blood side of the map, and
the sort gates occupy largely distinct territory within each tissue.

**How to read:** Three panels over ONE sampled frame of the same cells at the same
coordinates, so a cell sits in the same place in all three. Left is
tissue of origin, synovial fluid in vermillion and paired blood in
blue. Middle is the frozen FACS sort gate the compartment is built on,
Treg, Tcon and CD8. Right is donor, one hue per JIA patient. Points
are drawn in shuffled order so overlapping groups paint evenly, the
axes are UMAP coordinates without units, and all three panels share
one square bounding box, so the row is comparable panel to panel and
comparable to the score figures beside it. 60,000 of the 99,915 cells
are drawn with a fixed seed, and the source table gives the full per
cell state, tissue and donor counts next to the counts drawn. This is
annotation. Claims in this compartment rest on donor-level pseudobulk
differential expression within these frozen cell states, meta-analysed
as effect sizes with confidence intervals.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `figure_reference` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, point_size = 1.6, colours = colors.okabe_ito (tissue vermillion/blue, state green/orange/pink, donor 7 hues)` | `03_results/interactive/16_narrative_embedding.parquet` |

## figures/_overview/umap_full_arms.png

All three mouse 39 °C-derived up arms colour the synovial-fluid side
of the map brighter than the paired-blood side in every frozen sort
gate, and the per-cell AUCell means behind that colouring put the
WT_heat_up shift at Treg 0.0112 to 0.0190 in Treg, Tcon 0.0114 to
0.0183 in Tcon and CD8 0.0137 to 0.0179 in CD8, so the colouring
tracks tissue across the whole sorted compartment and reads as gate-
shared.

**How to read:** Three panels over the SAME sampled frame and square bounding box as
the reference figure. Panel titles carry the set identifier and its
size. WT_heat_up is the up arm of the mouse WT iTreg 39 versus 37 °C
contrast, 199 human symbols. KO_heat_up is the same contrast in cGAS-
knockout iTregs, 218 symbols. Interaction_up is the mouse genotype by
temperature up arm, 7 symbols, small enough that one gene moves the
score, so read it for location and treat its spread as noise.
WT_heat_up and KO_heat_up share 182 genes and are drawn on ONE colour
scale, 0.0023 to 0.0336, so the two panels can be compared pixel for
pixel. Interaction_up spans a range an order of magnitude wider, so a
common scale would flatten both, and it carries its own bar. Colour is
per-cell AUCell with the highest-scoring cells drawn last, and the
limits are the 2nd and 98th percentile of the values on that scale.
AUCell is bounded in [0, 1] and scales with set size, so the source
table carries mean, median, standard deviation and cell and donor
counts for any comparison the colour cannot make. Cells are pooled
across donors, making a tissue difference here pseudoreplicated and
descriptive. Claims in this compartment rest on donor-level pseudobulk
differential expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `figure_arms` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, point_size = 1.6, cmap = viridis, clip_percentiles = [2, 98], columns = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, shared_scale = WT_heat_up + KO_heat_up at 0.0023-0.0336` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_full_programs.png

The curated hypoxia lens and the generic type-I interferon axis colour
the synovial-fluid side brighter in all three sort gates, while the 21
published IFN-independent STING-activation genes sit far lower in the
Treg gate than in Tcon or CD8 in both tissues (per-cell AUCell mean
0.0186 in Treg blood against 0.0532 in Tcon and 0.0647 in CD8), so
that panel reports a sort-gate difference alongside a tissue one.

**How to read:** Three panels over the SAME sampled frame and bounding box as the
mouse-arm figure, on the same sequential colormap. Panel titles carry
the set identifier and its size. Each set here is curated, versioned
and derived without reference to the mouse anchor: HALLMARK_HYPOXIA
from MSigDB Hallmark, sting_specific_published the 21 published IFN-
independent STING-activation genes, and ifn_generic_axis a 200-gene
generic type-I interferon axis of which 116 symbols match this object,
the thinnest intersection in the panel. Colour is per-cell AUCell,
clipped to the 2nd and 98th percentile within each panel, highest-
scoring cells drawn last. These three sets are unrelated to each other
and their ranges differ, so each keeps its own scale and brightness
compares tissues within a panel while the source table carries the
cross-panel numbers. Hypoxia and temperature are both imposed by the
inflamed joint and stay entangled in cross-sectional human data, so
this hypoxia panel is one lens on that niche. Claims in this
compartment rest on donor-level pseudobulk differential expression
within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `figure_programs` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, point_size = 1.6, cmap = viridis, clip_percentiles = [2, 98], columns = HALLMARK_HYPOXIA_AUCell, sting_specific_published_AUCell, ifn_generic_axis_AUCell` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/arm_score_violins.png

Across the 99,915 sorted cells the median per-cell AUCell of the mouse
39 °C-derived up arm sits higher in synovial fluid than in paired
blood in all three cell states, and KO_heat_up follows the same
pattern, so the per-cell channel shows the same direction the donor-
level pseudobulk carries. The 7-gene interaction arm leaves 33% of
cells at exactly zero, which is what a 7-gene set does when none of
its genes reaches a cell's top-ranked genes, and its synovial-fluid
median exceeds its blood median in all three cell states. Every
distribution here is one vote per cell over 7 donors of unequal cell
yield (7,348 to 19,106 cells), so it carries shape and spread while
ranking the cell states stays with the donor-level panel.

**How to read:** Annotation tier. One panel per mouse-derived up arm, one violin pair
per frozen sort label: warm is synovial fluid, cool is paired
peripheral blood, black line at the median. AUCell is a rank-based
score in 0 to 1, the area under a cell's gene-recovery curve for that
arm, so it is robust to library size and composition. The comparison
the panel supports is the synovial-fluid-versus-blood offset inside
one cell state of one panel. Each panel has its own y axis because the
arms are 199, 218 and 7 genes and AUCell is computed against each
cell's own ranking, so a level in one panel means nothing against a
level in another. Both ends of every y axis carry headroom, so the
score itself is bounded in [0, 1] and the axis is not. Every one of
the 99,915 cells casts one vote, and the 7 donors contributed 7,348 to
19,106 cells each, so a panel-level average follows the donors that
contributed the most cells. Ranking the cell states is the job of the
donor-level pseudobulk panel 03_results/14_unbiased_enrichment/figures
/_overview/arm_nes_by_cell_state.png, where each donor carries one
vote inside a frozen label and the enrichment is tested; this panel
adds the shape and the spread the donor-level aggregate is built from.
The same-stem source table gives the cell count, mean, median,
quartiles and range of all 18 violins. Naming follows how each arm was
derived, from mouse iTreg 37 versus 39 °C contrasts, and the reading
stays correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring_arms_viz.py` | `build_figure` | `metric = AUCell (rank-based, 0 to 1); arms = WT_heat_up 199, KO_heat_up 218, Interaction_up 7 genes; y_pad_frac = 0.04` | `03_results/interactive/16_narrative_embedding.parquet` |

## tables/_overview/umap_full_reference.csv

The 39 populated (cell state x tissue x donor) strata behind the
reference map: all seven donors appear in both tissues, and the three
strata that are absent are single sort gates in one donor arm rather
than a missing donor.

**How to read:** One row per (`coarse_label` x `tissue` x `donor`) stratum present in
the frozen annotation. `n_cells` is the full-object count and
`n_cells_drawn` the count in the sampled frame the figure draws, so
the two say how faithful the drawn frame is to the object.
`frac_of_total` is `n_cells` over 99,915. Absent combinations carry no
row. Annotation tier, no test and no effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `reference_table` | `sample_n = 60000, sample_seed = 0` | `03_results/interactive/16_narrative_embedding.parquet` |

## tables/_overview/umap_full_arms.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_full_arms.png` (WT_heat_up, KO_heat_up,
Interaction_up), one row per cell state and tissue, so the colouring
can be read as numbers.

**How to read:** A restriction of the stage summary table to the sets this figure
draws. One row per (`set_name` x `coarse_label` x `tissue`) with the
mean, median and standard deviation of the per-cell AUCell score and
the cell and donor counts behind it. AUCell is bounded in [0, 1] and
its scale depends on set size, so values compare across tissue within
a `set_name` and the source of a cross-set comparison is the gene
lists, not these means. Cells are pooled across donors, so the unit of
replication is the cell and every tissue difference here is
pseudoreplicated. `evidence_tier` reads `secondary_percell`
throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `score_table` | `rows = 3 sets x 3 frozen cell states x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## tables/_overview/umap_full_programs.csv

Per-cell AUCell summaries of the 3 sets drawn in
`figures/_overview/umap_full_programs.png` (HALLMARK_HYPOXIA,
sting_specific_published, ifn_generic_axis), one row per cell state
and tissue, so the colouring can be read as numbers.

**How to read:** A restriction of the stage summary table to the sets this figure
draws. One row per (`set_name` x `coarse_label` x `tissue`) with the
mean, median and standard deviation of the per-cell AUCell score and
the cell and donor counts behind it. AUCell is bounded in [0, 1] and
its scale depends on set size, so values compare across tissue within
a `set_name` and the source of a cross-set comparison is the gene
lists, not these means. Cells are pooled across donors, so the unit of
replication is the cell and every tissue difference here is
pseudoreplicated. `evidence_tier` reads `secondary_percell`
throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `score_table` | `rows = 3 sets x 3 frozen cell states x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |
