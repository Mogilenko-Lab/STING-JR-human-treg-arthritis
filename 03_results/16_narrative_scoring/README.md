# 16_narrative_scoring: artifact captions

One per-cell substrate for the sorted JIA synovial-fluid/blood T-cell compartment
(GSE160097), built so the same embedding can be coloured by an empirical mouse-derived up
arm and by a curated program lens and the two colourings read side by side. Fourteen gene
sets are scored, in three kinds the manifest keeps apart:

- four mouse-derived, human-projected **up arms** — `WT_heat_up`, `KO_heat_up`,
  `Interaction_up`, `Interaction_fdrOnly_up` — which depend on the mouse experiment;
- nine **curated, anchor-independent lenses** — six frozen MSigDB Hallmark programs, the
  frozen `HSR_core` proteostasis lens, the 21 published IFN-independent STING genes, and the
  200-gene generic type-I interferon axis — versioned and independent of the anchor;
- one **project-derived lens**, `eTreg_up`, which is anchor-independent and carries a weaker
  provenance than the curated nine: this compartment derived it from the GSE161426
  supplementary log2 matrix by mean difference and Welch t over 4 synovial against 14 blood
  donors, gated at |log2FC| ≥ 1 and p < 0.05 and capped at 200 genes. It reads `kind =
  project_derived_lens` in the manifest and sits behind its own rule in every figure that
  draws it, so a reader never takes it for a curated set.

**Where this map comes from, and what correction it carries: none.** The embedding coordinates
are taken verbatim from `03_results/interactive/08_harvest_readout.parquet`, which inherits the
`01_qc_filter.py` recipe — `highly_variable_genes` → subset → `scale(max_value=10)` → `pca` →
`neighbors` → `umap` — run with no batch key. Harmony, scVI and every other correction stay
out of it. The donor structure that leaves is measured, not assumed: on the Treg cells of this
map, 42.0% of a cell's 30 nearest neighbours share its donor against 14.6% expected from the
donor proportions (`../17_treg_reembedding/tables/treg_reembedding_mixing.csv`, rows
`full_object_restricted`). Donor is crossed with tissue by design — all seven patients
contribute to both the synovial-fluid and the paired-blood arm — so that structure sits inside
each tissue rather than between the two. The Treg-only map in `../17_treg_reembedding/` is the
corrected one: there the same measure reads 66.1% uncorrected and 20.1% after Harmony over
donor, which is why that stage corrects and this one reports.

**One metric: AUCell.** This stage ships AUCell and nothing else, so that two colourings of
the same map differ by gene set and by nothing else. UCell is a by-product of the shared
rank-based scorer and is dropped before the substrate is written. Scores are computed on
log-normalised expression from the frozen annotation checkpoint.

**Six panels to a strip, one geometry.** Every figure here is a single row of six panels on
one canvas, drawn through `02_analysis/helpers/umap_grid.py`. That module fixes the panel size,
the panel aspect, the margins and the colourbar slots, and `../17_treg_reembedding/` draws its
strips through the same module, so a strip from either stage stacks against a strip from the
other with the panels landing at identical size and a column reading top to bottom. Each map's
bounding box is padded out to the shared panel aspect rather than cropped, so a UMAP keeps
equal aspect and reads undistorted wherever it appears.

**Panel titles give the size that was scored.** The second line of each score panel's title is
`n_genes_found_in_object` from the manifest — the count of symbols the score was actually
computed over in this object, after alias resolution. That number sits below the nominal set
size, so `WT_heat_up` reads 177 where its carried human set holds 205 symbols, and
`ifn_generic_axis` reads 119 of 203. Read the manifest for the nominal size and the fraction
found before reading any colouring.

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
interferon axis matches only 119 symbols (59%) — the thinnest intersection in the panel and the
one whose colouring carries the least support per gene.

**How to read:** One row per scored gene set. `kind` separates three provenances: the
anchor-dependent `mouse_derived_arm` sets, the versioned anchor-independent `curated_lens` sets,
and the single `project_derived_lens` (`eTreg_up`), anchor-independent while falling short of
curated or versioned.
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
coloured six ways, establishing the layout every score colouring of
this substrate is read against: all seven JIA donors contribute cells
to both the synovial-fluid and the paired-blood side, the sort gates
occupy largely distinct territory within each tissue, and the three
Treg identity genes corroborate the Treg gate from expression — FOXP3
mean 1.693 in the Treg gate against 0.059 in Tcon and 0.018 in CD8.
The map carries no batch correction, which leaves 0.420 same-donor
neighbours at k = 30 against 0.146 expected.

**How to read:** Six panels over one sampled frame at identical coordinates, so a cell
sits in the same place in all six. Panels 1 and 2 are the annotation
this compartment is built on: tissue of origin, synovial fluid in
vermillion and paired blood in blue, then the frozen FACS sort gate.
Panels 3 to 5 are log-normalised expression of FOXP3, the lineage
transcription factor, IL2RA, the CD25 chain, and CTLA4, the
suppressive effector, joined on barcode from the 07_embedding
substrate, so the gate is checked against expression; most cells read
zero, which is ordinary for single-cell counts, so the three pool onto
ONE clip at the 2nd and 98th percentile, 0.00 to 2.88, with the
highest-expressing cells drawn last. One bar serves all three and it
is in real units, so brightness compares between the genes as well as
within one. Panel 6 is donor. THESE COORDINATES CARRY NO BATCH
CORRECTION, so the donor panel shows real donor structure: 0.420 same-
donor neighbours at k = 30 against 0.146 expected. Donor is crossed
with tissue by design, every patient contributing to both arms, so
that structure sits inside each tissue. The source table gives the
per-stratum counts and marker means. Annotation tier: claims rest on
donor-level pseudobulk differential expression within these frozen
cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `draw_reference_row` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, panel = 3.24 x 3.60 in x 6 columns, point_size = 1.6, colours = colors.okabe_ito (tissue vermillion/blue, state green/orange/pink, donor 7 hues), markers = FOXP3, IL2RA, CTLA4 on viridis, clip_percentiles = [2, 98], marker_scale = one pooled clip at 0.0000-2.8761 over all 3 genes, one bar` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/07_embedding/tables/hook_factor_substrate.parquet, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## figures/_overview/umap_full_arms.png

All three mouse 39 °C-derived up arms colour the synovial-fluid side
of the map brighter than the paired-blood side in every frozen sort
gate, and so do all three anchor-independent lenses ruled off beside
them, which is the point of putting them in one row: the per-cell
AUCell means put the WT_heat_up shift at Treg 0.0111 to 0.0188 in
Treg, Tcon 0.0113 to 0.0180 in Tcon and CD8 0.0135 to 0.0176 in CD8,
while the curated hypoxia lens runs Treg 0.0718 to 0.0950 and the
curated proteostasis core Treg 0.0989 to 0.1075 over the same Treg
cells, so the anchor arm's tissue colouring is not distinctive of the
anchor.

**How to read:** Six panels over the reference strip's frame and bounding box. THE
VERTICAL RULES CARRY THE PROVENANCE. The left three are anchor-
dependent: WT_heat_up is the mouse WT iTreg 39-versus-37 °C up arm in
human projection, KO_heat_up the same contrast in cGAS-knockout
iTregs, Interaction_up the genotype-by- temperature arm at 7 genes,
small enough that one gene moves the score, so read it for location
and treat its spread as noise. The middle two are curated and
versioned: HALLMARK_HYPOXIA and the activation-free HSR_core
proteostasis lens. The right one is eTreg_up, this compartment's own
GSE161426 effector- Treg contrast over 4 synovial against 14 blood
donors, derived for exploration and ruled off for that reason. ONE BAR
SERVES THE ROW, and it reads 0 to 1, not AUCell: each panel is clipped
to its own 2nd and 98th percentile as before and then rescaled across
that clip, so the picture is unchanged and the panels merely stop
needing six separate keys. Brightness therefore compares tissue WITHIN
a panel and says nothing about level BETWEEN panels — the AUCell
limits each panel was rescaled over are in the config line below and
the values themselves in the source table. WT_heat_up and KO_heat_up
share 182 genes and one clip, 0.0023 to 0.0328, so those two do
compare pixel for pixel. Titles give the symbols scored. Cells pool
across donors, so a tissue difference here is pseudoreplicated.
Temperature and hypoxia are both imposed by the inflamed joint and
stay entangled in cross-sectional human data, so the hypoxia panel
carries no HIF claim, and similarity across a rule is a reason to
test. Claims rest on donor- level pseudobulk differential expression
within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `draw_score_row` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, panel = 3.24 x 3.60 in x 6 columns, point_size = 1.6, cmap = viridis, clip_percentiles = [2, 98], columns = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, HALLMARK_HYPOXIA_AUCell, HSR_core_AUCell, eTreg_up_AUCell, shared_scale = WT_heat_up + KO_heat_up at 0.0023-0.0328, colour = rescaled to panel clip onto [0, 1], one bar for the row; AUCell limits rescaled over: WT_heat_up_AUCell 0.0023-0.0328; KO_heat_up_AUCell 0.0023-0.0328; Interaction_up_AUCell 0.0000-0.3526; HALLMARK_HYPOXIA_AUCell 0.0397-0.1165; HSR_core_AUCell 0.0346-0.1562; eTreg_up_AUCell 0.0034-0.0984` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv, 03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## figures/_overview/umap_full_programs.png

Six curated lenses on one map separate two readings that a single
interferon panel would merge: the 21 published IFN-independent STING-
activation genes sit far lower in the Treg gate than in Tcon or CD8 in
both tissues (per-cell AUCell mean 0.0220 in Treg blood against 0.0547
in Tcon and 0.0672 in CD8), so that panel reports a sort-gate
difference alongside a tissue one, while the generic type-I interferon
axis and the three inflammation and activation programs brighten the
synovial-fluid side across all three gates (Treg 0.0172 to 0.0290 and
Treg 0.0503 to 0.0589 over Treg cells).

**How to read:** Six panels on the arm strip's frame and colormap; titles give the
symbols scored. The rule splits two families. Left of it:
sting_specific_published, the 21 published IFN-independent STING-
activation genes; ifn_generic_axis, a 200-gene generic type-I
interferon axis carrying the thinnest intersection in the strip; and
HALLMARK_INTERFERON_ALPHA_RESPONSE. Right of it:
HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_INFLAMMATORY_RESPONSE and
HALLMARK_IL2_STAT5_SIGNALING, the programs the first family has to be
distinguished from. A synovial-high colouring shared by both families
is generic inflammation; only a pattern the left family carries and
the right one lacks would be specific to STING or interferon. ONE BAR
SERVES THE ROW, and it reads 0 to 1, not AUCell: each panel is clipped
per set to the 2nd and 98th percentile with the highest drawn last,
exactly as before, and then rescaled across that clip. Brightness
therefore compares tissue WITHIN a panel and says nothing about level
BETWEEN panels; the limits each panel was rescaled over are in the
config line below and the values in the source table. The six sets
span 0.001-0.068 to 0.073-0.225 in AUCell, which is why a bar in real
units would flatten half of them. The published STING set is 21 genes
and its IFN-β validation in the positive-control compartment is
underpowered at three donors, so a bright or dim panel there is
consistent with STING pathway activity and never proof of it. Cells
pool across donors, so tissue differences are pseudoreplicated. Claims
rest on donor-level pseudobulk differential expression within the
frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `draw_score_row` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, panel = 3.24 x 3.60 in x 6 columns, point_size = 1.6, cmap = viridis, clip_percentiles = [2, 98], columns = sting_specific_published_AUCell, ifn_generic_axis_AUCell, HALLMARK_INTERFERON_ALPHA_RESPONSE_AUCell, HALLMARK_TNFA_SIGNALING_VIA_NFKB_AUCell, HALLMARK_INFLAMMATORY_RESPONSE_AUCell, HALLMARK_IL2_STAT5_SIGNALING_AUCell; per-panel limits, no pooling across sets, colour = rescaled to panel clip onto [0, 1], one bar for the row; AUCell limits rescaled over: sting_specific_published_AUCell 0.0000-0.1429; ifn_generic_axis_AUCell 0.0006-0.0685; HALLMARK_INTERFERON_ALPHA_RESPONSE_AUCell 0.0729-0.2253; HALLMARK_TNFA_SIGNALING_VIA_NFKB_AUCell 0.0296-0.0983; HALLMARK_INFLAMMATORY_RESPONSE_AUCell 0.0342-0.0905; HALLMARK_IL2_STAT5_SIGNALING_AUCell 0.0296-0.0957` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv, 03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## figures/_overview/arm_score_violins.png

Every one of the six sets drawn on the score map scores higher in
synovial fluid than in paired blood in all three sort labels, and the
mouse 39 °C-derived up arm's separation is Cliff's δ 0.623 in Treg,
0.508 in Tcon and 0.329 in CD8 — largest in Treg of the three, and the
three lenses beside it separate the same tissues at least as far
(eTreg_up reaches 0.958 in Treg). Every donor that carries both
tissues agrees in sign with the pooled value on that arm, in all three
sort labels. This per-cell ordering of the sort labels stands on its
own and the donor-level panel is the one that ranks them, where the
same arm reaches NES 2.68 in Tcon against 2.59 in Treg. The 7-gene
Interaction_up arm leaves 33% of cells at exactly zero, which is what
a set that thin does when none of its genes reaches a cell's top-
ranked genes.

**How to read:** Annotation tier. One panel per gene set, in the panel order of
figures/_overview/umap_full_arms.png in this directory, so the two
figures lay side by side and a column of that map has a distribution
here. Inside a panel the x axis is the three frozen sort labels and
the two violins of a label are the two tissues: warm is synovial
fluid, cool is paired peripheral blood, black line at the median.
AUCell is a rank-based score in 0 to 1, the area under a cell's gene-
recovery curve for that set, so it is robust to library size and
composition. A panel title counts the genes the score was really
computed over — the set's genes present in this object after symbol
resolution, which is the same count the map's own panel title carries
and is smaller than the set's nominal size. Each panel keeps its own y
axis because the sets range from 7 to 195 genes scored and AUCell is
computed against each cell's own ranking, so a level compares within a
panel and a shape compares anywhere. Both ends of every y axis carry
headroom, and the score is bounded in [0, 1]. The grey row under each
panel answers the question the violins are worst at, which is whether
a tissue difference is bigger in one sort label than another. It gives
Cliff's δ of synovial fluid against paired blood inside one sort
label: the probability that a randomly drawn synovial cell outscores a
randomly drawn blood cell of the same label, ties counted as half,
rescaled onto -1 to +1. It is unit-free and bounded, so ONE ±1 axis
serves every panel and a difference between sort labels in one panel
is legible against the same difference in any other. A label-selective
tissue effect is a δ that stands away from the other two in its own
panel and does not do so in the panels beside it. Behind each pooled δ
sit the per-donor δ values, each computed inside one donor's own
paired cells; a donor missing one tissue in a sort label contributes
none, which is why the count differs between labels. Read the donor
cloud first: a pooled marker whose donors straddle zero is one donor's
result. Every one of the 99,915 cells casts one vote and the 7 donors
contributed 7,348 to 19,106 cells each, so both rows are
pseudoreplicated and neither carries a p-value or an interval. Ranking
the sort labels, and testing any of this, is the job of the donor-
level panel 03_results/14_unbiased_enrichment/figures/_overview/arm_ne
s_by_cell_state.png, where each donor carries one vote inside a frozen
label. The same-stem source table gives the cell count, mean, median,
quartiles, range and zero fraction of all 36 violins together with
each pair's δ, its median shift and its donor range. Naming follows
how each set was derived and the reading stays correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring_arms_viz.py` | `build_figure` | `metric = AUCell (rank-based, 0 to 1); panels = percell_map_panels.arm_strip = WT_heat_up (177 genes scored), KO_heat_up (195 genes scored), Interaction_up (7 genes scored), HALLMARK_HYPOXIA (178 genes scored), HSR_core (48 genes scored), eTreg_up (173 genes scored); y_pad_frac = 0.04; delta axis = ±1 shared` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv` |

## tables/_overview/umap_full_reference.csv

The populated (cell state x tissue x donor) strata behind the
reference map, with each stratum's mean FOXP3, IL2RA and CTLA4: all
seven donors appear in both tissues, and the three Treg identity genes
are highest in the Treg strata, so the sort gate and the expression
agree.

**How to read:** One row per (`coarse_label` x `tissue` x `donor`) stratum present in
the frozen annotation. `n_cells` is the full-object count and
`n_cells_drawn` the count in the sampled frame the figure draws, so
the two say how faithful the drawn frame is to the object.
`frac_of_total` is `n_cells` over 99,915. Each `mean_<gene>` column is
that stratum's mean log-normalised expression, the numbers behind the
three marker panels; they are means over zero-inflated per-cell
values, so they rank strata and do not estimate a per-cell level.
Absent combinations carry no row. Annotation tier, no test and no
effect size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `stratum_table` | `sample_n = 60000, sample_seed = 0` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## tables/_overview/umap_full_arms.csv

Per-cell AUCell summaries of the 6 sets drawn in
`figures/_overview/umap_full_arms.png` — WT_heat_up, KO_heat_up,
Interaction_up, HALLMARK_HYPOXIA, HSR_core, eTreg_up — one row per
cell state and tissue, so the colouring reads as numbers.

**How to read:** A restriction of the stage summary table to the sets this figure
draws. One row per (`set_name` x `coarse_label` x `tissue`) with the
mean, median and standard deviation of the per-cell AUCell score and
the cell and donor counts behind it. AUCell is bounded in [0, 1] and
its scale depends on set size, so values compare across tissue within
a `set_name`; a cross-set comparison takes the gene lists themselves.
Cells are pooled across donors, so the unit of replication is the cell
and every tissue difference is pseudoreplicated. `evidence_tier` reads
`secondary_percell` throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `score_table` | `rows = 6 sets x 3 frozen cell states x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## tables/_overview/umap_full_programs.csv

Per-cell AUCell summaries of the 6 sets drawn in
`figures/_overview/umap_full_programs.png` — sting_specific_published,
ifn_generic_axis, HALLMARK_INTERFERON_ALPHA_RESPONSE,
HALLMARK_TNFA_SIGNALING_VIA_NFKB, HALLMARK_INFLAMMATORY_RESPONSE,
HALLMARK_IL2_STAT5_SIGNALING — one row per cell state and tissue, so
the colouring reads as numbers.

**How to read:** A restriction of the stage summary table to the sets this figure
draws. One row per (`set_name` x `coarse_label` x `tissue`) with the
mean, median and standard deviation of the per-cell AUCell score and
the cell and donor counts behind it. AUCell is bounded in [0, 1] and
its scale depends on set size, so values compare across tissue within
a `set_name`; a cross-set comparison takes the gene lists themselves.
Cells are pooled across donors, so the unit of replication is the cell
and every tissue difference is pseudoreplicated. `evidence_tier` reads
`secondary_percell` throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `score_table` | `rows = 6 sets x 3 frozen cell states x 2 tissues, metric = AUCell` | `03_results/16_narrative_scoring/tables/narrative_score_summary.csv` |

## figures/_overview/umap_full_patchwork.png

The reference layout and the signature colouring of the same
60,000-cell frame on one canvas, so the sort gate, the Treg identity
genes and the mouse-derived arms are read against each other without
turning a page: the arms brighten the synovial-fluid side of every
gate, and the gate that FOXP3 marks is not the gate where they
brighten most.

**How to read:** The two strips this stage ships separately, `umap_full_reference`
above `umap_full_arms`, on one canvas at identical panel size so a
column reads top to bottom. Nothing new is drawn, and both rows hold
the identical frame of cells at identical coordinates. The rows share
cells and coordinates, and their units differ: the top row is
categorical annotation plus log-normalised expression, the bottom row
per-cell AUCell of a gene set rescaled per panel onto a single 0-to-1
bar. Each row's own caption carries its full reading, the bottom row's
rules carry each signature's provenance, and the source table gives
every panel's channel as a column over the same strata. The
coordinates carry no batch correction, leaving 0.420 same-donor
neighbours at k = 30 against 0.146 expected. Annotation tier
throughout; claims rest on donor- level pseudobulk differential
expression within the frozen cell states.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `draw_reference_row` | `figures.dpi = 300, figures.rasterized_dpi = 600, sample_n = 60000, sample_seed = 0, panel = 3.24 x 3.60 in x 6 columns, point_size = 1.6, rows = 2 x 6, canvas = 24.0 x 12.6 in, top row = tissue, coarse_label, FOXP3, IL2RA, CTLA4, donor; bottom row = WT_heat_up_AUCell, KO_heat_up_AUCell, Interaction_up_AUCell, HALLMARK_HYPOXIA_AUCell, HSR_core_AUCell, eTreg_up_AUCell, rescaled to panel clip onto [0, 1] on one bar; AUCell limits rescaled over: WT_heat_up_AUCell 0.0023-0.0328; KO_heat_up_AUCell 0.0023-0.0328; Interaction_up_AUCell 0.0000-0.3526; HALLMARK_HYPOXIA_AUCell 0.0397-0.1165; HSR_core_AUCell 0.0346-0.1562; eTreg_up_AUCell 0.0034-0.0984` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_score_summary.csv, 03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv, 03_results/07_embedding/tables/hook_factor_substrate.parquet, 03_results/17_treg_reembedding/tables/treg_reembedding_mixing.csv` |

## tables/_overview/umap_full_patchwork.csv

Every channel the stacked layout draws, as one row per (cell state x
tissue x donor) stratum: the cell counts, the three Treg identity gene
means and the six per-cell AUCell means, so all twelve panels are
readable as numbers from one table.

**How to read:** One row per (`coarse_label` x `tissue` x `donor`) stratum, with the
full-object and drawn cell counts and then one `mean_<channel>` column
per panel of the stacked figure. The `mean_<gene>` columns are log-
normalised expression and the `mean_<set>_AUCell` columns are rank-
based scores in [0, 1] whose scale depends on set size, so a value
compares across strata within its own column and never across columns.
Means over cells within a stratum: the unit is the cell, so nothing
here is a donor-level effect. Annotation tier, no test and no effect
size.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_embedding_viz.py` | `stratum_table` | `sample_n = 60000, sample_seed = 0, channels = 3 marker genes + 6 AUCell sets` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/07_embedding/tables/hook_factor_substrate.parquet` |

## figures/_overview/program_score_violins.png

All six lenses score higher in synovial fluid than in paired blood in
all three sort labels, so a synovial-side colouring is shared by the
cGAS-STING family and by the inflammation programs it has to be told
apart from, and the separations span Cliff's δ 0.160
(HALLMARK_TNFA_SIGNALING_VIA_NFKB in CD8) to 0.739
(HALLMARK_IL2_STAT5_SIGNALING in Treg). Reading the tissue separation
by sort label, the widest Treg lead over both other labels belongs to
HALLMARK_IL2_STAT5_SIGNALING (δ 0.739 in Treg against 0.466 and
0.208), a gap of 0.273, which is the scale of Treg-selectivity this
per-cell channel carries at all. The 18-gene published STING panel
sits lowest of the six in the Treg gate (per-cell mean 0.0424 against
0.0809 in Tcon and 0.0816 in CD8) and leaves 40% of Treg blood cells
at exactly zero, so its Treg tissue difference rests on a zero-
inflated baseline.

**How to read:** Annotation tier. One panel per gene set, in the panel order of
figures/_overview/umap_full_programs.png in this directory, so the two
figures lay side by side and a column of that map has a distribution
here. Inside a panel the x axis is the three frozen sort labels and
the two violins of a label are the two tissues: warm is synovial
fluid, cool is paired peripheral blood, black line at the median.
AUCell is a rank-based score in 0 to 1, the area under a cell's gene-
recovery curve for that set, so it is robust to library size and
composition. A panel title counts the genes the score was really
computed over — the set's genes present in this object after symbol
resolution, which is the same count the map's own panel title carries
and is smaller than the set's nominal size. Each panel keeps its own y
axis because the sets range from 18 to 196 genes scored and AUCell is
computed against each cell's own ranking, so a level compares within a
panel and a shape compares anywhere. Both ends of every y axis carry
headroom, and the score is bounded in [0, 1]. The grey row under each
panel answers the question the violins are worst at, which is whether
a tissue difference is bigger in one sort label than another. It gives
Cliff's δ of synovial fluid against paired blood inside one sort
label: the probability that a randomly drawn synovial cell outscores a
randomly drawn blood cell of the same label, ties counted as half,
rescaled onto -1 to +1. It is unit-free and bounded, so ONE ±1 axis
serves every panel and a difference between sort labels in one panel
is legible against the same difference in any other. A label-selective
tissue effect is a δ that stands away from the other two in its own
panel and does not do so in the panels beside it. Behind each pooled δ
sit the per-donor δ values, each computed inside one donor's own
paired cells; a donor missing one tissue in a sort label contributes
none, which is why the count differs between labels. Read the donor
cloud first: a pooled marker whose donors straddle zero is one donor's
result. Every one of the 99,915 cells casts one vote and the 7 donors
contributed 7,348 to 19,106 cells each, so both rows are
pseudoreplicated and neither carries a p-value or an interval. Ranking
the sort labels, and testing any of this, is the job of the donor-
level panel 03_results/14_unbiased_enrichment/figures/_overview/progra
m_nes_by_cell_state.png, where each donor carries one vote inside a
frozen label. The same-stem source table gives the cell count, mean,
median, quartiles, range and zero fraction of all 36 violins together
with each pair's δ, its median shift and its donor range. Naming
follows how each set was derived and the reading stays correlative.
One reading this figure supports and its map does not: the zero
fraction. A thin set leaves many cells at exactly zero, which draws as
a violin body pinned to the axis and is given per violin in the source
table, and a tissue difference resting on that baseline tracks how
many cells score anything at all.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_narrative_scoring_programs_viz.py` | `build_figure` | `metric = AUCell (rank-based, 0 to 1); panels = percell_map_panels.program_strip = sting_specific_published (18 genes scored), ifn_generic_axis (119 genes scored), HALLMARK_INTERFERON_ALPHA_RESPONSE (97 genes scored), HALLMARK_TNFA_SIGNALING_VIA_NFKB (192 genes scored), HALLMARK_INFLAMMATORY_RESPONSE (182 genes scored), HALLMARK_IL2_STAT5_SIGNALING (196 genes scored); y_pad_frac = 0.04; delta axis = ±1 shared` | `03_results/interactive/16_narrative_embedding.parquet, 03_results/16_narrative_scoring/tables/narrative_scoring_manifest.csv` |
