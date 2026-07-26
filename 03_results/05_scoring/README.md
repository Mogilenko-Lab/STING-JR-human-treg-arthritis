# 05_scoring — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/wt_heat_nes_forest.png

The primary readout: whether the mouse 39 °C Treg up-program enriches
in JIA SF-vs-PB Tregs, and whether that enrichment is Treg-
preferential over Tcon/CD8.

**How to read:** Points = fgsea NES for WT_heat up (circle) / down (diamond), colored
by population; x=0 dashed = no enrichment; * = FDR < threshold. Read
Treg-preference as: is the Treg-up NES the largest and significant?
Correlative (consistent-with), not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.py` | `main` | `thresholds.gsea_fdr=0.05; gsea_min_size=5` | `03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv` |

## figures/_overview/score_violins.png

Corroborative per-cell view: donor-mean WT_heat_up AUCell score in SF
vs PB across Treg/Tcon/CD8 — is the SF-vs-PB shift Treg-preferential?

**How to read:** Each dot is one donor's mean WT_heat_up AUCell score for that
state×tissue, and the violins summarise across donors. AUCell is a
rank-based score in [0,1], the area under each cell's gene-recovery
curve for the up-set, robust to library size and composition. Read the
RELATIVE SF-vs-PB shift within each population (Treg SF vs Treg PB),
not the absolute level. This is a different estimand from the forest
NES (fgsea on the donor-pseudobulk ranked list). Secondary tier
(percell), NEVER pooled with the pseudobulk NES. Down arm omitted
because up and down co-shift in SF. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.py` | `main` | `signature = WT_heat_up (AUCell, rank-based [0,1])` | `03_results/05_scoring/tables/donor_label_score_means.csv` |

## figures/_overview/wt_heat_running_sum_treg.png

Per-population leading-edge view: where the mouse 39 °C WT_heat up-
and down-programs concentrate along each population's SF-vs-PB
pseudobulk ranking. The Treg up-curve carries the claim; Tcon and CD8
test whether it is Treg-selective.

**How to read:** Top panel = weighted running enrichment score (ES) walking the ranked
list from SF-enriched (left) to PB-enriched (right); a positive,
left-shifted peak = SF enrichment. Middle rug = gene-set member
positions; bottom = the signed Wald ranking metric. Two curves per
panel, same colour in curve and rug: WT_heat up = warm brown, WT_heat
down = cool blue. ES y clamped to [-1, 1] for cross-population
comparability. Display of compute output (clusterProfiler
gseaResult); correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_treg.rds` |

## figures/_overview/wt_heat_running_sum_tcon.png

Per-population leading-edge view: where the mouse 39 °C WT_heat up-
and down-programs concentrate along each population's SF-vs-PB
pseudobulk ranking. The Treg up-curve carries the claim; Tcon and CD8
test whether it is Treg-selective.

**How to read:** Top panel = weighted running enrichment score (ES) walking the ranked
list from SF-enriched (left) to PB-enriched (right); a positive,
left-shifted peak = SF enrichment. Middle rug = gene-set member
positions; bottom = the signed Wald ranking metric. Two curves per
panel, same colour in curve and rug: WT_heat up = warm brown, WT_heat
down = cool blue. ES y clamped to [-1, 1] for cross-population
comparability. Display of compute output (clusterProfiler
gseaResult); correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_tcon.rds` |

## figures/_overview/wt_heat_running_sum_cd8.png

Per-population leading-edge view: where the mouse 39 °C WT_heat up-
and down-programs concentrate along each population's SF-vs-PB
pseudobulk ranking. The Treg up-curve carries the claim; Tcon and CD8
test whether it is Treg-selective.

**How to read:** Top panel = weighted running enrichment score (ES) walking the ranked
list from SF-enriched (left) to PB-enriched (right); a positive,
left-shifted peak = SF enrichment. Middle rug = gene-set member
positions; bottom = the signed Wald ranking metric. Two curves per
panel, same colour in curve and rug: WT_heat up = warm brown, WT_heat
down = cool blue. ES y clamped to [-1, 1] for cross-population
comparability. Display of compute output (clusterProfiler
gseaResult); correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_cd8.rds` |

## tables/runsum_interactive_{treg,tcon,cd8}_WT_heat_{up,down}.csv

Interactive-widget data contract — the exact-genes substrate the report's
interactive running-sum will read (collaborators hover the ranked genes). One tidy
CSV per population × gene set (6 files). Every row is one gene at its position in
that population's SF-vs-PB pseudobulk ranked list; the `running_es` column is the
DOSE weighted running enrichment score computed identically to the plotted curve
(so the widget curve overlays the static figure exactly). Written by the compute
step (`helpers/fgsea_prerank.R`) off the clusterProfiler `gseaResult`.

**Columns:**

| Column | Type | Meaning |
|---|---|---|
| `rank` | int | 1-based position in the ranked list (1 = most SF-enriched; N = most PB-enriched). |
| `gene` | str | HGNC symbol at this rank. |
| `stat` | float | Signed Wald statistic — the ranking metric (positive = SF-up, negative = PB-up). |
| `running_es` | float | Weighted running enrichment score at this rank (DOSE `gseaScores`, exponent 1). The curve to plot vs `rank`. |
| `hit` | bool | TRUE if `gene` is a member of this gene set (a rug tick / step-up in the curve). |
| `leading_edge` | bool | TRUE if `gene` is a core / leading-edge gene (member of the object's `core_enrichment`); the genes driving the ES. |
| `gene_set` | str | `WT_heat_up` or `WT_heat_down`. |
| `population` | str | `treg` / `tcon` / `cd8`. |
| `contrast` | str | `SF_vs_PB_<Pop>`. |

**How to read (widget author):** plot `running_es` against `rank` as a line; draw a rug
at `rank` where `hit` is TRUE; highlight rows where `leading_edge` is TRUE as the ES-driving
core. The ES peak's rank + sign is the enrichment (positive left-shifted peak = SF-enriched).
The NES / p-value summarising each curve live in the sibling `gsea_pseudobulk_{tag}.csv` /
`.rds`. Primary evidence tier = `primary_pseudobulk`; correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/helpers/fgsea_prerank.R` | (top-level) | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` + frozen `WT_heat_{up,down}.txt` |

## tables/gsea_pseudobulk_{treg,tcon,cd8}.csv

The mouse 39 °C `WT_heat_up` program enriches toward the SF end of every
population's donor-pseudobulk ranked list — NES 2.5146 in Treg (padj 9.65e-14, 105
members matched), 2.5717 in Tcon (padj 4.47e-15, 111), 2.0503 in CD8 (padj
8.42e-07, 115) — while `WT_heat_down` reaches significance nowhere (Treg NES 1.011,
padj 0.435; Tcon 1.337, padj 0.074; CD8 1.224, padj 0.148). The axis is up-arm only,
and it is pan-T rather than Treg-exclusive, with Treg and Tcon carrying it more
strongly than CD8.

**How to read:** One CSV per sorted population, one row per gene set (`WT_heat_up`,
`WT_heat_down`). `nes` is the normalized enrichment score of that set against the
population's signed-Wald SF-vs-PB ranked list: positive = enriched toward the SF-up
end of the list, negative = toward PB-up. `set_size` counts set members actually
present in that ranked list, so it varies by population; `padj` is BH FDR across the
sets in the file; `core_enrichment` is the `/`-delimited leading edge;
`database=mouse_projection` marks the set's provenance. These are the PRIMARY,
confirmatory numbers — they become the `primary_pseudobulk` rows of
`03_results/master/effect_sizes_treg_arthritis.csv`. Never pool them with the
per-cell AUCell score means, which estimate a different quantity on a secondary
tier. Correlative (consistent-with), not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures.py` | `run_fgsea` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000; gsea_fdr=0.05; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` + frozen `WT_heat_{up,down}.txt` |

## tables/donor_label_score_means.csv

Donor-mean `WT_heat_up` AUCell sits higher in SF than in paired PB in all three
populations — Treg 0.0193 vs 0.0112, Tcon 0.0178 vs 0.0114, CD8 0.0178 vs 0.0137 —
the largest relative lift (~1.7x) being Treg, so the per-cell view shadows the
pseudobulk NES without independently confirming it.

**How to read:** 39 rows, one per donor x tissue x frozen label; three strata are
absent because those donors lack that population. Four score columns hold the AUCell
and UCell donor means of `WT_heat_up` and `WT_heat_down`. Both scorers are unsigned
and rank-based (robust to depth and composition), so only the SF-vs-PB contrast
*within* a population is interpretable — absolute levels do not compare across gene
sets of different size, and neither score carries a direction of its own. `n_cells`
is the number of cells averaged into the row. SECONDARY annotation tier: this table
is the substrate for the donor-level AUCell standardized mean difference, a
different estimand from the fgsea NES that must never be pooled with it.
Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures.py` | `main` | `percell_score_ncores=8; signature=WT_heat_{up,down} (AUCell + UCell, rank-based [0,1])` | `03_results/objects/02_annotation.h5ad`, `../mouse_anchor/03_results/human_projection/signatures/WT_heat/WT_heat_{up,down}.txt` |

