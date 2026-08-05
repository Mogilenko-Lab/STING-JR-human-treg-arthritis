# 05_scoring: artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The cohort contains 7 JIA donors, of whom 6 span both arms in each analyzed population. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/wt_heat_nes_forest.png

The mouse 39 °C-derived up arm separates synovial fluid from paired
blood in every sorted population: NES 2.5914 in Treg (120 of 202 genes
ranked), 2.6826 in Tcon (131) and 2.0614 in CD8 (114), all at FDR
below 7e-07. The result is pan-T, and Tregs are one of the three
populations carrying it. The down arm reaches NES 1.4322 at FDR 0.035
in Tcon, the same sign as the up arm, and carries no direction in Treg
(1.0386) or CD8 (1.1331).

**How to read:** This is the confirmatory tier: donor-level pseudobulk within frozen
sort labels, limma-voom then fgsea, on the 6 donors present in both
arms. Points are NES for the up (circle) and down (diamond) arms,
coloured by population; the asterisk marks FDR below 0.05. Labels give
effective and nominal set sizes plus FDR. Effective size tracks the
NES ordering, so the ordering is a size effect. Ordered NES dot plot
with FDR encoding; the rows carry no interval. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.py` | `main` | `thresholds.gsea_fdr=0.05; gsea_min_size=5` | `03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv` |

## figures/_overview/score_violins.png

Donor-mean WT_heat_up AUCell sits higher in synovial fluid than in
paired blood in all three sorted populations, so the per-cell channel
shadows the pseudobulk result in the same direction. This
corroborates. A per-cell score is a different estimand on a secondary
tier, and the shift it shows spans all three populations.

**How to read:** This panel corroborates; the confirmatory answer is the pseudobulk NES
dot plot. Each dot is one donor's mean WT_heat_up AUCell score for
that state × tissue, and the violins summarise across donors. AUCell
is a rank-based score in [0, 1], the area under each cell's gene-
recovery curve for the up set, robust to library size and composition.
Read the SF-versus-PB shift within a population (Treg SF against Treg
PB); the absolute level carries no reading. This is a different
estimand from the pseudobulk NES dot plot and shares no axis with it.
The down arm is omitted because up and down co-shift in synovial
fluid. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.py` | `main` | `signature = WT_heat_up (AUCell, rank-based [0,1])` | `03_results/05_scoring/tables/donor_label_score_means.csv` |

## figures/_overview/wt_heat_running_sum_treg.png

In Treg the up arm reaches NES +2.5914 at FDR 1e-14 with 120 of its
202 genes in the ranked list, and the down arm reaches NES +1.0386 at
FDR 0.385 with 59 of its 96 genes in the ranked list. The curve gives
the place along this population's synovial-fluid-versus-blood ranking
where each arm concentrates. The cross-population comparison —
whether one sorted population separates more than another — is read
off the ordered NES dot plot.

**How to read:** One population per panel, showing the donor-pseudobulk fgsea result
behind the confirmatory answer. The top trace walks from SF-enriched
to PB-enriched genes; a positive left peak indicates SF enrichment.
The middle rug marks set members and the bottom shows the signed
moderated-t ranking. Warm brown is the up arm and cool blue the down
arm. Legends report effective and nominal size, NES, and FDR. The
shared [-1, 1] enrichment-score range supports shape comparison. Read
the cross-population result from the ordered NES dot plot, which
establishes the pan-T pattern. Display of compute output;
correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_treg.rds` |

## figures/_overview/wt_heat_running_sum_tcon.png

In Tcon the up arm reaches NES +2.6826 at FDR 1e-16 with 131 of its
202 genes in the ranked list, and the down arm reaches NES +1.4322 at
FDR 0.035 with 64 of its 96 genes in the ranked list. The curve gives
the place along this population's synovial-fluid-versus-blood ranking
where each arm concentrates. The cross-population comparison —
whether one sorted population separates more than another — is read
off the ordered NES dot plot.

**How to read:** One population per panel, showing the donor-pseudobulk fgsea result
behind the confirmatory answer. The top trace walks from SF-enriched
to PB-enriched genes; a positive left peak indicates SF enrichment.
The middle rug marks set members and the bottom shows the signed
moderated-t ranking. Warm brown is the up arm and cool blue the down
arm. Legends report effective and nominal size, NES, and FDR. The
shared [-1, 1] enrichment-score range supports shape comparison. Read
the cross-population result from the ordered NES dot plot, which
establishes the pan-T pattern. Display of compute output;
correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_tcon.rds` |

## figures/_overview/wt_heat_running_sum_cd8.png

In CD8 the up arm reaches NES +2.0614 at FDR 7e-07 with 114 of its
202 genes in the ranked list, and the down arm reaches NES +1.1331 at
FDR 0.256 with 60 of its 96 genes in the ranked list. The curve gives
the place along this population's synovial-fluid-versus-blood ranking
where each arm concentrates. The cross-population comparison —
whether one sorted population separates more than another — is read
off the ordered NES dot plot.

**How to read:** One population per panel, showing the donor-pseudobulk fgsea result
behind the confirmatory answer. The top trace walks from SF-enriched
to PB-enriched genes; a positive left peak indicates SF enrichment.
The middle rug marks set members and the bottom shows the signed
moderated-t ranking. Warm brown is the up arm and cool blue the down
arm. Legends report effective and nominal size, NES, and FDR. The
shared [-1, 1] enrichment-score range supports shape comparison. Read
the cross-population result from the ordered NES dot plot, which
establishes the pan-T pattern. Display of compute output;
correlative.

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
| `stat` | float | Signed moderated t-statistic (limma-voom) — the ranking metric (positive = SF-up, negative = PB-up). |
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

The mouse 39 °C-derived `WT_heat_up` set enriches toward the SF end of every
population's donor-pseudobulk ranked list — NES 2.5915 in Treg (padj 3.23e-14, 119
members matched), 2.6809 in Tcon (padj 8.09e-17, 130), 2.0710 in CD8 (padj
3.61e-07, 113). `WT_heat_down` also leans positive and reaches significance in Tcon
(NES 1.4718, padj 0.026) but not in Treg (0.9676, padj 0.512) or CD8 (1.0943, padj
0.308), so the up arm carries the axis in Treg and CD8 while both arms move together
in Tcon. The enrichment is pan-T rather than Treg-exclusive, with Treg and Tcon
carrying it more strongly than CD8.

**How to read:** One CSV per sorted population, one row per gene set (`WT_heat_up`,
`WT_heat_down`). `nes` is the normalized enrichment score of that set against the
population's signed moderated-t SF-vs-PB ranked list: positive = enriched toward the SF-up
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

## tables/source_hash_manifest.csv

The stage-05 mouse-signature read is pinned to the current mouse-anchor projection files.

**How to read:** Each row is one cross-compartment source the scoring and running-sum scripts may
read. `sha256` is checked before the file is consumed, so a changed mouse projection stops the stage
instead of silently changing the JIA result.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures.py` | `verify_source_hashes()` | pinned SHA-256 | `../mouse_anchor/03_results/human_projection/signatures/WT_heat/` |

