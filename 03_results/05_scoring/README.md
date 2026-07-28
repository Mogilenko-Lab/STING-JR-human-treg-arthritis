# 05_scoring — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/wt_heat_nes_forest.png

The mouse 39 °C-derived up arm separates synovial fluid from paired
blood in every sorted population — NES 2.5915 in Treg (119 of 199
genes ranked), 2.6809 in Tcon (130) and 2.0710 in CD8 (113), all at
FDR below 1e-6 — so the answer to Treg preference is NO: the result is
pan-T with Tcon the largest, and Tregs are in it rather than
privileged in it. The down arm is not silent either, reaching NES
1.4718 at FDR 0.026 in Tcon, the same sign as the up arm, while
carrying no direction in Treg (0.9676) or CD8 (1.0943).

**How to read:** ANSWERS, at the only tier that may: donor-level pseudobulk within
frozen sort labels, limma-voom then fgsea, paired within each of the 7
donors. Points = fgsea NES for the WT_heat up (circle) and down
(diamond) arm, coloured by sorted population; x = 0 dashed = no
enrichment; * = FDR below 0.05. Beside each point is the effective set
size — members present in that population's ranked list — against the
nominal arm size, and the FDR. Read the answer as written on the face,
not as a Treg-preference check: Tcon has the largest up-arm NES, all
three are significant, the result is pan-T. Read the down arm too —
significant in Tcon at the up arm's sign, so the up arm is not the
only informative one. Do NOT read the ordering of the three NES as a
biological ranking; effective set size tracks it across these rows. An
ordered NES dot plot, not a forest: no pseudobulk NES row here carries
an interval. Correlative, not causal.

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

In Treg the up arm reaches NES +2.5915 at FDR 3e-14 with 119 of its
199 genes in the ranked list, and the down arm reaches NES +0.9676 at
FDR 0.512 with 56 of its 94 genes in the ranked list. The curve shows
WHERE along this population's synovial-fluid-versus-blood ranking
each arm concentrates; whether one sorted population separates more
than another is a cross-population comparison and is read off the
ordered NES dot plot, not off this panel.

**How to read:** One population per panel. Shows HOW the confirmatory answer arises
rather than adding one: the same donor-level pseudobulk fgsea result,
one sorted population at a time. Top panel = weighted running
enrichment score walking the ranked list from SF-enriched (left) to
PB-enriched (right); a positive, left-shifted peak = SF enrichment.
Middle rug = member positions; bottom = the signed moderated-t
ranking metric. Same colour in curve and rug: WT_heat up = warm
brown, down = cool blue. Each legend entry carries that arm's
effective set size against its nominal size, its NES and its FDR, so
no enrichment score sits here without the size it was computed on. ES
y clamped to [-1, 1] so heights compare across the three panels, but
the NES ordering ACROSS populations is not read here — the ordered
NES dot plot carries that, and carries the answer that the result is
pan-T rather than Treg-preferential. Display of compute output;
correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_treg.rds` |

## figures/_overview/wt_heat_running_sum_tcon.png

In Tcon the up arm reaches NES +2.6809 at FDR 8e-17 with 130 of its
199 genes in the ranked list, and the down arm reaches NES +1.4718 at
FDR 0.026 with 61 of its 94 genes in the ranked list. The curve shows
WHERE along this population's synovial-fluid-versus-blood ranking
each arm concentrates; whether one sorted population separates more
than another is a cross-population comparison and is read off the
ordered NES dot plot, not off this panel.

**How to read:** One population per panel. Shows HOW the confirmatory answer arises
rather than adding one: the same donor-level pseudobulk fgsea result,
one sorted population at a time. Top panel = weighted running
enrichment score walking the ranked list from SF-enriched (left) to
PB-enriched (right); a positive, left-shifted peak = SF enrichment.
Middle rug = member positions; bottom = the signed moderated-t
ranking metric. Same colour in curve and rug: WT_heat up = warm
brown, down = cool blue. Each legend entry carries that arm's
effective set size against its nominal size, its NES and its FDR, so
no enrichment score sits here without the size it was computed on. ES
y clamped to [-1, 1] so heights compare across the three panels, but
the NES ordering ACROSS populations is not read here — the ordered
NES dot plot carries that, and carries the answer that the result is
pan-T rather than Treg-preferential. Display of compute output;
correlative, not causal.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; running_sum_ylim=[-1,1]; engine=clusterProfiler::GSEA(by=fgsea)` | `03_results/05_scoring/tables/gsea_pseudobulk_tcon.rds` |

## figures/_overview/wt_heat_running_sum_cd8.png

In CD8 the up arm reaches NES +2.0710 at FDR 4e-07 with 113 of its
199 genes in the ranked list, and the down arm reaches NES +1.0943 at
FDR 0.308 with 57 of its 94 genes in the ranked list. The curve shows
WHERE along this population's synovial-fluid-versus-blood ranking
each arm concentrates; whether one sorted population separates more
than another is a cross-population comparison and is read off the
ordered NES dot plot, not off this panel.

**How to read:** One population per panel. Shows HOW the confirmatory answer arises
rather than adding one: the same donor-level pseudobulk fgsea result,
one sorted population at a time. Top panel = weighted running
enrichment score walking the ranked list from SF-enriched (left) to
PB-enriched (right); a positive, left-shifted peak = SF enrichment.
Middle rug = member positions; bottom = the signed moderated-t
ranking metric. Same colour in curve and rug: WT_heat up = warm
brown, down = cool blue. Each legend entry carries that arm's
effective set size against its nominal size, its NES and its FDR, so
no enrichment score sits here without the size it was computed on. ES
y clamped to [-1, 1] so heights compare across the three panels, but
the NES ordering ACROSS populations is not read here — the ordered
NES dot plot carries that, and carries the answer that the result is
pan-T rather than Treg-preferential. Display of compute output;
correlative, not causal.

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
