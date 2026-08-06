# 05_scoring — Does the mouse 39 °C signature separate the niche?

The rankings are frozen. This stage hands them the mouse-derived arms and asks the compartment's
central question: **does the mouse 39 °C-derived up arm separate synovial fluid from paired blood
within a frozen sort label?**

It does, in every population. NES **2.5914** in Treg on 120 of 202 arm genes, **2.6826** in Tcon
on 131 and **2.0614** in CD8 on 114, all at FDR below 7e-07
(`tables/gsea_pseudobulk_{treg,tcon,cd8}.csv`). The Treg score sits between the Tcon and CD8
scores, so **the result is pan-T and Tregs are one of the three populations carrying it.**

The down arm complicates the reading and is reported alongside. It reaches NES 1.4322 at FDR
0.0354 in Tcon — the same sign as the up arm — and carries no direction in Treg (1.0386) or CD8
(1.1331). Both arms move the same way, so the pattern is a shared shift. Recapitulating the mouse
contrast would take the two arms apart.

**Two tiers ship here and they are kept apart.** The pre-ranked enrichment on donor-level
pseudobulk is confirmatory and writes the `primary_pseudobulk` rows of
[`../master/effect_sizes_treg_arthritis.csv`](../master/). The per-cell AUCell scores are a
different estimand on a secondary tier: they shadow the same direction and are never pooled with
the enrichment scores.

**Source pinning.** The mouse projection files are read across the compartment boundary under a
SHA-256 pin, so a changed anchor halts the run before it can silently change the JIA result.

---

## Figures

### `figures/_overview/wt_heat_nes_forest.png`

**The confirmatory answer, both arms and all three populations on one axis.**
Ordered dot plot. x, normalised enrichment score for synovial fluid over paired blood; y, arm ×
population. Circles give the up arm and diamonds the down arm, coloured by population; an
asterisk marks FDR below 0.05. Labels give the effective and nominal set sizes and the FDR.

Effective size tracks the NES ordering, so the ordering across populations is a size effect as
much as a biology one. This is donor-level pseudobulk within frozen sort labels, limma-voom then
pre-ranked enrichment, over the six donors present in both arms. The rows carry no interval,
because an enrichment score of this kind has none.
*Source* `tables/gsea_pseudobulk_{treg,tcon,cd8}.csv` ·
`02_analysis/scripts/05_score_signatures_viz.py`.

### `figures/_overview/wt_heat_running_sum_{treg,tcon,cd8}.png` — three panels

**Where along each population's ranking the two arms concentrate.**
One population per figure, three stacked panels. Top, the weighted running enrichment score as
the ranked list is walked from synovial-up to blood-up, so a positive left-shifted peak marks
synovial enrichment. The y range is pinned to [−1, 1], the one range every running sum in this
project uses. Middle, a rug marking each set member's rank. Bottom, the signed moderated-t
ranking the curve was computed on. Warm brown gives the up arm and cool blue the down arm, and
the legend reports effective and nominal size, NES and FDR.

| Population | Up arm | Down arm |
|---|---|---|
| Treg | NES +2.5914, FDR 1e-14, 120 of 202 genes | NES +1.0386, FDR 0.385, 59 of 96 |
| Tcon | NES +2.6826, FDR 1e-16, 131 of 202 | NES +1.4322, FDR 0.035, 64 of 96 |
| CD8 | NES +2.0614, FDR 7e-07, 114 of 202 | NES +1.1331, FDR 0.256, 60 of 96 |

A curve gives the place along one population's own ranking. The cross-population comparison —
whether one population separates further than another — is read off the ordered dot plot.
*Source* `tables/gsea_pseudobulk_{treg,tcon,cd8}.rds` ·
`02_analysis/scripts/05_score_signatures_viz.R`.

### `figures/_overview/score_violins.png`

**The per-cell channel shadows the same direction.**
Each dot is one donor's mean `WT_heat_up` AUCell score for that state × tissue, and the violins
summarise across donors. AUCell is a rank-based score in [0, 1], the area under each cell's
gene-recovery curve for the up set, robust to library size and composition.

Donor-mean scores sit higher in synovial fluid than in paired blood in all three populations.
Read the shift within a population. The absolute level carries no reading, because AUCell's scale
depends on set size. This is a different estimand from the enrichment dot plot and shares no axis
with it. The down arm is omitted because both arms shift the same way in synovial fluid.
*Source* `tables/donor_label_score_means.csv` ·
`02_analysis/scripts/05_score_signatures_viz.py`.

---

## Tables

### `tables/gsea_pseudobulk_{treg,tcon,cd8}.csv` — the primary numbers

One file per sorted population, one row per gene set (`WT_heat_up`, `WT_heat_down`). `nes` is the
normalised enrichment score against that population's signed moderated-t ranking: positive means
enriched toward the synovial-up end. `set_size` counts set members present in that ranked list, so
it varies by population. `padj` is BH across the sets in the file, `core_enrichment` is the
slash-delimited leading edge, and `database = mouse_projection` marks the set's provenance.

`WT_heat_up` reaches NES 2.5915 in Treg (padj 3.23e-14, 119 members matched), 2.6809 in Tcon
(8.09e-17, 130) and 2.0710 in CD8 (3.61e-07, 113). `WT_heat_down` leans positive throughout and
reaches significance in Tcon alone (1.4322, padj 0.0354) against Treg (1.0386, 0.385) and CD8
(1.0943, 0.308).

**These are the confirmatory rows** and they become the `primary_pseudobulk` entries of
`../master/effect_sizes_treg_arthritis.csv`. Keep them apart from the per-cell score means, which
estimate a different quantity on a secondary tier.

### `tables/donor_label_score_means.csv`

39 rows, one per donor × tissue × frozen label. Three strata are absent because those donors lack
that population. Four score columns hold the AUCell and UCell donor means of `WT_heat_up` and
`WT_heat_down`, and `n_cells` is the number of cells averaged into the row.

Donor-mean `WT_heat_up` AUCell runs Treg 0.0193 against 0.0112, Tcon 0.0178 against 0.0114, CD8
0.0178 against 0.0137 — synovial first in each pair — with the largest relative lift, about 1.7×,
in Treg.

Both scorers are unsigned and rank-based, so only the synovial-versus-blood contrast **within** a
population is interpretable. An absolute level compares only within one gene set, and neither
score carries a direction of its own. This is the substrate for the donor-level standardised mean
difference, a different estimand from the enrichment score.

### `tables/runsum_interactive_{treg,tcon,cd8}_WT_heat_{up,down}.csv`

Six files, one per population × gene set: the gene-by-gene substrate the interactive running-sum
widget reads. Every row is one gene at its position in that population's ranked list, and
`running_es` is the weighted running enrichment score computed identically to the plotted curve,
so a widget curve overlays the static figure exactly.

| Column | Meaning |
|---|---|
| `rank` | 1-based position in the ranked list; 1 is most synovial-enriched. |
| `gene` | HGNC symbol at this rank. |
| `stat` | Signed moderated t — the ranking metric; positive is synovial-up. |
| `running_es` | Weighted running enrichment score at this rank. Plot against `rank`. |
| `hit` | TRUE where the gene is a member of this set — a rug tick and a step up in the curve. |
| `leading_edge` | TRUE where the gene is in the core enrichment, the genes driving the score. |
| `gene_set` · `population` · `contrast` | `WT_heat_up` or `WT_heat_down`; `treg` / `tcon` / `cd8`; `SF_vs_PB_<Pop>`. |

The NES and p-value summarising each curve live in the sibling `gsea_pseudobulk_{tag}.csv`.

### `tables/per_cell_scores.csv` · `tables/source_hash_manifest.csv`

`per_cell_scores.csv` holds one row per cell barcode with its donor, tissue, frozen label and the
four AUCell and UCell score columns. The donor means aggregate from it, and so do the
colocalisation analyses in [`../09_heat_hypoxia/`](../09_heat_hypoxia/) and
[`../10_hsr_lens/`](../10_hsr_lens/).

`source_hash_manifest.csv` pins the mouse-anchor projection files this stage reads: one row per
cross-compartment source, with the SHA-256 checked before the file is consumed.
