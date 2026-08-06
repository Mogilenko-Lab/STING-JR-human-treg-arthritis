# 01_qc — Which cells survive, and what the stressed pocket is

Adaptive per-library quality control on median-absolute-deviation windows, doublet flagging, and
the first unsupervised embedding. QC retains **99,915 of 108,414 cells, 92.2%**. One library
drops entirely — SF-Treg p5 (GSM4859852), whose median depth is around 14 UMIs — leaving six of
seven donors with paired synovial and blood Tregs for the donor-level contrast.

**The mitochondrial gate is a 40% ceiling alone**, with no per-library MAD applied to %mt, so
stressed high-mitochondrial cells survive by design. That choice needs its own audit, and the
`mthi_*` artifacts are it: they characterise the pocket of such cells, establish that it is a
discrete and reproducibly definable region, and test what it carries.

**What the pocket is.** Two of the eighteen Treg leiden clusters carry a median near 20%
mitochondrial content against 3.8–8.7% in the other sixteen, and they split on effector identity.
Cluster 3 is effector-high and synovial-restricted; cluster 16 is equally mitochondrial and
effector-low. The pocket retains Treg identity, holds real cells at lower depth (median 1,338
genes against a 200-gene floor), and carries no doublet flag. It is 98% synovial fluid and 69%
one donor, and its blood arm is essentially empty, so it supports no paired contrast of its own.
It is therefore retained inside the main Treg pseudobulk and carved out as a separate stratum
nowhere.

Everything `mthi_*` is `secondary_percell` tier: descriptive, and never pooled with the
donor-pseudobulk statistics.

---

## Figures

### `figures/_overview/qc_violins_per_gsm.png`

**Depth and mitochondrial fraction per library.**
One violin per GSM, coloured by sorted population. x, GSM; rows, UMIs, genes (log y) and %mt.
Read for a grossly degraded library; the source table lists per-GSM medians.
*Source* `tables/qc_metrics_per_cell.csv` · `02_analysis/scripts/01_qc_filter_viz.py`.

### `figures/_overview/cells_kept_dropped.png`

**What QC removes, per stratum.**
Stacked bars per stratum, blue kept and red dropped, where a drop is a MAD outlier, a low-gene
cell, a doublet call, or the hard library exclusion. Read it to confirm every synovial and blood
Treg stratum retains enough cells for pseudobulk. The SF-Treg p5 column drops entirely.
*Source* `tables/cells_kept_dropped.csv` · `02_analysis/scripts/01_qc_filter_viz.py`.

### `figures/_overview/unsupervised_umap.png`

**The first embedding, and the sort gates on it.**
Top row places every cell on the unsupervised UMAP coloured by sort population, tissue, donor and
leiden cluster. Bottom row gives FOXP3, IL2RA, CTLA4 and IKZF2 expression on magma. The three
sorted populations occupy largely distinct territory and the four Treg markers concentrate in the
Treg gate, which supports the sort fidelity. The source table is the leiden × population
cross-tab, a contamination check. This map is a usability lens; the biology is read from
donor-level pseudobulk.
*Source* `../objects/01_qc.h5ad` (`X_umap_unsupervised`) ·
`02_analysis/scripts/01_qc_filter_viz.py`.

### `figures/_overview/mthi_cluster_mt_etreg.png`

**Two Treg clusters carry high mitochondrial content, and they differ in effector identity.**
Left panel: one dot per Treg leiden cluster, size scaling with cell count. x, median %mt; y,
median `score_eTreg`. Orange marks the mitochondrial-high effector pocket (leiden 3) and blue the
mitochondrial-high effector-low cluster (leiden 16). Right panel: a per-Treg %mt against eTreg
hexbin on log density with pocket cells overlaid in orange. Both panels label the two clusters
with the numbering `mthi_cluster_enrichment.csv` uses. High mitochondrial content and effector
identity are separable here, so %mt alone defines no pocket.
*Source* `tables/mthi_cluster_enrichment.csv` + `tables/mthi_treg_membership.csv` ·
`02_analysis/scripts/01_qc_mthi_characterize_viz.py`.

### `figures/_overview/mthi_identity_retention.png`

**The pocket retains Treg identity.**
Grouped bars over five canonical Treg markers: orange the mitochondrial-high effector pocket,
blue the mitochondrial-high non-effector cluster, grey normal Treg. Left panel, median
log-normalised expression, with rank-biserial correlation against normal Treg (Mann-Whitney,
BH-FDR; asterisk marks FDR < 0.05). Right panel, the fraction expressing.

IKZF2 (rbc +0.60) and CTLA4 (+0.37) sit above normal Treg, IL2RA and TIGIT are comparable, and
FOXP3 is modestly lower — expressed in 55% of pocket cells against 80% of normal Tregs, which
tracks the lower sequencing depth of high-mitochondrial cells. These cells are Tregs.
*Source* `tables/mthi_identity_retention.csv` ·
`02_analysis/scripts/01_qc_mthi_characterize_viz.py`.

### `figures/_overview/mthi_qc_discrimination.png`

**The pocket holds real cells at lower depth.**
One panel per metric, three bars each: `eff` the effector pocket, `non` the non-effector cluster,
`norm` normal Treg. Rank-biserial correlation against normal Treg prints on each (Mann-Whitney,
all BH-FDR < 0.05), and the red dashed line marks the 200-gene QC floor.

Median gene count in the pocket is 1,338, far above the floor, which is the expected corollary of
a high mitochondrial fraction. `score_apoptosis` (rbc −0.14) and `score_HSP` (−0.22) both sit
below normal Treg. No pocket cell carries a `predicted_doublet` flag. `doublet_score` went
unpopulated in this run, a Scrublet gap, so the doublet evidence rests on the flag.
`n_genes_by_counts` and `total_counts` are objective QC measures; the two score columns are hand
marker modules and are QC-descriptive.
*Source* `tables/mthi_qc_discrimination.csv` ·
`02_analysis/scripts/01_qc_mthi_characterize_viz.py`.

### `figures/_overview/mthi_heat_honesty.png`

**The mouse 39 °C-derived signature is quiet in the pocket.**
Point gives the group median and bar the interquartile range, for `WT_heat_up` on the left and
the balanced `WT_heat_updown` on the right, across the two mitochondrial-high Treg groups, normal
Treg, Tcon and CD8.

The balanced channel is essentially flat: median −0.075 in the effector pocket against −0.060 in
normal Treg. The one-sided `WT_heat_up` channel shifts up, and that shift co-varies with the
effector and depth axis. The same mouse up arm enriches the donor-pseudobulk synovial-versus-blood
contrast at NES 2.59 in Treg, 2.68 in Tcon and 2.07 in CD8, unchanged when these high-mitochondrial
cells were recovered; those values live in [`../05_scoring/`](../05_scoring/).
*Source* `tables/mthi_heat_percell.csv` ·
`02_analysis/scripts/01_qc_mthi_characterize_viz.py`.

### `figures/_overview/mthi_donor_tissue.png`

**The pocket is one tissue and largely one donor.**
Left panel stacks pocket cells per donor, synovial orange and blood blue, with the black dashed
line at the per-stratum cell floor. Right panel places all Tregs on the unsupervised UMAP with
the effector cluster in orange and the non-effector cluster in blue.

The pocket is 98% synovial fluid (284 of 291 cells) and 69% one donor (p6, 202 cells). It clears
the 20-cell floor in three synovial donors and has an essentially empty blood arm — p6
contributes 4 blood cells and no donor reaches the floor — so it supports no paired
synovial-versus-blood pseudobulk contrast of its own.
*Source* `tables/mthi_donor_tissue.csv` + `tables/mthi_treg_membership.csv` ·
`02_analysis/scripts/01_qc_mthi_characterize_viz.py`.

---

## Tables

### `tables/qc_thresholds_per_gsm.csv`

One row per GSM, carrying the window QC applied to it.
`log1p_total_counts_{lo,hi}` and `log1p_n_genes_by_counts_{lo,hi}` give that library's MAD window
(median ± 5 × MAD) on the log1p scale — exponentiate to read them as UMIs and genes — and a cell
outside either window is a `mad_outlier`. `pct_mt_hi` is the effective %mt cutoff and
`pct_mt_ceiling` the hard ceiling; they are equal in every row, because `mt_policy` reads
`ceiling_only` in all forty. `excluded_gsm` marks the hard-dropped library and
`n_flagged_outlier` counts the cells the MAD rule caught before the gene floor, doublet and
exclusion filters union in.

The windows are genuinely per-library: across the 39 retained GSMs the lower UMI bound runs 189
to 2,639 and the upper 11,168 to 104,597. GSM4859852 alone collapses, its window spanning 7 to 27
UMIs and under 30 genes per cell. Across the cohort the MAD rule flags 4,992 of 108,414 cells
(4.6%), between 0.1% and 14.9% within a GSM.

### `tables/mthi_cluster_enrichment.csv`

One row per unsupervised leiden cluster, restricted to sorted Treg.
`fisher_or` / `fisher_p` / `fisher_fdr` test one-sided over-representation of the interactively
saved pocket cells in that cluster, BH across clusters; `frac_saved_captured` is the cluster's
share of all saved Tregs. `is_mthi_cluster` is True at FDR < 0.05 with `median_pct_mt` above twice
the global Treg median; `is_mthi_effector` adds `median_score_eTreg` above the global Treg median.

Cluster 3 (n = 291, median %mt 20.13, median `score_eTreg` +0.0790, OR 1550.3, FDR 0) is
effector-high and captures 77% of the saved pocket cells. Cluster 16 (n = 216, %mt 19.65,
`score_eTreg` −0.0738, OR 55.9, FDR 2.1e-84) is equally mitochondrial and effector-low. A cluster
holding no saved cell carries OR 0 and p 1 by construction; read those rows for their %mt and
depth.

### `tables/mthi_population_rule.csv`

One row per candidate definition of the pocket, testing whether a quantitative rule reproduces
the hand-drawn one. It does: `%mt ≥ P97.5 (10.03)` and `score_eTreg ≥ the within-Treg median
(−0.0031)` selects 340 Tregs, overlapping the cluster definition at Jaccard 0.626 and recovering
0.646 of the 305 saved Tregs, against 291 cells and 0.767 recovery for the cluster rule itself.

`n_cells` is how many sorted Tregs a rule selects; `jaccard_vs_A_effector` is intersection over
union against `A_cluster_effector`, blank on that row; `frac_saved_recovered` is the fraction of
the 305 saved Tregs the rule recaptures. `A_cluster_all_mthi` recovers every saved cell (1.0) at
the lower Jaccard 0.574, which follows from its folding in the effector-low cluster.

### `tables/mthi_heat_percell.csv`

One row per group: the two mitochondrial-high Treg subsets, normal Treg, and the Tcon and CD8
sorts for scale. `WT_heat_up_*` is the one-sided up-gene score and `WT_heat_updown_*` the balanced
up-minus-down score, each as a median with quartiles. `*_rbc_vs_normal` is the rank-biserial
correlation against normal Treg, populated on the `mt_hi_effector` row.

One-sided `WT_heat_up` sits above normal Treg (median −0.036 against −0.060; rbc +0.48, p
4.2e-45) while the balanced channel sits slightly lower (−0.075 against −0.060; rbc −0.17, p
6.0e-7), and both mitochondrial-high groups fall below Tcon (−0.014) and CD8 (−0.037) on that
balanced channel. **Read the balanced channel:** with tens of thousands of cells the p-values are
trivially small, so sign and rbc magnitude carry the meaning. The `caveat` column travels with the
table on purpose.

### The remaining tables

| File | What it holds |
|---|---|
| `qc_metrics_per_cell.csv` | Per-cell UMI count, gene count and %mt, with the keep or drop decision and its reason. |
| `cells_kept_dropped.csv` | Per stratum, cells in and cells out, split by which filter removed them. |
| `mthi_donor_tissue.csv` | Pocket cells per donor × tissue, against the per-stratum floor. |
| `mthi_treg_membership.csv` | Per Treg cell, its leiden cluster and its pocket membership flags. |
| `mthi_identity_retention.csv` · `mthi_qc_discrimination.csv` | Per marker or metric, the three group values with the rank-biserial correlation and its FDR. |
| `_overview/*.csv` | The same-stem source of each figure above. |
