# 01_qc — Quality Control

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood; Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory; rbc =
rank-biserial correlation; MAD = median absolute deviation. The cohort holds 7 JIA donors, of whom 6
span both arms in each analyzed population after QC._

Adaptive per-GSM QC on MAD windows, doublet flagging, and the first unsupervised embedding. The
mitochondrial gate is a 40% ceiling alone, so stressed high-mito cells survive by design; the `mthi_*`
artifacts characterise the pocket of such cells.

## figures/_overview/qc_violins_per_gsm.png

Per-GSM UMI and gene depth is adequate across every stratum and %mito stays low, so no GSM is
grossly degraded.

**How to read:** One violin per GSM on x, coloured by sorted population; rows are UMIs, genes (log
y) and %mito. The source table lists per-GSM medians. QC diagnostic — no biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.qc_n_mads = 5; qc_pct_mt_max = 40.0` | `03_results/01_qc/tables/qc_metrics_per_cell.csv` |

## figures/_overview/cells_kept_dropped.png

QC retains 99,915 of 108,414 cells, 92.2% overall. The SF-Treg p5 library (GSM4859852) is near-empty
at a median around 14 UMIs and drops entirely, leaving 6 of 7 donors with paired SF and PB Tregs for
the donor-level contrast.

**How to read:** Stacked bars per stratum, blue kept and red dropped, where a drop is a MAD outlier,
a low-gene cell, a doublet call, or the hard GSM exclusion. Read it to confirm every SF and PB Treg
stratum retains enough cells for pseudobulk. QC diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.qc_min_genes = 200; scrublet_expected_doublet_rate = 0.06` | `03_results/01_qc/tables/cells_kept_dropped.csv` |

## figures/_overview/unsupervised_umap.png

Sorted Treg, Tcon and CD8 cells occupy largely distinct transcriptomic territory, and FOXP3, IL2RA,
CTLA4 and IKZF2 concentrate in the Treg gate, supporting sort fidelity.

**How to read:** Top row places cells on the unsupervised UMAP coloured by sort population, tissue, donor
and leiden cluster; bottom row gives Treg-marker expression on magma. The source table is the leiden ×
population cross-tab, a contamination check. This UMAP is a usability lens; biology is read from
donor-level pseudobulk.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.hvg_n_top = 2000; n_pcs = 30; leiden_resolution = 1.0` | `03_results/objects/01_qc.h5ad (X_umap_unsupervised)` |

## figures/_overview/mthi_cluster_mt_etreg.png

Two Treg leiden clusters carry high mitochondrial content, a median near 20% against 3.8% to 8.7%
in the remaining sixteen. They split on effector identity: one is the effector-like, SF-restricted
pocket, the other is equally mito-high and effector-low. The pocket is a discrete, reproducibly
defined region of the embedding.

**How to read:** Left panel plots one dot per Treg leiden cluster, size scaling with cell count, x the
median %mt and y the median `score_eTreg`; orange is the mt-hi effector pocket, purple the mt-hi
non-effector cluster. Right panel is a per-Treg %mt against eTreg hexbin on log density, pocket cells
overlaid in green. The panel legend labels the clusters cl6 and cl16; `mthi_cluster_enrichment.csv`
numbers them leiden 3 and 16. `secondary_percell` / EDA tier — never pseudobulk evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.qc_pct_mt_max = 40.0; qc_n_mads_mt = null (ceiling-only)` | `03_results/01_qc/tables/mthi_cluster_enrichment.csv + mthi_treg_membership.csv` |

## figures/_overview/mthi_identity_retention.png

The mt-hi effector pocket retains Treg identity. IKZF2 (rbc +0.60) and CTLA4 (+0.37) sit above
normal Treg, IL2RA and TIGIT are comparable, and FOXP3 is modestly lower (expressed in 55% of
pocket cells against 80% of normal Tregs), which tracks the lower sequencing depth of high-mito
cells. These cells are Tregs.

**How to read:** Grouped bars over 5 canonical Treg markers: mt-hi effector in orange, mt-hi
non-effector in purple, normal Treg in grey. Left panel is median log-normalised expression, where
rbc is the rank-biserial correlation against normal Treg (Mann-Whitney, BH-FDR; asterisk marks FDR
< 0.05). Right panel is the fraction expressing. `secondary_percell` tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `decisions.treg_gate.basis = sorting (refined by FOXP3/IL2RA/CTLA4/IKZF2)` | `03_results/01_qc/tables/mthi_identity_retention.csv` |

## figures/_overview/mthi_qc_discrimination.png

The pocket holds real cells at lower depth: median 1,338 genes, far above the 200-gene QC floor,
which is the expected corollary of a high mitochondrial fraction. `score_apoptosis` (rbc −0.14) and
`score_HSP` (−0.22) both sit below normal Treg. No cell in the pocket is flagged
`predicted_doublet`. `doublet_score` was not populated in this run, a Scrublet gap, so the doublet
evidence rests on the flag alone.

**How to read:** Four bars per metric across groups — eff is the mt-hi effector pocket, non the mt-hi
non-effector cluster, norm normal Treg. rbc is the rank-biserial correlation against normal Treg
(Mann-Whitney, all BH-FDR < 0.05); the red dashed line is the 200-gene QC floor. `n_genes_by_counts` and
`total_counts` are objective QC measures; `score_apoptosis` and `score_HSP` are Tier-3 hand marker
modules, QC-descriptive. `secondary_percell` tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.qc_min_genes = 200; qc_pct_mt_max = 40.0` | `03_results/01_qc/tables/mthi_qc_discrimination.csv` |

## figures/_overview/mthi_heat_honesty.png

The mouse 39 °C-derived signature is quiet in the pocket. The balanced `WT_heat_updown` channel is
essentially flat, median −0.075 in the mt-hi effector pocket against −0.060 in normal Treg. The
one-sided `WT_heat_up` channel shifts up, and that shift co-varies with the effector/depth axis.
This is `secondary_percell` tier, and it is a per-cell descriptive reading. The donor-pseudobulk
enrichment of the same mouse up arm was unchanged when these high-mito cells were recovered; those
NES values and their FDRs live in `03_results/05_scoring/`.

**How to read:** Point is the group median and bar the IQR, for `WT_heat_up` on the left and the
balanced `WT_heat_updown` on the right, across the two mt-hi Treg groups, normal Treg, Tcon and CD8.
`secondary_percell` / EDA tier — descriptive, never pooled with the donor-pseudobulk NES.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `decisions.go_no_go.primary_signature = WT_heat (pseudobulk, NOT per-cell)` | `03_results/01_qc/tables/mthi_heat_percell.csv` |

## figures/_overview/mthi_donor_tissue.png

The pocket is 98% synovial fluid (284 of 291 cells) and 69% one donor (p6, 202 cells). It clears the
20-cell floor in 3 SF donors and has an essentially empty PB arm — p6 contributes 4 PB cells and no
donor reaches the floor — so it cannot support its own paired SF-versus-PB pseudobulk contrast. It
is therefore retained inside the main SF-versus-PB Treg pseudobulk and shown on the embedding, and
it is carved out as its own DE stratum nowhere.

**How to read:** Left panel stacks pocket cells per donor, SF orange and PB blue, with the black
dashed line the per-stratum `min_cells`. Right panel places all Tregs on the unsupervised UMAP with
the mt-hi effector cluster in orange and the non-effector cluster in purple. `secondary_percell`
tier; the UMAP is a usability lens.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.pseudobulk_min_cells = 20; pseudobulk_min_donors = 3` | `03_results/01_qc/tables/mthi_donor_tissue.csv + mthi_treg_membership.csv` |

## tables/qc_thresholds_per_gsm.csv

The mitochondrial gate is the 40% ceiling alone in all 40 GSMs (`mt_policy = ceiling_only`), so
stressed high-mito cells survive QC by design. The depth windows are genuinely per-library: across
the 39 retained GSMs the lower UMI bound runs 189 to 2,639 and the upper bound 11,168 to 104,597.
GSM4859852 alone collapses, its MAD window spanning 7 to 27 UMIs and under 30 genes per cell, which
is why it is the single hard-excluded library. Across the cohort the MAD rule flags 4,992 of 108,414
cells (4.6%), between 0.1% and 14.9% within a GSM.

**How to read:** One row per GSM. `log1p_total_counts_{lo,hi}` and `log1p_n_genes_by_counts_{lo,hi}` give
that GSM's MAD window (median ± 5 × MAD) on the log1p scale — exponentiate to read them as UMIs and genes;
a cell outside either window is a `mad_outlier`. `pct_mt_hi` is the effective %mt cutoff and
`pct_mt_ceiling` the hard ceiling, equal in every row because no per-GSM MAD is applied to %mt.
`excluded_gsm` marks the hard-dropped library, and `n_flagged_outlier` counts the cells the MAD rule
caught in that GSM before the gene floor, doublet and exclusion filters union in. QC diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter.py` | `main` | `thresholds.qc_n_mads = 5; qc_n_mads_mt = null (ceiling-only); qc_pct_mt_max = 40.0; decisions.qc.drop_gsms = [GSM4859852]` | `03_results/objects/00_build.h5ad` |

## tables/mthi_cluster_enrichment.csv

2 of the 18 Treg leiden clusters carry high mitochondrial content, a median near 20% against 3.8% to
8.7% elsewhere, and the two split on effector identity. Cluster 3 (n = 291, median %mt 20.13, median
`score_eTreg` +0.0790, Fisher OR 1550.3, FDR 0) is effector-high and captures 77% of the saved
pocket cells. Cluster 16 (n = 216, %mt 19.65, `score_eTreg` −0.0738, OR 55.9, FDR 2.1e-84) is
equally mito-high and effector-low. High mitochondrial content and effector identity are therefore
separable, and %mt alone does not define the pocket.

**How to read:** One row per unsupervised leiden cluster, restricted to sorted Treg.
`fisher_or`/`fisher_p`/`fisher_fdr` test one-sided over-representation of the interactively saved
pocket cells in that cluster, BH across clusters; `frac_saved_captured` is the cluster's share of
all saved Tregs. `is_mthi_cluster` is True when FDR < 0.05 and `median_pct_mt` exceeds twice the
global Treg median; `is_mthi_effector` adds `median_score_eTreg` above the global Treg median, so a
positive `median_score_eTreg` means effector-like relative to the Treg pool. Clusters holding no
saved cells carry OR 0 and p 1 by that construction — read those rows for their %mt and depth alone.
`secondary_percell` tier — descriptive, never pseudobulk evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | `thresholds.qc_pct_mt_max = 40.0; qc_n_mads_mt = null (ceiling-only)`; script rule `fisher_fdr < 0.05 AND median_pct_mt > 2x global Treg median` | `03_results/interactive/01_qc_explore.parquet + 02_analysis/notebooks/01_qc_explore/eda/selection_sf_ctla4_pocket.csv` |

## tables/mthi_population_rule.csv

A purely quantitative rule reproduces the hand-drawn pocket. `%mt ≥ P97.5 (10.03)` and `score_eTreg
≥ the within-Treg median (−0.0031)` selects 340 Tregs, overlapping the cluster definition at Jaccard
0.626 and recovering 0.646 of the 305 saved Tregs, against 291 cells and 0.767 recovery for the
cluster rule itself. The pocket is definable without the lasso.

**How to read:** One row per candidate definition of the mt-hi effector pocket. `n_cells` is how many
sorted Tregs that rule selects. `jaccard_vs_A_effector` is intersection over union against
`A_cluster_effector`, blank on that row itself, where 1.0 would be identity. `frac_saved_recovered` is the
fraction of the 305 Tregs inside the interactively saved 335-cell lasso that the rule recaptures.
`mt_threshold` and `global_etreg_median` are the cut points behind `B_threshold`, both computed within
Treg. Agreement across rules is what this table establishes: `A_cluster_all_mthi` recovers every saved
cell (1.0) at Jaccard 0.574, the lower agreement following from its folding in the effector-low cluster.
`secondary_percell` tier — a definition audit.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | script constant `MT_PCTILE = 97.5`; cluster rule inherited from `mthi_cluster_enrichment.csv` | `03_results/interactive/01_qc_explore.parquet + 02_analysis/notebooks/01_qc_explore/eda/selection_sf_ctla4_pocket.csv` |

## tables/mthi_heat_percell.csv

The pocket does not carry the mouse 39 °C-derived program. One-sided `WT_heat_up` does sit above
normal Treg (median −0.036 against −0.060; rbc +0.48, p 4.2e-45), while the balanced
`WT_heat_updown` sits slightly lower (−0.075 against −0.060; rbc −0.17, p 6.0e-7). Both mt-hi Treg
groups fall below Tcon (−0.014) and CD8 (−0.037) on that balanced channel. The up-only shift tracks
the effector/depth axis.

**How to read:** One row per group: the two mt-hi Treg subsets, normal Treg, and the Tcon and CD8 sorts
for scale. `WT_heat_up_*` is the one-sided up-gene score and `WT_heat_updown_*` the balanced up-minus-down
score, each as a median with q25 and q75. `*_rbc_vs_normal` is the rank-biserial correlation of that group
against normal Treg (Mann-Whitney; positive means higher), populated only on the `mt_hi_effector` row. Read
the balanced channel: with tens of thousands of cells the p-values are trivially small, so sign and rbc
magnitude carry the meaning. The `caveat` column travels with the table on purpose. `secondary_percell`
tier — never pooled with the donor-pseudobulk NES.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | script constant `TIER = secondary_percell`; groups from `is_mthi_effector`/`is_mthi_cluster` in `mthi_cluster_enrichment.csv` | `03_results/interactive/01_qc_explore.parquet` |
