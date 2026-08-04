# 01_qc — Quality Control

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The cohort contains 7 JIA donors, of whom 6 span both arms in each analyzed population after QC. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

Figure legend sheet for `01_qc/`. MAD-based adaptive QC per GSM, doublet flagging, and the
first unsupervised embedding (usability review only) — feeds the QC usability review.

---

## figures/_overview/qc_violins_per_gsm.png

Per-GSM UMI/gene depth is adequate across strata; %mito stays low, so
no GSM is grossly degraded.

**How to read:** One violin per GSM (x), colored by sorted population; rows = UMIs,
genes (log y), %mito. Table lists per-GSM medians. QC diagnostic — no
biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.qc_n_mads = 5; qc_pct_mt_max = 40.0` | `03_results/01_qc/tables/qc_metrics_per_cell.csv` |

## figures/_overview/cells_kept_dropped.png

QC retains ~86% of cells overall, but the SF-Treg p5 library
(GSM4859852) is near-empty (median ~14 UMIs) and drops entirely,
leaving 6 of 7 donors with paired SF+PB Tregs for the donor-level
contrast.

**How to read:** Stacked bars per stratum: blue = kept, red = dropped (MAD outlier /
low-gene / doublet). Confirm every SF+PB Treg stratum retains enough
cells for pseudobulk. QC diagnostic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.qc_min_genes = 200; scrublet_expected_doublet_rate = 0.06` | `03_results/01_qc/tables/cells_kept_dropped.csv` |

## figures/_overview/unsupervised_umap.png

Sorted Treg/Tcon/CD8 occupy largely distinct transcriptomic territory;
FOXP3/IL2RA/CTLA4/IKZF2 concentrate in the Treg gate, supporting sort
fidelity.

**How to read:** Top row: cells on the unsupervised UMAP colored by sort population,
tissue, donor, leiden. Bottom row: Treg-marker expression (magma). The
UMAP is a usability lens — biology is NOT read off it. Table = leiden
x population cross-tab (contamination check).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter_viz.py` | `main` | `thresholds.hvg_n_top = 2000; n_pcs = 30; leiden_resolution = 1.0` | `03_results/objects/01_qc.h5ad (X_umap_unsupervised)` |

## figures/_overview/mthi_cluster_mt_etreg.png

Two Treg leiden clusters carry high %mt (~20% vs ~4% rest): cl6 is the
effector-like pocket (eTreg-high, SF-restricted) and cl16 is mt-hi but
eTreg-low. The pocket is a discrete, reproducibly-defined region, not
scattered noise.

**How to read:** Left: one dot per Treg leiden cluster (size ~ n cells), x = median
%mt, y = median score_eTreg; orange = mt-hi effector pocket (cl6),
purple = mt-hi non-effector (cl16). Right: per-Treg %mt vs eTreg
hexbin (log density), pocket cells overlaid green. secondary_percell /
EDA descriptive — not pseudobulk evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.qc_pct_mt_max = 40.0; qc_n_mads_mt = null (ceiling-only)` | `03_results/01_qc/tables/mthi_cluster_enrichment.csv + mthi_treg_membership.csv` |

## figures/_overview/mthi_identity_retention.png

The mt-hi effector pocket (cl6) retains Treg identity: IKZF2 (rbc
+0.59) and CTLA4 (+0.38) are UP vs normal Treg, IL2RA/TIGIT
comparable; FOXP3 is modestly lower (frac 0.55 vs 0.80) consistent
with the lower sequencing depth of high-mito cells, NOT loss of Treg
identity. These are still Tregs.

**How to read:** Grouped bars over 5 canonical Treg markers: mt-hi effector (orange),
mt-hi non-effector (purple), normal Treg (grey). Left = median log-
norm expression (rbc = rank-biserial vs normal Treg, Mann-Whitney BH-
FDR; * = FDR<0.05); right = fraction expressing. secondary_percell
tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `decisions.treg_gate.basis = sorting (refined by FOXP3/IL2RA/CTLA4/IKZF2)` | `03_results/01_qc/tables/mthi_identity_retention.csv` |

## figures/_overview/mthi_qc_discrimination.png

The pocket is not debris/dying: median n_genes 1339 (>>200 QC floor),
i.e. real cells with lower depth (the expected corollary of high mito
fraction). score_apoptosis and score_HSP are LOWER, not higher, than
normal Treg (rbc -0.09, -0.22). Doublets: no cell in the pocket is
flagged predicted_doublet, but doublet_score was not populated this
run (Scrublet gap — see reasoning note).

**How to read:** Four bars per metric across groups (eff = mt-hi effector cl6, non =
cl16, norm = normal Treg). rbc = rank-biserial vs normal Treg (Mann-
Whitney, all BH-FDR<0.05). Red dashed = 200-gene QC floor.
n_genes/total_counts are objective QC; score_apoptosis/HSP are Tier-3
hand markers (QC-descriptive, NOT evidence). secondary_percell tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.qc_min_genes = 200; qc_pct_mt_max = 40.0` | `03_results/01_qc/tables/mthi_qc_discrimination.csv` |

## figures/_overview/mthi_heat_honesty.png

WT_heat is quiet in the pocket. The balanced WT_heat_updown is
essentially flat (mt-hi effector median -0.082 vs normal -0.092);
WT_heat_up shifts up but is confounded (co-varies with the
effector/depth axis). This is secondary_percell and must NOT be read
as the pocket carrying the mouse 39 °C-derived signature. The
pseudobulk NES (Treg 2.53 / Tcon 2.59 / CD8 2.07, pan-T) was unchanged
when these high-mito cells were recovered.

**How to read:** Point = group median, bar = IQR, for WT_heat_up (left) and balanced
WT_heat_updown (right), across the two mt-hi Treg groups, normal Treg,
Tcon, CD8. secondary_percell / EDA tier — descriptive only, never
pooled with pseudobulk NES.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `decisions.go_no_go.primary_signature = WT_heat (pseudobulk, NOT per-cell)` | `03_results/01_qc/tables/mthi_heat_percell.csv` |

## figures/_overview/mthi_donor_tissue.png

The pocket is 97% synovial fluid and 69% one donor (p6); it clears
min_cells (20) in 3 SF donors but has an essentially empty PB arm
(p6=5 cells, 0 donors at floor), so it cannot support its own paired
SF-vs-PB pseudobulk contrast. It is therefore retained WITHIN the main
SF-vs-PB Treg pseudobulk and shown on the embedding, never carved out
as its own DE stratum.

**How to read:** Left: stacked bars of pocket cells per donor (SF orange / PB blue);
black dashed = per-stratum min_cells. Right: all Tregs on the
unsupervised UMAP with the mt-hi effector (orange) and non-effector
(purple) clusters highlighted. secondary_percell tier; the UMAP is a
usability lens, not biology.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize_viz.py` | `main` | `thresholds.pseudobulk_min_cells = 20; pseudobulk_min_donors = 3` | `03_results/01_qc/tables/mthi_donor_tissue.csv + mthi_treg_membership.csv` |

## tables/qc_thresholds_per_gsm.csv

The mitochondrial gate is the 40% ceiling alone in all 40 GSMs
(`mt_policy = ceiling_only`), so stressed high-mito cells survive QC by
design, while the depth windows are genuinely per-library (lower bound
189-2639 UMIs, upper bound 11k-105k). Only GSM4859852 collapses — its
MAD window spans 7-27 UMIs and under 30 genes per cell — which is why it
is the single hard-excluded library. Across the cohort the MAD rule
flags 4,992/108,414 cells (4.6%), 0.1-14.9% within a GSM.

**How to read:** One row per GSM. `log1p_total_counts_{lo,hi}` and
`log1p_n_genes_by_counts_{lo,hi}` are that GSM's MAD window (median +/-
5*MAD) on the log1p scale — exponentiate to read them as UMIs/genes; a
cell outside either window is a `mad_outlier`. `pct_mt_hi` is the
effective %mt cutoff and `pct_mt_ceiling` the hard ceiling; they are
equal here because no per-GSM MAD is applied to %mt. `excluded_gsm`
marks the hard-dropped library and `n_flagged_outlier` counts the cells
the MAD rule caught in that GSM, before the gene floor, doublet, and
exclusion filters union in. QC diagnostic — no biological claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_filter.py` | `main` | `thresholds.qc_n_mads = 5; qc_n_mads_mt = null (ceiling-only); qc_pct_mt_max = 40.0; decisions.qc.drop_gsms = [GSM4859852]` | `03_results/objects/00_build.h5ad` |

## tables/mthi_cluster_enrichment.csv

Only 2 of the 18 Treg leiden clusters carry high mito (~20% median vs
~4-5% everywhere else), and the two split on effector identity: cluster
3 (n=291, median %mt 20.13, median score_eTreg +0.0790, OR 1550.3,
FDR 0) is effector-high and captures 77% of the saved pocket cells,
whereas cluster 16 (n=216, %mt 19.65, score_eTreg -0.0738, OR 55.9,
FDR 2.1e-84) is just as mito-high but effector-LOW. High mito and
effector identity are separable, so "high %mt" alone does not define
the pocket.

**How to read:** One row per unsupervised leiden cluster, restricted to
sorted Treg. `fisher_or`/`fisher_p`/`fisher_fdr` test one-sided
over-representation of the interactively saved pocket cells in that
cluster (BH across clusters); `frac_saved_captured` is the cluster's
share of all saved Tregs. `is_mthi_cluster` is True when FDR < 0.05 AND
`median_pct_mt` exceeds twice the global Treg median;
`is_mthi_effector` adds `median_score_eTreg` above the global Treg
median, so a positive `median_score_eTreg` means effector-like relative
to the Treg pool. Clusters holding no saved cells get OR 0 / p 1 by
construction — read those rows for their %mt and depth only.
secondary_percell tier: descriptive, never pseudobulk evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | `thresholds.qc_pct_mt_max = 40.0; qc_n_mads_mt = null (ceiling-only)`; script rule `fisher_fdr < 0.05 AND median_pct_mt > 2x global Treg median` | `03_results/interactive/01_qc_explore.parquet + 02_analysis/notebooks/01_qc_explore/eda/selection_sf_ctla4_pocket.csv` |

## tables/mthi_population_rule.csv

A purely quantitative rule reproduces the hand-drawn pocket: %mt >=
P97.5 (10.03) AND score_eTreg >= the within-Treg median (-0.0031)
selects 340 Tregs that overlap the cluster definition at Jaccard 0.626
and recover 0.646 of the 305 saved Tregs, against 291 cells and 0.767
recovery for the cluster rule itself. The pocket is definable without
the lasso.

**How to read:** One row per candidate definition of the mt-hi effector
pocket. `n_cells` is how many sorted Tregs that rule selects;
`jaccard_vs_A_effector` is intersection-over-union against
`A_cluster_effector` (blank on that row itself; 1.0 would be identity);
`frac_saved_recovered` is the fraction of the 305 Tregs inside the
interactively saved 335-cell lasso that the rule recaptures.
`mt_threshold` and `global_etreg_median` are the numeric cut points
behind `B_threshold`, both computed within Treg. Agreement across rules
is the point, not any one rule: `A_cluster_all_mthi` recovers
everything (1.0) but only at Jaccard 0.574, because it folds in the
effector-LOW cluster. secondary_percell tier — a definition audit, not
evidence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | script constant `MT_PCTILE = 97.5`; cluster rule inherited from `mthi_cluster_enrichment.csv` | `03_results/interactive/01_qc_explore.parquet + 02_analysis/notebooks/01_qc_explore/eda/selection_sf_ctla4_pocket.csv` |

## tables/mthi_heat_percell.csv

The pocket does not carry the mouse 39C program. One-sided WT_heat_up
does sit above normal Treg (median -0.036 vs -0.060; rank-biserial
+0.48, p 4.2e-45), but the balanced WT_heat_updown is if anything LOWER
(-0.075 vs -0.060; rbc -0.17, p 6.0e-7), and both mt-hi Treg groups sit
below Tcon (-0.014) and CD8 (-0.037) on that balanced channel. The
up-only shift tracks the effector/depth axis, not heat.

**How to read:** One row per group: the two mt-hi Treg subsets, normal
Treg, and the Tcon and CD8 sorts for scale. `WT_heat_up_*` is the
one-sided up-gene score, `WT_heat_updown_*` the balanced up-minus-down
score, each given as median with q25/q75. `*_rbc_vs_normal` is the
rank-biserial correlation of that group against normal Treg
(Mann-Whitney; positive = higher than normal Treg) and is populated
only on the `mt_hi_effector` row. Read the balanced channel, not the
up-only one: with tens of thousands of cells the p-values are trivially
small, so sign and rbc magnitude are what carry meaning. The `caveat`
column travels with the table on purpose. secondary_percell tier —
descriptive, never pooled with the donor-pseudobulk NES.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_qc_mthi_characterize.py` | `main` | script constant `TIER = secondary_percell`; groups from `is_mthi_effector`/`is_mthi_cluster` in `mthi_cluster_enrichment.csv` | `03_results/interactive/01_qc_explore.parquet` |
