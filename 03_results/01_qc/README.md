# 01_qc — Quality Control

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

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
leaving 6 of 7 donors with paired SF+PB Tregs for the forest.

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
as the pocket carrying the mouse 39C signal — the pseudobulk primary-result
NES (Treg 2.53 / Tcon 2.59 / CD8 2.07, pan-T) was unchanged when these
high-mito cells were recovered.

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
