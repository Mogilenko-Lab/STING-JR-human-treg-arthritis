# 08_harvest_readout — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint), PB = peripheral blood. The SF-vs-PB comparison is within the paired JIA donors of GSE160097. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory, Tcon = CD4⁺CD25⁻ conventional, CD8 = CD8⁺CD45RO⁺ memory._

I added two MSigDB Hallmark per-cell readouts — HALLMARK_HYPOXIA and HALLMARK_UNFOLDED_PROTEIN_RESPONSE — to the frozen JIA SF/PB annotation, scored with the same rank-based AUCell+UCell engine the SAVI STING-reference compartment uses, so the two compartments read on one scale. Hypoxia here is a **readout** consistent with a low-oxygen, metabolically stressed state, and makes no HIF-causality claim. These per-cell scores are a **secondary, annotation-tier** view, never pooled with the confirmatory donor-pseudobulk enrichment, and they carry alongside the compartment's already-derived score_HSP, score_eTreg and mouse-anchor WT_heat readouts for one tidy comparison. The per-cell feed for the reactive review lives at `03_results/interactive/08_harvest_readout.parquet` (regenerable, not committed).

The mouse 39 °C anchor now enters as **three** readouts — `WT_heat_up`, `WT_heat_down`, and the balanced `WT_heat_updown` composite — rather than the composite alone. The composite is an up-minus-down difference, so it cancels whenever both arms shift the same way and reports flat regardless of what either arm did. Carrying the arms separately is what makes the anchor measurable here, and it matches how the donor-pseudobulk enrichment treats the two sets: scored one at a time, never differenced.

Every score in this directory is carried in verbatim or scored once here and is read correlatively. Nothing in it enters the cross-dataset effect-size accumulator.

## tables/harvest_readout_summary.csv

Split into its arms, the mouse 39 °C anchor is higher in synovial than in paired blood T cells — the up arm by d ≈ +0.98 in Treg, its largest value across the three sorted states — but the down arm rises too (d ≈ +0.68 in Treg), which is precisely why the balanced composite reads flat at d ≈ −0.10; the hypoxia readout meanwhile rises in every state and is strongest in Treg (d ≈ +1.67), with ~30% of SF Tregs above the whole-dataset 90th percentile against ~1.7% in blood, while UPR stays broadly flat.

**How to read:** One row per (`coarse_label` × `tissue` × `readout`). The seven readouts are the two Hallmark AUCell scores, the heat-shock module `score_HSP`, the effector-Treg module `score_eTreg`, and the mouse anchor's `WT_heat_up`, `WT_heat_down` and `WT_heat_updown` channels. AUCell scores are rank-based in [0,1]; the carried module scores are mean-centred, so their zero is arbitrary and only differences between strata are interpretable.

`mean`/`sd`/`median` summarise the stratum, `n_cells` its size, `mean_pct_counts_mt` its mean mitochondrial fraction as a viability cross-check. `frac_above_p90` is the fraction above the readout's whole-dataset 90th percentile (`global_p90`), flagging a high-expressing pocket rather than a whole-distribution shift.

`sf_minus_pb_mean` and `sf_minus_pb_smd` give the SF-minus-PB shift for that cell-state as a raw mean difference and a standardized Cohen's d, **positive meaning higher in synovial fluid**. These pool cells across donors, so the unit of replication is the cell, not the donor — read them as descriptive per-cell magnitudes, not as tests, and use the donor-level companion below for anything compared against donor-pseudobulk results. Correlative, secondary tier; claim tier L3 per-cell statistic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/08_harvest_readout.py` | `main` → `summarise` | `percell_score_ncores=8` | `03_results/objects/02_annotation.h5ad`, `00_data/references/msigdb_hallmark/{HALLMARK_HYPOXIA,HALLMARK_UNFOLDED_PROTEIN_RESPONSE}.txt`, `03_results/interactive/{05_gonogo_explore,01_qc_explore}.parquet` |

## tables/harvest_readout_donor_means.csv

Each of the 7 JIA donors carries every readout in both tissues except two strata absent by design — no SF Treg for patient 5 and no PB Tcon/CD8 for patient 3 — leaving 6 donors that pair SF against blood in each sorted state, and the thinnest surviving stratum still holds 266 cells, so no donor is lost to the cell-count floor.

**How to read:** One row per (`donor` × `tissue` × `coarse_label`), with one column per readout holding that donor's mean over its cells in that stratum and `n_cells` the stratum size. Strata below the donor-stratum floor `pseudobulk_min_cells` are dropped so this level sees the same donors the donor-pseudobulk aggregation does; here nothing is dropped. This is the audit substrate for the donor contrast below — it is the file to read if a donor-level effect looks driven by one patient. Annotation tier, correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/08_harvest_readout.py` | `main` → `donor_means` | `pseudobulk_min_cells=20` | `03_results/objects/02_annotation.h5ad`, `00_data/references/msigdb_hallmark/{HALLMARK_HYPOXIA,HALLMARK_UNFOLDED_PROTEIN_RESPONSE}.txt`, `03_results/interactive/{05_gonogo_explore,01_qc_explore}.parquet` |

## tables/harvest_readout_donor_contrast.csv

Moving the contrast from cells to donors does not weaken the anchor's up arm — it strengthens it: paired within donor across n = 6, `WT_heat_up` is higher in SF in every sorted state (Treg dz = +2.18, 95% CI 0.71–3.65, p = 0.003; Tcon +3.16; CD8 +3.26), the down arm rises in parallel (Treg +3.50), and the composite is the one channel that finds nothing in Treg (dz = −0.38, CI −1.21–0.45, p = 0.39), confirming that the flat composite reflects cancellation between two rising arms rather than an absent per-cell signal.

**How to read:** One row per (`coarse_label` × `readout`). Each donor contributes one mean per tissue from the donor-means table, and the contrast is taken **across donors**, paired within donor — GSE160097 is a paired SF/PB design, and `n_donors_paired` records how many donors actually carry both tissues, with `n_donors_sf`/`n_donors_pb` giving the per-arm counts.

`sf_minus_pb_mean` and `sf_minus_pb_sd` are the mean and spread of the within-donor SF-minus-PB differences; **positive means higher in synovial fluid**. `sf_minus_pb_dz` is the paired standardized effect (mean difference over its own SD), with `dz_se`, `dz_ci_low`/`dz_ci_high` a 95% interval and `dz_pvalue` a paired-t p-value on 5 degrees of freedom. `sf_minus_pb_smd_unpaired` recomputes the same donor means ignoring the pairing, for readers comparing against unpaired donor-level effect sizes. `min_cells_per_donor_stratum` is the thinnest donor stratum behind that cell-state.

Read this one, not the per-cell shift, next to donor-pseudobulk results: it shares their unit of replication. It remains **annotation tier** — six donors, uncorrected across readouts, on scores that were carried in rather than re-derived, so it corroborates the pseudobulk enrichment and never substitutes for it, and it contributes no row to the cross-dataset effect-size accumulator. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/08_harvest_readout.py` | `main` → `donor_contrast` → `paired_effect`, `unpaired_smd` | `pseudobulk_min_cells=20`, `tissue_levels.synovial_fluid=synovial_fluid`, `tissue_levels.peripheral_blood=peripheral_blood` | `03_results/08_harvest_readout/tables/harvest_readout_donor_means.csv` (in-memory from `main`) |
