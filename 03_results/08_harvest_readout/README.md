# 08_harvest_readout — Per-cell readouts on one scale

Two MSigDB Hallmark programs, HALLMARK_HYPOXIA and HALLMARK_UNFOLDED_PROTEIN_RESPONSE, are scored
per cell on the frozen annotation with the same rank-based AUCell and UCell engine the rest of the
project uses. They sit beside the readouts this compartment already carries — `score_HSP`,
`score_eTreg`, and the mouse anchor — so seven channels read on one scale in one tidy table.

**The mouse anchor enters as three channels.** `WT_heat_up`, `WT_heat_down` and the balanced
`WT_heat_updown` composite are all carried. The composite is an up-minus-down difference, so it
cancels whenever both arms shift the same way and reports flat regardless of what either arm did.
Carrying the arms separately is what makes the anchor measurable here, and it matches how the
donor-pseudobulk enrichment treats the two sets: scored one at a time, never differenced.

This stage publishes no figure. Its three tables are a compute resource, every score is carried
in verbatim or scored once here, and nothing enters the cross-dataset effect-size accumulator.

**Hypoxia here is a transcriptional readout** consistent with a low-oxygen, metabolically
stressed state. It carries no HIF-causality claim.

## Signature provenance

The seven channels draw on four external sources.

| Channel | Origin | How the list was derived |
|---|---|---|
| `HALLMARK_HYPOXIA` (200), `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` (113) | **MSigDB Hallmark collection H**, *Homo sapiens*, **v2026.1.Hs**, retrieved offline through **msigdbr 26.1.0**. | Frozen one symbol per line under `00_data/references/msigdb_hallmark/`, each with a validated expected size, and used whole. |
| `WT_heat_up`, `WT_heat_down`, `WT_heat_updown` | **GSE329522**, this project's own mouse anchor. No paper reference recorded. | Bulk RNA-seq of induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype (WT, cGAS-KO) × temperature (37 °C, 39 °C), 5 biological replicates per group. The WT 39 °C-against-37 °C up and down arms, projected to human orthologs with pinned offline babelgene. The composite is their up-minus-down difference. |
| `score_eTreg` | **GSE161426**, 26 bulk RNA-seq samples of sorted CD4 populations. Mijnheer / Lutter et al. 2021, *Nature Communications*, PMID 33976194, doi 10.1038/s41467-021-22975-7. | Derived here as a synovial-fluid Treg against peripheral-blood Treg contrast on the deposited log2-normalised matrix `GSE161426_Gene_expression_table_log2.xlsx`. |
| `score_HSP` | This compartment's own hand marker module. | A small curated heat-shock marker panel. |

---

## Tables

### `tables/harvest_readout_summary.csv`

One row per `coarse_label` × `tissue` × `readout`, over the seven channels. AUCell scores are
rank-based in [0, 1]. The carried module scores are mean-centred, so their zero is arbitrary and
only differences between strata are interpretable.

`mean` / `sd` / `median` summarise the stratum, `n_cells` its size, and `mean_pct_counts_mt` its
mean mitochondrial fraction as a viability cross-check. `frac_above_p90` is the fraction above
the readout's whole-dataset 90th percentile (`global_p90`), so it picks out a high-expressing
pocket where `mean` reads the whole distribution. `sf_minus_pb_mean` and `sf_minus_pb_smd` give
the synovial-minus-blood shift as a raw mean difference and a standardised Cohen's d, **positive
meaning higher in synovial fluid**.

Split into its arms, the mouse anchor sits higher in synovial than in paired blood T cells. The
up arm rises by d ≈ +0.98 in Treg, its largest value across the three states, and the down arm
rises too, by d ≈ +0.68 in Treg. That is precisely why the balanced composite reads flat at
d ≈ −0.10. The hypoxia readout rises in every state and is strongest in Treg (d ≈ +1.67), with
about 30% of synovial Tregs above the whole-dataset 90th percentile against about 1.7% in blood.
The unfolded-protein readout stays broadly flat.

**These pool cells across donors**, so the unit of replication is the cell. Read them as
descriptive per-cell magnitudes and use the donor-level companion for anything compared against
donor-pseudobulk results.

### `tables/harvest_readout_donor_means.csv`

One row per `donor` × `tissue` × `coarse_label`, with one column per readout holding that donor's
mean over its cells in that stratum, and `n_cells` the stratum size. Strata below
`pseudobulk_min_cells` are dropped so this level sees the same donors the pseudobulk aggregation
does. Here nothing is dropped, and the thinnest surviving stratum still holds 266 cells.

Each of the seven donors carries every readout in both tissues except the two strata absent by
design — no synovial Treg for patient 5, no blood Tcon or CD8 for patient 3 — leaving six donors
that pair synovial against blood in each sorted state.

**This is the audit substrate for the donor contrast below.** It is the file to read when a
donor-level effect looks driven by one patient.

### `tables/harvest_readout_donor_contrast.csv`

One row per `coarse_label` × `readout`. Each donor contributes one mean per tissue from the table
above, and the contrast is taken **across donors, paired within donor** — GSE160097 is a paired
design, and `n_donors_paired` records how many donors carry both tissues, with `n_donors_sf` and
`n_donors_pb` giving the per-arm counts.

`sf_minus_pb_mean` and `sf_minus_pb_sd` are the mean and spread of the within-donor differences,
positive meaning higher in synovial fluid. `sf_minus_pb_dz` is the paired standardised effect,
with `dz_se`, `dz_ci_low` / `dz_ci_high` a 95% interval and `dz_pvalue` a paired-t p on five
degrees of freedom. `sf_minus_pb_smd_unpaired` recomputes the same donor means ignoring the
pairing, for readers comparing against unpaired donor-level effect sizes.
`min_cells_per_donor_stratum` is the thinnest donor stratum behind that cell state.

Paired within donor across n = 6, `WT_heat_up` is higher in synovial fluid in every sorted state
(Treg dz = +2.18, 95% CI 0.71–3.65, p = 0.003; Tcon +3.16; CD8 +3.26) — **and so is the down arm,
harder in Treg and CD8** (+3.50; +4.27), which leaves the directional composite flat in Treg
(dz = −0.38, p = 0.39) and negative in CD8 (−4.63).

**Read that as a property of the scoring channel.** The two arms are anti-correlated by
construction, so their moving together on the same cells reports a shared non-directional
component — general activation or expression breadth lifting large rank-based scores. The
donor-pseudobulk enrichment keeps the arms apart, and where the two disagree it has standing. The
hypoxia, effector-Treg and heat-shock readouts carry no such internal contradiction and are read
as written.

**Read this table next to donor-pseudobulk results**, because it shares their unit of
replication. It stays annotation tier — six donors, uncorrected across readouts, on carried-in
scores — and the pseudobulk enrichment stays the result that carries a claim.

---

The per-cell feed for the reactive review notebook is
[`../interactive/08_harvest_readout.parquet`](../interactive/), regenerable and untracked.
