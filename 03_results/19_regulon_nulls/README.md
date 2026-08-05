# 19_regulon_nulls — artifact captions

_**Abbreviations:** TF = transcription factor, ULM = univariate linear model, HGNC = HUGO Gene
Nomenclature Committee symbol, SF = synovial fluid, PB = peripheral blood._

The neighbouring stage established that the CollecTRI HIF1A regulon's ULM activity on the sorted-Treg
SF-versus-PB donor-pseudobulk contrast holds its rank across four network variants and three
estimators, that activity scales with regulon size, and that the *direction* of the score comes from
targets many regulons share. It named two gaps in its own nulls. This stage closes one each.

**Gap one: none of the three nulls there holds target promiscuity fixed.** The random-regulon null
matches size, repressing-edge fraction and expression decile, then draws from the whole universe, so it
lands near zero and every real regulon beats it. The size-conditional residual compares against real
regulons, conditioning on size alone.

**Gap two: all three permute the annotation** — which genes a regulon claims, or which gene carries which
statistic — **and none permutes the design.** Gene-label permutation additionally destroys the gene-to-gene
correlation that makes a real contrast's statistics dependent, which is what makes it anti-conservative.

Everything here is annotation tier. No row reaches `effect_sizes_treg_arthritis.csv` or any
`03_results/master/` accumulator. Both nulls are statements about the CollecTRI HIF1A regulon's ULM
activity on this contrast. HIF1A protein activity is untested.

## Two gates run before either null

Neither is a finding and either one stops the run.

The closed-form ULM that makes 1,064 null fits affordable reproduces decoupleR's `run_ulm` on the
observed network at Spearman 1.000000, largest absolute difference 2.1e-14, so the nulls score the
statistic the headline was read off. The sign-flip's identity configuration — the refit that flips no
donor — reproduces the committed ranked list at Spearman 1.000000 over all 13,999 genes, so the null is
measured against the published contrast. The first figure is committed in `ulm_engine_validation.csv`;
the second is a run-time gate message under `regulon_nulls.signflip_identity_min_spearman = 0.99`.

## What the two nulls returned

**Holding promiscuity fixed moves the null's centre.** Rewiring the CollecTRI graph
by curveball trades preserves every regulon's size *and* every target's in-degree, so a drawn regulon
oversamples the jointly-owned high-|t| genes exactly as heavily as the observed one does. HIF1A's null
mean is therefore **+5.54** against an observed +8.30. HIF1A clears that null at empirical p 0.005, z
2.58, clearing the 95th percentile of +7.29 by about one unit, where the same score cleared the
size-and-expression-matched random regulon's 95th percentile of +2.10 by more than six. Roughly two
thirds of the score is what any regulon of that size and promiscuity profile earns on this contrast.

`null_ladder.csv` is the compact form: as each null holds more of the network's real structure fixed,
HIF1A's null centre climbs from −0.03 through +0.52 to +5.54 while the observed +8.30 never moves.

**The rewiring null does not single HIF1A out.** All eight focus factors, ordered by rewiring p:

| Factor | Targets | Observed | Null mean | p | z |
|---|---|---|---|---|---|
| REL | 83 | 7.56 | 3.40 | 0.001 | 3.26 |
| HIF1A | 293 | 8.30 | 5.54 | 0.005 | 2.58 |
| NFKB1 | 280 | 8.29 | 5.45 | 0.006 | 2.47 |
| HSF1 | 67 | 5.60 | 2.83 | 0.011 | 2.18 |
| ATF3 | 72 | 3.96 | 1.68 | 0.047 | 1.71 |
| STAT3 | 285 | 6.89 | 4.98 | 0.055 | 1.63 |
| CREB1 | 276 | 6.87 | 5.43 | 0.106 | 1.23 |
| EPAS1 | 62 | 2.79 | 2.79 | 0.48 | 0.003 |

STAT3 and CREB1 fail this null, and NFKB1 remains indistinguishable from HIF1A. The smaller regulons
doing relatively better than the large ones is the ordering the size-conditional residual produced by a
different route. EPAS1 is the clean negative and the useful one: observed +2.794 against a null mean of
+2.790. HIF2A's regulon here is exactly what its size and promiscuity predict, which is what a factor
with no signal looks like.

**Permuting the design says the contrast is real and says nothing about which factor.** Swapping the SF
and PB labels within donor and refitting end to end is the exchangeability this paired design licenses.
Six Treg donors carry both arms, so all 2⁶ = 64 configurations are enumerated and the test is exact —
no seed enters it. The observed labelling gives HIF1A the largest score of all 64 (+8.30 against a null
maximum of +6.63), which is p = 1/64 = 0.0156, the finest this many donors can resolve. Seven of the
eight focus factors reach that same floor. ATF3 is the one that does not, with seven permuted
configurations beating it (p 0.125).

So the sign-flip answers *is there signal in this contrast at all* — yes, and essentially every regulon's
activity is maximal at the true labels — with no power to separate one factor from another. The rewiring
null is the one with resolving power, and it resolves that HIF1A's margin over a promiscuity-matched
regulon is real but small, shared with its neighbours, and absent for EPAS1.

## Two ways to misread this stage

**A floor is a weak p-value.** `p_exact_one_sided = 0.0156` for seven factors at once is one bit of
information about the contrast. A smaller p-value here takes more paired donors; more computation
cannot reach it.

**The rewiring null is the strongest available and remains short of a claim.** It holds the network's
degree structure fixed, the sharpest constraint on offer, and stays a question about one curated
network's target assignments. This stage is silent on hypoxia, on temperature, and on HIF1A protein.

## tables/source_hash_manifest.csv

The CollecTRI human regulon table read across the compartment
boundary is pinned at sha256 4473c918..., so a change to that network
halts this stage. Without the pin, both nulls would move in silence.

**How to read:** One row per cross-compartment source: `source_label` is the name this
stage refers to it by, `source_path` is repository-root-relative,
`sha256` the digest of the bytes actually read. The first run writes
the pin; every later run verifies against it and stops on a mismatch.
This stage keeps its own pin, so a moved network stops both stages.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `verify_source_hash` | `unbiased_enrichment.tf_network.path = ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` | `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/ulm_engine_validation.csv

The closed-form univariate-linear-model score used for all 1,064 null
draws reproduces decoupleR's `run_ulm` on the observed network
exactly: Spearman 1.000000 and a largest absolute difference of
2.13e-14 over 601 factors, which is machine precision, so the nulls
measure the same statistic the headline was read off.

**How to read:** One row. `score_closed_form` against `score_decoupler` is a
like-for-like comparison on the observed signed network, summarised
by `spearman`, `pearson` and `max_abs_diff`. decoupleR regresses
every gene's contrast statistic on the regulon's mode of regulation,
zero for a non-target, so the fit has a closed form costing one
sparse matrix-vector product for all factors at once, which is what
makes tens of thousands of null fits feasible. `gate` reads `pass`
only when `spearman` reaches `gate_min_spearman`; the run stops
otherwise. This is a guard, and it satisfies the compartment's
standing requirement that an engine swap be shown by rank correlation
to be a method change.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `ulm_fast` | `regulon_nulls.ulm_validation_min_spearman = 0.9999; thresholds.gsea_min_size = 5` | `03_results/03_pseudobulk/tables/ranked_treg.tsv, ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/rewiring_null.csv

Held against random regulons that share its size and every one of its
targets' in-degrees, HIF1A's CollecTRI-ULM score of 8.30 clears the
null at empirical p 0.005, z 2.58, and the null's centre sits at
5.54, so most of the score is what any regulon of that size and
promiscuity profile earns on this contrast. The null does not single
HIF1A out: REL clears it harder (p 0.001, z 3.26) on 83 targets,
NFKB1 is indistinguishable (p 0.006), and EPAS1 sits exactly on its
own null (obs 2.79 against a null mean of 2.79, p 0.48).

**How to read:** One row per focus factor. `obs_score` and `obs_rank` are the observed
signed-network ULM values; the `null_*` columns describe the draws. A
curveball trade holds the targets two regulons share and
redistributes the rest between them, so every regulon's size and
every target's in-degree stay invariant and only the assignment of
targets to factors is randomised — which is why `null_mean` sits far
above zero and why beating this null means more than beating a random
gene set of matched size. `mean_target_indeg` is the observed
regulon's average target promiscuity, the property this null holds
fixed. `p_empirical` is (draws at or above obs + 1)/(draws + 1), so
it never returns exactly zero, and `z_vs_null` standardises the gap.
Because size is preserved exactly, `null_rank_median` and
`null_rank_best` share the `n_tfs_scored` denominator with
`obs_rank`. Annotation tier: this is a statement about the regulon's
target set, and protein activity is untested.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `curveball_trade + run_trades + ulm_fast` | `regulon_nulls.rewiring_draws = 1000; regulon_nulls.rewiring_burn_in_trades_per_edge = 5; regulon_nulls.rewiring_trades_between_draws_per_tf = 2; tf_activity.focus_tfs = 8; thresholds.gsea_min_size = 5` | `03_results/03_pseudobulk/tables/ranked_treg.tsv, ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/rewiring_null_draws.csv

The per-draw substrate behind the rewiring null: 1000
degree-preserving rewirings scored for each of the 8 focus factors
present in this configuration.

**How to read:** One row per (draw, factor). `score` is that draw's closed-form ULM
value and `rank` its position by descending score among every factor
scored in the same draw. Draws come from one Markov chain advanced
between samples by trades, so consecutive rows are decorrelated and
dependent by construction. This is the distribution
`rewiring_null.csv` summarises; read it when a summary statistic
needs checking against the shape it came from.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `run_trades + ulm_fast` | `regulon_nulls.rewiring_draws = 1000` | `in-memory curveball chain` |

## tables/signflip_null.csv

Permuting the experimental design with the annotation held fixed, the
observed synovial-fluid-versus-blood labelling gives HIF1A the
largest score of all 64 within-donor configurations (8.30 against a
null maximum of 6.63), which is exact p = 0.0156, the finest this
design can resolve. It carries no information about HIF1A
specifically: 7 of the 8 focus factors reach the same floor, ATF3
being the one that does not (p 0.125). What it establishes is that
the contrast itself survives relabelling.

**How to read:** One row per focus factor. The two tissue labels are swapped within a
donor and the whole limma-voom contrast is refitted, which is the
exchangeability a paired design licenses and which leaves the
correlation between genes intact, where permuting gene labels would
destroy it. `n_paired_donors` donors carry both arms, so
`n_configurations` = 2^that, and every one is enumerated: the test is
exact and no seed enters it. `n_ge_obs` counts configurations scoring
at or above the observed, the observed one included, and
`p_exact_one_sided` is that count over `n_configurations`. `p_floor`
is 1/`n_configurations`, and `at_resolution_floor` is TRUE when only
the observed configuration reaches the observed score, i.e. the
p-value is as small as this many donors can make it. Flipping every
donor negates the contrast, so the configurations come in
sign-symmetric pairs and `null_mean` sits near zero. The gene set is
filtered once on the observed design and held fixed, so scores are
comparable across configurations; voom is recomputed per
configuration, because its weights depend on the design.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `flip_design + fit_contrast + ulm_fast` | `design.donor_key = donor; design.tissue_key = tissue; regulon_nulls.signflip_max_exact_donors = 20; regulon_nulls.signflip_identity_min_spearman = 0.99` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv, 03_results/03_pseudobulk/tables/pseudobulk_coldata.csv, 03_results/03_pseudobulk/tables/gene_symbols.csv, ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/signflip_null_draws.csv

The per-configuration substrate behind the sign-flip null: all 64
within-donor label swaps of the 6 paired donors, scored for each
focus factor.

**How to read:** One row per (configuration, factor). `donors_flipped` names the
donors whose two tissue labels were swapped, pipe-delimited, and
`is_observed` marks the configuration that flips none — the published
contrast. `score` and `rank` are that configuration's closed-form ULM
value and its position among every factor scored. Configurations are
enumerated in a fixed order, so this table is byte-stable across
runs. Read it for the null's shape, and for the observed row being
the extreme one for most factors.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `flip_design + fit_contrast + ulm_fast` | `design.donor_key = donor; design.tissue_key = tissue` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv + pseudobulk_coldata.csv` |

## tables/null_ladder.csv

Set side by side, the 4 nulls order by what they hold fixed, and the
conclusion weakens along that order: HIF1A clears a
size-and-expression-matched random regulon by a wide margin, clears a
size-and-promiscuity-matched real-degree rewiring by a narrow one,
and is one of 7 focus factors that max out the exact design
permutation. Each rung is published beside the others because no
single rung supports a claim alone.

**How to read:** One row per (factor, null), ordered weakest-to-strongest by what the
null holds constant. `holds_fixed` names that and `permutes` names
what is randomised, which together bound what a p-value from that row
can mean. `obs` is the same observed CollecTRI-ULM score in every row
for a given factor; only the reference distribution changes, so the
columns are comparable down a factor's block and the movement in
`p_empirical` is the whole content of the table. `null_q95` is empty
for the sign-flip rung, whose 64 configurations do not support a
stable upper percentile. The two random-regulon rungs are re-read
from their own published table.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `bind_rows` | `tf_activity.focus_tfs = 8; tf_activity.null_draws; regulon_nulls.rewiring_draws` | `03_results/18_tf_activity/tables/size_matched_null.csv, 03_results/19_regulon_nulls/tables/rewiring_null.csv + signflip_null.csv` |

