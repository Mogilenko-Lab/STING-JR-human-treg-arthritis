# 19_regulon_nulls — artifact captions

_**Abbreviations:** TF = transcription factor, ULM = univariate linear model, HGNC = HUGO Gene
Nomenclature Committee symbol, SF = synovial fluid, PB = peripheral blood._

The neighbouring stage established that the CollecTRI HIF1A regulon's ULM activity on the
sorted-Treg SF-versus-PB donor-pseudobulk contrast holds its rank across four network variants
and three estimators, that activity scales with regulon size, and that the *direction* of the
score is carried by targets many regulons share. It also named two gaps in its own nulls, and
this stage closes one each.

The first gap: none of the three nulls there holds **target promiscuity** fixed. The
random-regulon null matches size, repressing-edge fraction and expression decile but draws from
the whole universe, so it lands near zero and every real regulon beats it. The size-conditional
residual compares against real regulons but conditions on size alone. The second gap: all three
permute the **annotation** — which genes a regulon claims, or which gene carries which statistic
— and none permutes the **design**. Gene-label permutation additionally destroys the
gene-to-gene correlation that makes a real contrast's statistics dependent, which is what makes
it anti-conservative.

Everything here is annotation tier. No row reaches `effect_sizes_treg_arthritis.csv` or any
`03_results/master/` accumulator, and both nulls are statements about the CollecTRI HIF1A
regulon's ULM activity on this contrast — never about HIF1A protein.

## Two gates run before either null

Neither is a finding and either stops the run. The closed-form ULM that makes 1,064 null fits
affordable reproduces decoupleR's `run_ulm` on the observed network at Spearman 1.000000 with a
largest absolute difference of 2.1e-14, so the nulls score the same statistic the headline was
read off. And the sign-flip's identity configuration — the refit that flips no donor —
reproduces the committed ranked list at Spearman 1.000000 over all 13,999 genes, so the null is
measured against the published contrast rather than a near-miss of it.

## What the two nulls returned

**Holding promiscuity fixed moves the null's centre, not just its tail.** Rewiring the CollecTRI
graph by curveball trades preserves every regulon's size *and* every target's in-degree, so a
drawn regulon oversamples the jointly-owned high-|t| genes exactly as heavily as the observed one
does. The consequence is that the null mean for HIF1A is **+5.54**, not zero, against an observed
+8.30. HIF1A still clears it — empirical p 0.005, z 2.58 — but it clears the 95th percentile of
+7.29 by about one unit, where the same score cleared the size-and-expression-matched random
regulon's 95th percentile of +2.10 by more than six. Roughly two thirds of the score is what any
regulon of that size and promiscuity profile earns on this contrast.

`null_ladder.csv` is the compact form of that: as each null holds more of the network's real
structure fixed, its centre climbs from -0.03 through +0.52 to +5.54 while the observed +8.30
never moves.

**And it does not single HIF1A out.** REL clears the same null harder on 83 targets (p 0.001, z
3.26), HSF1 clears it at p 0.011 on 67, NFKB1 remains indistinguishable at p 0.006 against
HIF1A's 0.005, and ATF3 reaches p 0.047 on 72. The smaller regulons doing relatively better than
the large ones is the same ordering the size-conditional residual produced by a different route.
EPAS1 is the clean negative and the useful one: observed +2.794 against a null mean of +2.790,
p 0.48, z 0.003. HIF2A's regulon on this contrast is *exactly* what its size and promiscuity
predict, which is what a factor with no signal here should look like.

**Permuting the design says the contrast is real and says nothing about which factor.** Swapping
the SF and PB labels within donor and refitting end to end is the exchangeability this paired
design licenses. Six Treg donors carry both arms, so all 2^6 = 64 configurations are enumerated
and the test is exact — no seed enters it. The observed labelling gives HIF1A the largest score
of all 64 (+8.30 against a null maximum of +6.63), which is p = 1/64 = 0.0156, the finest this
many donors can resolve. Seven of the eight focus factors reach that same floor; ATF3 is the only
one that does not, with seven permuted configurations beating it (p 0.125).

So the sign-flip answers *is there signal in this contrast at all* — yes, and essentially every
regulon's activity is maximal at the true labels — while having no power to separate one factor
from another. The rewiring null is the one with resolving power, and what it resolves is that
HIF1A's margin over a promiscuity-matched regulon is real but small, shared with its neighbours,
and absent for EPAS1.

## Two ways to misread this stage

**A floor is not a strong p-value.** `p_exact_one_sided = 0.0156` for seven factors at once is
one bit of information about the contrast, not seven findings about seven factors. A smaller
p-value here needs more paired donors, not more computation.

**The rewiring null is the strongest available and still not a claim.** It holds the network's
degree structure fixed, which is the sharpest constraint on offer, but it is still a question
about one curated network's target assignments. Nothing in this stage is evidence about hypoxia,
temperature, or HIF1A protein, and the reading stays correlative.

## tables/source_hash_manifest.csv

The CollecTRI human regulon table read across the compartment
boundary is pinned at sha256 4473c918..., so a change to that network
stops this stage instead of quietly moving both nulls.

**How to read:** One row per cross-compartment source: `source_label` is the name this
stage refers to it by, `source_path` is repository-root-relative, and
`sha256` is the digest of the bytes actually read. The first run
writes the pin; every later run verifies against it and stops on a
mismatch. This stage keeps its own pin rather than trusting the
sibling stage's, so a moved network stops both. Verification is the
only gate.

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
zero for a non-target, so the fit has a closed form that costs one
sparse matrix-vector product for all factors at once — which is what
makes tens of thousands of null fits feasible. `gate` is `pass` only
when `spearman` reaches `gate_min_spearman`; the run stops otherwise.
This is a guard, not a measurement, and it is the compartment's
standing requirement that an engine swap be shown by rank correlation
to be a method change rather than a result change.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `ulm_fast` | `regulon_nulls.ulm_validation_min_spearman = 0.9999; thresholds.gsea_min_size = 5` | `03_results/03_pseudobulk/tables/ranked_treg.tsv, ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/rewiring_null.csv

Held against random regulons that share its size AND every one of its
targets' in-degrees, HIF1A's CollecTRI-ULM score of 8.30 still clears
the null (empirical p 0.005, z 2.58) but the null's centre is 5.54
rather than zero, so most of the score is what any regulon of that
size and promiscuity profile earns on this contrast. The null does
not single HIF1A out: REL clears it harder (p 0.001, z 3.26) on 83
targets, NFKB1 is indistinguishable (p 0.006), and EPAS1 sits exactly
on its own null (obs 2.79 against a null mean of 2.79, p 0.48).

**How to read:** One row per focus factor. `obs_score` and `obs_rank` are the observed
signed-network ULM values; the `null_*` columns describe the draws. A
curveball trade holds the targets two regulons share and
redistributes the rest between them, so every regulon's size and
every target's in-degree are invariant and only the assignment of
targets to factors is randomised — which is why `null_mean` is far
above zero and why beating this null means more than beating a random
gene set of matched size. `mean_target_indeg` is the observed
regulon's average target promiscuity, the property this null holds
fixed. `p_empirical` is (draws at or above obs + 1)/(draws + 1) so it
never returns exactly zero, and `z_vs_null` standardises the gap.
Because size is preserved exactly, the same factors clear minsize in
every draw, so `null_rank_median` and `null_rank_best` share the
`n_tfs_scored` denominator with `obs_rank`. Annotation tier: this is
a statement about the regulon's target set, never about the
transcription factor's protein activity.

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
between samples rather than restarted, so consecutive rows are
decorrelated by trades and not independent by construction. This is
the distribution `rewiring_null.csv` summarises; read it when a
summary statistic needs checking against the shape it came from.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `run_trades + ulm_fast` | `regulon_nulls.rewiring_draws = 1000` | `in-memory curveball chain` |

## tables/signflip_null.csv

Permuting the experimental design rather than the annotation, the
observed synovial-fluid-versus-blood labelling gives HIF1A the
largest score of all 64 within-donor configurations (8.30 against a
null maximum of 6.63), which is exact p = 0.0156, the finest this
design can resolve. It carries no information about HIF1A
specifically: 7 of the 8 focus factors reach the same floor, ATF3
being the only one that does not (p 0.125). What it establishes is
that the contrast itself is not an artefact of labelling.

**How to read:** One row per focus factor. The two tissue labels are swapped within a
donor and the whole limma-voom contrast is refitted, which is the
exchangeability a paired design licenses and which — unlike permuting
gene labels — leaves the correlation between genes intact.
`n_paired_donors` donors carry both arms, so `n_configurations` =
2^that, and ALL of them are enumerated: the test is exact and no seed
enters it. `n_ge_obs` counts configurations scoring at or above the
observed, the observed one included, and `p_exact_one_sided` is that
count over `n_configurations`. `p_floor` is 1/`n_configurations` and
`at_resolution_floor` is TRUE when only the observed configuration
reaches the observed score, i.e. the p-value is as small as this many
donors can make it and a smaller one would need more donors, not more
computation. Flipping every donor negates the contrast, so the
configurations come in sign-symmetric pairs and `null_mean` sits near
zero by construction. The gene set is filtered once on the observed
design and held fixed, so scores are comparable across
configurations, while voom is recomputed per configuration because
its weights depend on the design.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `flip_design + fit_contrast + ulm_fast` | `design.donor_key = donor; design.tissue_key = tissue; regulon_nulls.signflip_max_exact_donors = 20; regulon_nulls.signflip_identity_min_spearman = 0.99` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv, 03_results/03_pseudobulk/tables/pseudobulk_coldata.csv, 03_results/03_pseudobulk/tables/gene_symbols.csv, ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/signflip_null_draws.csv

The per-configuration substrate behind the sign-flip null: all 64
within-donor label swaps of the 6 paired donors, scored for each
focus factor.

**How to read:** One row per (configuration, factor). `donors_flipped` names the
donors whose two tissue labels were swapped, pipe-delimited, and
`is_observed` marks the single configuration that flips none — the
published contrast. `score` and `rank` are that configuration's
closed-form ULM value and its position among every factor scored.
Configurations are enumerated in a fixed order, so this table is
byte-stable across runs. Read it to see the null's shape, and in
particular that the observed row is the extreme one for most factors
rather than merely a high one.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `flip_design + fit_contrast + ulm_fast` | `design.donor_key = donor; design.tissue_key = tissue` | `03_results/03_pseudobulk/tables/pseudobulk_counts.csv + pseudobulk_coldata.csv` |

## tables/null_ladder.csv

Set side by side, the four nulls order by what they hold fixed and
the conclusion weakens along that order: HIF1A clears a
size-and-expression-matched random regulon by a wide margin, clears a
size-and-promiscuity-matched real-degree rewiring by a narrow one,
and is one of seven focus factors that max out the exact design
permutation. No single rung supports a claim on its own, which is why
they are published as one table.

**How to read:** One row per (factor, null), ordered weakest-to-strongest by what the
null holds constant. `holds_fixed` names exactly that and `permutes`
names what is randomised, which together are what a p-value from that
row can and cannot mean. `obs` is the same observed CollecTRI-ULM
score in every row for a given factor; only the reference
distribution changes, so the columns are comparable down a factor's
block and the movement in `p_empirical` is the whole content of the
table. `null_q95` is empty for the sign-flip rung because its 64
configurations do not support a stable upper percentile. The two
random-regulon rungs are re-read from their own published table
rather than recomputed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_regulon_nulls.R` | `bind_rows` | `tf_activity.focus_tfs = 8; tf_activity.null_draws; regulon_nulls.rewiring_draws` | `03_results/18_tf_activity/tables/size_matched_null.csv, 03_results/19_regulon_nulls/tables/rewiring_null.csv + signflip_null.csv` |

