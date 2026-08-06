# 19_regulon_nulls — Two harder nulls for the same activity score

The neighbouring stage established that the CollecTRI HIF1A regulon's ULM activity on the
sorted-Treg synovial-versus-blood contrast holds its rank across four network variants and three
estimators, that activity scales with regulon size, and that the *direction* of the score comes
from targets many regulons share. It named two gaps in its own nulls. This stage closes one each.

**Gap one: none of the three nulls there holds target promiscuity fixed.** The random-regulon
null matches size, repressing-edge fraction and expression decile, then draws from the whole
universe, so it lands near zero and every real regulon beats it. The size-conditional residual
compares against real regulons, conditioning on size alone.

**Gap two: all three permute the annotation** — which genes a regulon claims, or which gene
carries which statistic — **and none permutes the design.** Gene-label permutation additionally
destroys the gene-to-gene correlation that makes a real contrast's statistics dependent, which is
what makes it anti-conservative.

**What was computed.** A degree-preserving rewiring null at 1,000 draws, an exact sign-flip null
enumerating all 64 within-donor label swaps, a four-rung ladder placing them beside the earlier
nulls, and two engine gates. This stage publishes no figure — its six tables are the deliverable.

**Tier.** Annotation. No row reaches [`../master/`](../master/). Both nulls are statements about
the CollecTRI HIF1A regulon's ULM activity on this contrast, and HIF1A protein activity stays
untested.

## Two gates run before either null

Neither is a finding, and either one stops the run.

The closed-form ULM that makes 1,064 null fits affordable reproduces decoupleR's `run_ulm` on the
observed network at Spearman 1.000000, largest absolute difference 2.1e-14, so the nulls score the
statistic the headline was read off (`tables/ulm_engine_validation.csv`).

The sign-flip's identity configuration — the refit that flips no donor — reproduces the committed
ranked list at Spearman 1.000000 over all 13,999 genes, so the null is measured against the
published contrast. That second figure is a run-time gate message.

## What the two nulls returned

**Holding promiscuity fixed moves the null's centre.** Rewiring the CollecTRI graph by curveball
trades preserves every regulon's size *and* every target's in-degree, so a drawn regulon
oversamples the jointly-owned high-|t| genes exactly as heavily as the observed one does. HIF1A's
null mean is therefore **+5.54** against an observed +8.30. It clears that null at empirical p
0.005, z 2.58, passing the 95th percentile of +7.29 by about one unit, where the same score
cleared the size-and-expression-matched random regulon's 95th percentile of +2.10 by more than
six. Roughly two thirds of the score is what any regulon of that size and promiscuity profile
earns on this contrast.

**The rewiring null singles out no factor.** All eight, ordered by rewiring p
(`tables/rewiring_null.csv`):

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

STAT3 and CREB1 fail this null, and NFKB1 remains indistinguishable from HIF1A. The smaller
regulons doing relatively better than the large ones is the ordering the size-conditional
residual produced by a different route. EPAS1 is the clean negative and the useful one: observed
+2.794 against a null mean of +2.790. That regulon is exactly what its size and promiscuity
predict, which is what a factor with no signal looks like.

**Permuting the design says the contrast is real and says nothing about which factor.** Swapping
the tissue labels within donor and refitting end to end is the exchangeability this paired design
licenses. Six Treg donors carry both arms, so all 2⁶ = 64 configurations are enumerated and the
test is exact — no seed enters it. The observed labelling gives HIF1A the largest score of all 64
(+8.30 against a null maximum of +6.63), which is p = 1/64 = 0.0156, the finest this many donors
can resolve. Seven of the eight focus factors reach that same floor, ATF3 being the exception at
p 0.125 (`tables/signflip_null.csv`).

So the sign-flip answers *is there signal in this contrast at all* — yes, and essentially every
regulon's activity is maximal at the true labels — with no power to separate one factor from
another. The rewiring null is the one with resolving power, and it resolves that HIF1A's margin
over a promiscuity-matched regulon is real, small, shared with its neighbours, and absent for
EPAS1.

## Two ways to misread this stage

**A floor is a weak p-value.** `p_exact_one_sided = 0.0156` for seven factors at once is one bit
of information about the contrast. A smaller p-value here takes more paired donors, and more
computation reaches it in no case.

**The rewiring null is the strongest available and remains short of a claim.** It holds the
network's degree structure fixed, the sharpest constraint on offer, and stays a question about
one curated network's target assignments. This stage is silent on hypoxia, on temperature, and on
HIF1A protein.

---

## Tables

### `tables/rewiring_null.csv` — the null with resolving power

One row per focus factor. `obs_score` and `obs_rank` are the observed signed-network ULM values,
and the `null_*` columns describe the draws. A curveball trade holds the targets two regulons
share and redistributes the rest between them, so every regulon's size and every target's
in-degree stay invariant and only the assignment of targets to factors is randomised. That is why
`null_mean` sits far above zero, and why beating this null means more than beating a random gene
set of matched size.

`mean_target_indeg` is the observed regulon's average target promiscuity, the property this null
holds fixed. `p_empirical` is (draws at or above obs + 1)/(draws + 1), so it returns exactly zero
in no case, and `z_vs_null` standardises the gap. Because size is preserved exactly,
`null_rank_median` and `null_rank_best` share the `n_tfs_scored` denominator with `obs_rank`.

### `tables/rewiring_null_draws.csv` — the per-draw substrate

One row per (draw, factor): 1,000 degree-preserving rewirings scored for each of the 8 focus
factors. `score` is that draw's closed-form ULM value and `rank` its position among every factor
scored in the same draw. Draws come from one Markov chain advanced between samples by trades, so
consecutive rows are decorrelated and dependent by construction. Read it when a summary statistic
needs checking against the shape it came from.

### `tables/signflip_null.csv` — the exact design permutation

One row per focus factor. The two tissue labels are swapped within a donor and the whole
limma-voom contrast is refitted, which is the exchangeability a paired design licenses and which
leaves the correlation between genes intact. `n_paired_donors` donors carry both arms, so
`n_configurations` = 2^that, and every one is enumerated.

`n_ge_obs` counts configurations scoring at or above the observed, the observed one included, and
`p_exact_one_sided` is that count over `n_configurations`. `p_floor` is 1/`n_configurations`, and
`at_resolution_floor` reads TRUE when the p-value is as small as this many donors can make it.
Flipping every donor negates the contrast, so the configurations come in sign-symmetric pairs and
`null_mean` sits near zero. The gene set is filtered once on the observed design and held fixed,
so scores compare across configurations, and voom is recomputed per configuration because its
weights depend on the design.

### `tables/signflip_null_draws.csv` — the per-configuration substrate

One row per (configuration, factor): all 64 within-donor label swaps of the 6 paired donors,
scored for each focus factor. `donors_flipped` names the donors whose labels were swapped,
pipe-delimited, and `is_observed` marks the configuration that flips none. Configurations are
enumerated in a fixed order, so this table is byte-stable across runs. Read it for the null's
shape, and for the observed row being the extreme one for most factors.

### `tables/null_ladder.csv` — the four nulls side by side

One row per (factor, null), ordered weakest-to-strongest by what the null holds constant.
`holds_fixed` names that and `permutes` names what is randomised, which together bound what a
p-value from that row can mean. `obs` is the same observed score in every row for a given factor,
so only the reference distribution changes and the movement in `p_empirical` is the whole content
of the table. `null_q95` is empty for the sign-flip rung, whose 64 configurations support no
stable upper percentile.

Read down a factor's block and the conclusion weakens along the order. HIF1A's null centre climbs
from −0.03 through +0.52 to +5.54 while the observed +8.30 never moves. Each rung is published
beside the others because no single rung supports a claim alone.

### `tables/ulm_engine_validation.csv` — the engine gate

One row. `score_closed_form` against `score_decoupler` on the observed signed network, summarised
by `spearman`, `pearson` and `max_abs_diff`. decoupleR regresses every gene's contrast statistic
on the regulon's mode of regulation, zero for a non-target, so the fit has a closed form costing
one sparse matrix-vector product for all factors at once, which is what makes tens of thousands
of null fits feasible. `gate` reads `pass` only when `spearman` reaches its threshold.

This satisfies the compartment's standing requirement that an engine swap be shown by rank
correlation to be a method change.

### `tables/source_hash_manifest.csv`

The SHA-256 pin on the CollecTRI human regulon table this stage reads. This stage keeps its own
pin, so a moved network stops both this stage and its neighbour.
