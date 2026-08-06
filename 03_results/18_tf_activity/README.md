# 18_tf_activity — Attacking the HIF1A activity score from several angles

`HIF1A` ranks 4 of 592 by NES in the committed unsigned-regulon transcription-factor sweep of the
sorted-Treg synovial-fluid-versus-paired-blood contrast: NES 2.329, pooled FDR 4.93e-14, over 293
targets. The three regulons immediately around it carry 276 to 285 targets. That coincidence is
the reason for this stage.

**The question.** How much of that rank survives when the network, the estimator, the regulon
size, the edge signs and the symbol vocabulary are each varied in turn?

The mouse anchor ran the same forensics on its own contrast, where the equivalent result moved
from rank 1 to 12 on a network swap and from 12 to 142 on a swap of estimator, and the score
turned out to be carried by generic stress and glycolytic genes sitting in many regulons at
once. These tables repeat that examination on the human contrast.

**What was computed.** Thirteen network-by-estimator configurations for eight focus factors, a
size-conditional calibration over every scored regulon, a per-target decomposition of two
regulons, two matched random-regulon nulls, a signed-against-unsigned comparison, and a
three-way symbol ledger with an alias-recovery audit.

**What was drawn.** Six panels: two rank cascades, the size calibration, the promiscuity strip,
the selective-target bars and their volcano.

**Tier.** Annotation throughout. No row reaches [`../master/`](../master/), and the language
stays correlative.

## What an inferred activity is here

**An inferred transcription-factor activity is a statistic computed over target-gene
expression.** decoupleR-ULM regresses every gene's contrast statistic on the mode of regulation
a network assigns it and reports the t-statistic of the slope. The number describes how the genes
the network assigns to a factor behave on this contrast. HIF1A protein, its nuclear localisation
and its occupancy are all untested here.

The statistic inherits every property of the regulon it was computed over — how many targets,
which genes, how many other regulons claim the same genes, and whether the recorded per-edge
signs are used. The object under test throughout is "the CollecTRI HIF1A regulon's CollecTRI-ULM
activity on the Treg synovial-versus-blood contrast", a name that stays checkable.

## The five forensics and what they returned

**The rank cascade holds, and HIF1A is the steadiest of the eight.** Across thirteen
configurations — four network variants crossed with three estimators, plus the committed
unsigned-regulon fgsea column — HIF1A places in the top twelve twelve times
(`tables/hif1a_rank_cascade.csv`).

| Factor | Top-12 placements of 13 | Worst rank | Worst configuration |
|---|---|---|---|
| HIF1A | 12 | 42 of 388 | literature_signed / MLM |
| STAT3 | 10 | 25 of 601 | signed / MLM |
| NFKB1 | 9 | 252 of 601 | unsigned / MLM |
| CREB1 | 9 | 72 of 388 | literature_signed / MLM |
| ATF3 | 4 | 364 of 388 | literature_signed / MLM |
| REL | 4 | 298 of 388 | literature_signed / MLM |
| HSF1 | 3 | 70 of 601 | unsigned / MLM |
| EPAS1 | 0 | 355 of 601 | unsigned / MLM |

The multivariate estimator, the axis that collapsed the murine result, *improves* HIF1A's rank:
6 to 4 under the signed network and 6 to 2 under the unsigned one. On that same axis NFKB1 falls
from ULM rank 7 to MLM rank 138 signed and 252 unsigned, and REL reaches 298 of 388 under
literature-signed MLM. The rank instability the anchor documented is present in this contrast,
falling on the NF-κB and AP-1 members.

**Activity scales with regulon size, and that is the binding constraint.** Spearman ρ between
size and NES is 0.541 in Treg, 0.582 in Tcon, 0.431 in CD8, and between size and ULM score
0.471, 0.482, 0.352 (`tables/regulon_size_spearman.csv`). Permuting the gene labels of the same
ranked list drops the ULM correlation to 0.084, which places the gradient in the breadth of the
synovial-side shift that a larger regulon samples more thoroughly.

Removing the size-conditional expectation taken over the real regulons leaves HIF1A's NES
residual at +0.269, ranking **252 of 592**, with NFKB1 at 232, STAT3 at 272 and CREB1 at 287
(`tables/regulon_size_calibration.csv`). The four largest regulons in the headline table are all
mid-pack on the part of their score that size leaves unexplained, and the smaller regulons rise:
ATF3 to 65, REL to 97, HSF1 to 129. On the signed ULM score, whose size dependence is weaker,
HIF1A's residual rank is 20 of 601.

**The headline's denominator excludes the most promiscuous regulons.** The sweep applies its size
cap of 500 to the raw CollecTRI regulon before intersecting with the ranked list, so nine
regulons sit outside the 592-member family a headline rank is read against. Five of those nine
are the top five of the same ranked list under ULM: SP1 (11.51), JUN (10.80), RELA (9.98), NFKB
(9.86) and AP1 (8.99), all above HIF1A's 8.30 (`tables/fgsea_family_size_cap.csv`).

**The score is carried by jointly-owned genes.** Of HIF1A's 293 tested targets, 27 are claimed by
no other CollecTRI regulon and sum to 0.14% of its signed contribution. The 73 targets sitting in
more than 25 regulons carry 35%, and the mean target sits in 23.6 regulons. Its two largest
contributors are BHLHE40 (19 regulons) and CDKN1A (254).

NFKB1 decomposes the same way — 2 exclusive targets at 0.07%, 30% from targets in more than 25
regulons, mean 26.2 — and three of the ten largest contributors are the same genes for both
factors (NAMPT, SDC4, CCL4). On every axis measured the two are interchangeable: ULM 8.299
against 8.286, NES 2.329 against 2.348, size residual rank 20 against 17.

**Edge signs reshuffle the network and leave HIF1A where it was.** Forcing every edge positive
gives Spearman ρ 0.708 between the signed and unsigned rank vectors over 601 factors, with 536 of
them moving more than ten places — ATF3, at 26% repressing edges among its 72 present targets,
moving from rank 60 to rank 9. HIF1A carries 21 repressing edges among its 293 present targets,
and the swap moves its score from +8.30 to +9.51 while leaving its rank at 6
(`tables/signed_vs_unsigned.csv`). The committed unsigned-regulon NES pools activating and
repressing targets into one gene set, so that headline uses no sign at all.

## That 0.14% is a cancellation

The 27 exclusively-claimed targets are 15% of HIF1A's signed total in *magnitude* — 72.6 of 474.6
units — and net to 0.65, because 13 go up on the synovial-fluid side and 14 go down. Twelve of
the 27 clear the compartment's FDR cut, and those twelve split six each way, so the set carries
evidence and carries no direction. Counted on fold change the split reads 12 and 15: TM9SF4 sits
on a repressing edge, so it contributes positively while going down.

**16 of the 27 carry `default activation`.** CollecTRI records no edge direction for them, so
activation is the direction the arithmetic assumes for those sixteen, and it cites literature
evidence for the remaining 11.

What this licenses is narrow and is a statement about the regulon: the *direction* of HIF1A's
ULM score comes from targets other regulons also claim. That is the bound the size calibration
reaches by another route.

## The canonical targets move in both directions

All eleven named HIF1A-selective canonical targets carry a recorded activating edge
(`tables/canonical_hif1a_targets.csv`). The glycolytic members go up on the synovial-fluid side
(PGK1 t = +9.76, LDHA +6.58, SLC2A1 +3.69, HK2 +3.16), BNIP3L (+3.18) and VEGFA (+2.97) go up,
BNIP3 (+1.29), CA9 (+1.25) and EGLN3 (+0.14) are flat, and PDK1 goes **down** at t = −4.65 with a
gene-level FDR of 0.0065. ADM is absent from the ranked list, dropped by the expression filter.

## The symbol-vocabulary guard, and why it needed its own guard

The JIA count matrix carries pre-2019 HGNC symbols while CollecTRI carries current ones, so a
renamed target fails to join and leaves the regulon smaller with no error raised. All four named
probes confirm the vocabulary: the matrix holds `MB21D1`, `TMEM173`, `MARCH5` and `MRE11A`, and
`CGAS`, `STING1`, `MARCHF5` and `MRE11` are absent from it
(`tables/symbol_vocabulary_probes.csv`).

Of HIF1A's 463 network targets, 293 match, 89 sit in the count matrix and were dropped by the
expression filter, 81 are absent from the count matrix, and 3 are lost to a resolvable rename
(`MMUT`→`MUT`, `ATP5IF1`→`ATPIF1`, `TIGAR`→`C12orf5`). Recovering those across the whole Treg
network adds 394 edges over 124 renamed symbols and 184 factors, moves HIF1A's regulon from 293
to 296 targets, and moves its ULM rank from 6 of 601 to 7 of 603.

**A naive alias join would be worse than the gap it closes.** Many pre-2019 symbols were
reassigned as the official symbol of a *different* gene, so `PGF` resolves to `PIGF`, `THPO` to
`TPO` and `KCNA5` to `HK2`, each of which would attach one gene's expression to another gene's
edge. The recovery therefore rejects any candidate that is the official symbol of another Entrez
id, discarding 495 Treg edges. It rejects a further 100 whose reference symbol reaches more than
one present candidate and 6 ambiguous in `org.Hs.eg.db`, for 601 rejections of 995 candidate
rows. Every accepted and rejected resolution is published in `tables/alias_recovery.csv`.

## Three ways to misread these tables

**A high rank is evidence about the regulon.** Twelve of thirteen configurations put HIF1A in the
top twelve, which says the estimator and the network variant decide no part of the answer. How
much the answer means is decided by the size calibration, where HIF1A sits 252nd of 592 on the
statistic the headline was read off.

**Two regulons this similar cannot be separated by this contrast.** HIF1A and NFKB1 agree to
three significant figures on both statistics and share their largest contributors. A reading that
singles out one of them reads noise between them.

**A random-gene-set null is the weaker of the two nulls here.** HIF1A's observed +8.30 sits far
above the 95th percentile of 1,000 random regulons matched on size and repressing-edge fraction
(+1.62), and above the stricter draw that also reproduces its expression-decile composition
(+2.10), so its targets are more than an arbitrary bag of 293 genes. The informative comparison
is against other *real* regulons of the same size, and there HIF1A is ordinary. Both draws
permute the annotation and leave the design fixed, and neither holds target promiscuity fixed.
[`../19_regulon_nulls/`](../19_regulon_nulls/) closes both gaps, and against regulons that
preserve every target's in-degree the null's centre sits at +5.54.

---

## Figures

### `figures/_overview/tf_rank_cascade.png`

**Eight factors across thirteen network-by-estimator configurations.**
One line per factor, labelled at the right edge, each with its own point shape so two lines of
similar hue stay separable. y, rank by descending activity within a configuration, inverted on a
log scale, so rank 1 sits at the top. The four variants: `signed` uses CollecTRI's recorded
per-edge mode of regulation, `unsigned` forces every edge positive, `literature_signed` keeps
evidence-signed edges only, and `alias_recovered` adds targets resolved to the pre-2019 symbol
this matrix carries. The leftmost column is the committed unsigned-regulon fgsea rank.

HIF1A sits in the top twelve in twelve of the thirteen configurations, its remaining placement
rank 42 of 388 under literature-signed MLM. Ordered by the span of ranks a factor traverses,
HIF1A places 2 of 8 at 40 places, with STAT3 narrowest at 23. The same axes move its neighbours
much further, so the rank instability the mouse anchor documented falls here on the NF-κB and
AP-1 members. Denominators differ between configurations and sit in the source table.
*Source* `tables/_overview/tf_rank_cascade.csv` · `02_analysis/scripts/18_tf_activity_viz.R`.

### `figures/_overview/hif1a_rank_cascade_linear.png`

**The same cascade for HIF1A alone, on a linear rank axis.**
One factor, one line, on the configuration axis and order `tf_rank_cascade` uses, so a column
means the same in both. y, HIF1A's rank among the factors scored in that configuration, linear
and inverted, so a step's height is the size of the rank move. Labels give the rank, the factors
scored and the score behind it. Point colour is that score, and the four estimators share no
scale, so colour compares within an estimator alone.

The rank stays between 2 and 12 in twelve configurations and reaches 42 of 388 in the
thirteenth. On the linear axis the mouse anchor uses for its own cascade, that traverse is nearly
flat, where the murine one runs rank 1 to 12 to 142 and back to 8. The two panels are comparable
in shape alone, since the ranked lists differ in length and the factors scored differ between
configurations.
*Source* `tables/_overview/hif1a_rank_cascade_linear.csv` ·
`02_analysis/scripts/18_tf_activity_viz.R`.

### `figures/_overview/tf_activity_vs_regulon_size.png`

**Activity against regulon size, with the size-conditional expectation drawn.**
One grey point per factor, faceted by statistic. x, targets present in the ranked list on a log
scale. The dashed curve is the size-conditional expectation fitted over the real regulons
themselves, so a point above it is more active than its size alone accounts for. Coloured
labelled points are the headline factors. The open triangle below each is the 95th percentile of
random regulons matched on size, repressing-edge fraction and average-expression decile, and the
stalk spans that percentile to the observed value, so a short stalk means a factor barely beats
a matched bag of genes.

Inferred activity rises with regulon size across every factor tested: ρ = 0.47 between size and
ULM score over 601 factors and 0.54 between size and unsigned-regulon NES over 592 sets, falling
to 0.08 when gene labels are permuted. Every large-regulon factor in the headline table sits on
that gradient. The facets use free y axes, so position relative to the curve and the triangle is
the comparable quantity.
*Source* `tables/_overview/tf_activity_vs_regulon_size.csv` ·
`02_analysis/scripts/18_tf_activity_viz.R`.

### `figures/_overview/tf_target_promiscuity.png`

**Where each regulon's score comes from, target by target.**
One point per regulon target, faceted by factor. x, how many CollecTRI regulons contain that
target, on a log scale, so points at x = 1 belong to this regulon alone. y, the target's signed
contribution — its moderated t multiplied by the edge sign — so positive means the target moves
with the synovial-fluid side and the zero rule separates the directions. Orange marks the
exclusively-claimed targets. The ten largest positive contributors per facet are named in black,
and the three largest exclusive ones in each direction in orange. In-panel text gives the
exclusive share twice, in magnitude and net, then the share from targets in more than 25
regulons.

HIF1A's 27 exclusive targets hold 15% of its signed total in magnitude and net 0.14%, because 13
go up and 14 go down, while its 73 targets in more than 25 regulons carry 35% net. NFKB1's 2
exclusive targets split the same way at 0.07% net, so joint ownership of the directional high-t
genes bounds both regulons equally.
*Source* `tables/_overview/tf_target_promiscuity.csv` ·
`02_analysis/scripts/18_tf_activity_viz.R`.

### `figures/_overview/tf_selective_targets.png`

**The exclusively-claimed targets, named one by one.**
One row per target that no other CollecTRI regulon contains, faceted by factor, ordered by signed
contribution. The bar runs from zero to that contribution, so length is magnitude and side is
direction, with colour restating the direction. Row pitch is equal in both panels, so a bar
length means the same in each. A filled point means CollecTRI records literature evidence for the
edge direction, and an open point means activation was assumed by default. A dashed bar marks a
repressing edge, which flips the contribution's sign away from the gene's own direction.

The targets HIF1A alone claims carry magnitude without direction: 13 of the 27 go up on the
synovial-fluid side and 14 go down. Glycolytic members fall on both sides (PGAM1 +5.68, GBE1
+3.93 up; PFKL −2.55, TKTL1 −4.25 down). 16 of the 27 carry no recorded evidence for the edge
direction. NFKB1's 2 such targets split one each way (GCA +2.80, BST1 −2.47). TM9SF4 is the only
repressing edge here, so its positive contribution comes from a gene that goes down in synovial
fluid.
*Source* `tables/_overview/tf_selective_targets.csv` ·
`02_analysis/scripts/18_tf_activity_viz.R`.

### `figures/_overview/tf_selective_targets_volcano.png`

**The same targets returned to the contrast they were scored on.**
The standard volcano of the committed donor-pseudobulk contrast: x is log2 fold change, synovial
fluid over paired blood, and y is raw p on a −log10 scale, with colour the significance category
decided on FDR. The dashed horizontal rule is the raw p that realises the FDR cut and the
vertical rules the fold-change cut. Only the exclusively-claimed targets are ringed and named,
circles for HIF1A and triangles for NFKB1, the ring unfilled so the point keeps its category
colour.

The spread of the named genes across both halves is the point, so nothing else is labelled. 12
go up in synovial fluid and 15 go down, 12 clear FDR 0.05, and those 12 split 6 each way. That is
what makes the 0.14% share a cancellation: in magnitude the same 27 targets are 15% of HIF1A's
signed total. NFKB1's 2 exclusive targets split one each way and clear the FDR cut in neither
case. The split here is by fold change; by signed contribution it reads 13 up and 14 down.
*Source* `tables/_overview/tf_selective_targets_volcano.csv` ·
`02_analysis/scripts/18_tf_selective_volcano_viz.R`.

---

## Tables

### The activity surface

**`tables/tf_activity_all.csv`** — one row per (population, variant, method, factor), 36
configurations over three populations. `score` is the estimator's own statistic, comparable
within a configuration alone: ULM and MLM report regression t-statistics on different fits while
consensus reports a mean of folded z-scores, and the regulons differ across the `variant` axis
too. `padj` is BH within a configuration, `rank` is by descending score, `pct_rank` normalises by
`n_tfs_scored`, and `regulon_size` is the targets present. Use rank and `n_tfs_scored` together
for any cross-configuration reading.

The top of the signed ULM ranking for Treg is a crowd of promiscuous regulators — SP1 11.51, JUN
10.80, RELA 9.98, NFKB 9.86, AP1 8.99 — with HIF1A sixth at 8.30 and NFKB1 seventh at 8.29.

**`tables/hif1a_rank_cascade.csv`** — one row per (population, factor, configuration) for the
eight focus factors. `configuration` is `variant / method`, and the `unsigned_geneset / fgsea`
row is the committed sweep value re-read here. `n_tfs_scored` is that configuration's
denominator, which differs between variants, so read the pair together and use `pct_rank` across
them.

**`tables/network_variants.csv`** — one row per (population, variant). `n_edges` counts edges
surviving the intersection with that ranking, `n_tfs_ge_minsize` how many factors clear the size
floor, and `n_repressing_edges` is zero for `unsigned` by construction. `hif1a_size` and
`nfkb1_size` are given per row so a rank reads against the regulon it was computed over.

The four variants differ enough to be a real test axis: the signed network scores 601 factors
over 25,135 edges with HIF1A at 293 targets, keeping only evidence-signed edges halves that to
138 targets and 388 factors, and alias recovery raises it to 296 targets and 603 factors.

### The calibration and the nulls

**`tables/regulon_size_calibration.csv`** — one row per factor. `*_size_expected` is the
size-conditional expectation from a loess fit of the statistic on log10 regulon size taken over
the real regulons themselves, `*_size_residual` is observed minus expected, and
`*_size_residual_rank` ranks the residual descending. The two statistics differ in denominator,
because the unsigned-regulon family excludes the nine regulons above the raw-size cap. A residual
is a calibration device and carries no p-value.

**`tables/regulon_size_spearman.csv`** — one row per population plus one label-permuted Treg
control. Each ρ is computed over every factor scored in that configuration. The permuted row
keeps the same ranked-list values and the same network and shuffles only which gene carries which
statistic.

**`tables/size_matched_null.csv`** — the two random-regulon draws behind the figure's triangles,
matched on size and repressing-edge fraction, and additionally on average-expression decile.
Both permute the annotation and hold the design fixed, and neither holds target promiscuity
fixed.

**`tables/fgsea_family_size_cap.csv`** — one row per excluded regulon. `raw_targets` is the size
before intersecting, the quantity the size filter acts on, and `targets_present` the size after.
`in_fgsea_family` is FALSE on every row by construction, and `ulm_rank` places the excluded
factor in the uncapped ULM ranking. The exclusion is skewed in composition, so a rank read
against the capped family is a rank among the survivors.

### The per-target decomposition

**`tables/target_decomposition.csv`** — one row per regulon target for HIF1A and NFKB1. `stat` is
the limma-voom moderated t, `mor` the recorded edge sign, and `contrib = sign(mor) * stat` the
signed contribution. `n_regulons` counts every CollecTRI regulon containing that target,
`selective` is TRUE when no other regulon claims it, and `promiscuity_band` bins `n_regulons`.
`sign_decision` records whether CollecTRI cites a PMID for the direction or applied `default
activation`. Contributions are arithmetic and carry no separate test.

CDKN1A (t = +20.84), second only to BHLHE40, sits in 254 regulons, and of HIF1A's ten largest
contributors only S100A11 sits in fewer than ten.

**`tables/target_decomposition_summary.csv`** — one row per (factor, promiscuity band).
`sum_contrib` is the band's total signed contribution and `pct_of_total_contrib` its share, so
percentages add to 100 across bands within a factor. `mean_n_regulons` makes the top band
interpretable: for HIF1A its members sit in 67 regulons on average. Bands are cumulative upper
bounds by name and disjoint in membership.

**`tables/canonical_hif1a_targets.csv`** — the eleven named canonical targets returned to the
contrast, each with its moderated t, its edge sign and its gene-level FDR.

**`tables/signed_vs_unsigned.csv`** — the rank vectors of the two sign conventions side by side,
with the per-factor move and the repressing-edge counts behind it.

### The vocabulary ledger

| File | What it holds |
|---|---|
| `tables/symbol_vocabulary_probes.csv` | One row per probe pair. `matrix_symbol` is the pre-2019 spelling and `current_symbol` the present HGNC symbol, with four booleans for presence in the ranked list and in the full pre-filter matrix vocabulary. A row reading TRUE then FALSE is a gene any current-symbol reference set drops in silence. |
| `tables/symbol_vocabulary_check.csv` | One row per (population, factor). `n_targets_in_network` against `n_matched`, with three counts partitioning the unmatched remainder: `n_expression_filtered`, `n_absent_from_count_matrix`, `n_alias_recoverable`. The first two are facts about the dataset and the third is a join failure. Every focus regulon loses 27% to 41% of its targets before any statistic is computed. |
| `tables/alias_recovery.csv` | One row per edge whose target failed to join directly. `resolution` takes four values: `accepted`, plus three rejection classes — candidate owned by another gene, multiple candidates present, or a reference symbol ambiguous in `org.Hs.eg.db`. Only accepted rows enter the `alias_recovered` variant. |
| `tables/ranked_list_keycheck.csv` | One row per population. `n_ensembl_like` counts names matching `^ENSG[0-9]{6,}`, and the run stops unless `key` reads `hgnc_symbol`. An Ensembl-keyed list intersects every network at approximately zero and both fgsea and decoupleR report that as an empty result. |
| `tables/source_hash_manifest.csv` | The SHA-256 pin on the CollecTRI human regulon table this stage reads. The first run writes it and every later run verifies, so a changed network halts the analysis. |

### `tables/_overview/<figure stem>.csv`

Six same-stem sources, one per figure, each holding exactly the rows its panel draws.
