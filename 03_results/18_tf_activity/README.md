# 18_tf_activity: artifact captions

_**Abbreviations:** TF = transcription factor, ULM = univariate linear model, MLM = multivariate
linear model, NES = normalised enrichment score, HGNC = HUGO Gene Nomenclature Committee symbol, SF =
synovial fluid, PB = peripheral blood, BH = Benjamini-Hochberg._

`HIF1A` ranks 4 of 592 by NES in the committed unsigned-regulon TF sweep of the sorted-Treg
synovial-fluid-versus-paired-blood donor-pseudobulk contrast: NES 2.329, pooled BH FDR 4.93e-14, over
293 targets tested. The three regulons immediately around it carry 276 to 285 targets. The pooled FDR
moves whenever the pooled family changes size, so it is quoted against the family the run actually saw
— 11,345 tests in Treg.

The mouse anchor investigated this shape in `mouse_anchor/03_results/04_tf`, where the equivalent
murine result moved from rank 1 to 12 on a network swap and from 12 to 142 on a swap of estimator, and
the score turned out to be carried by generic stress and glycolytic genes sitting in many regulons at
once. The tables here run the same forensics on the human contrast.

## What an inferred TF activity is here

**An inferred transcription-factor activity is a statistic computed over target-gene expression.**
decoupleR-ULM regresses every gene's contrast statistic on the mode of regulation a network assigns
it, and reports the t-statistic of the slope. The number describes how the genes the network assigns
to a factor behave on this contrast. Protein activity is untested: nothing here measures HIF1A
protein, its nuclear localisation, or its occupancy.

The statistic inherits every property of the regulon it was computed over — how many targets that
regulon has, which genes they are, how many other regulons claim the same genes, and whether the
recorded per-edge signs are used. The object under test throughout is "the CollecTRI HIF1A regulon's
CollecTRI-ULM activity on the Treg SF-versus-PB contrast", a name that stays checkable.

Everything here is annotation tier. No row reaches `effect_sizes_treg_arthritis.csv` or any
`03_results/master/` accumulator, nothing pools with the donor-pseudobulk claim spine, and the
language stays correlative.

## The four forensics and what they returned

**The rank cascade holds, and HIF1A is the steadiest of the eight.** Thirteen configurations: four
network variants crossed with three estimators, plus the committed unsigned-regulon fgsea column.

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

The multivariate estimator, the axis that collapsed the murine result, *improves* HIF1A's rank: 6 to 4
under the signed network and 6 to 2 under the unsigned one. On that same axis NFKB1 falls from ULM rank
7 to MLM rank 138 signed and 252 unsigned, and REL — which ranks 2 on the committed fgsea column — reaches
rank 298 of 388 under literature-signed MLM. The rank instability the anchor documented is present in
this contrast, on the NF-kB and AP-1 members.

**Activity scales with regulon size, and that is the binding constraint.** Spearman rho between size
and NES is 0.541 in Treg, 0.582 in Tcon, 0.431 in CD8; between size and CollecTRI-ULM score, 0.471,
0.482, 0.352. Permuting the gene labels of the same Treg ranked list drops the ULM correlation to
0.084, which places the gradient in the breadth of the synovial-fluid-side shift that a larger regulon
samples more thoroughly.

Removing the size-conditional expectation taken over the real regulons leaves HIF1A's NES residual at
+0.269, ranking **252 of 592**, with NFKB1 at 232, STAT3 at 272 and CREB1 at 287. The four largest
regulons in the headline table are all mid-pack on the part of their score that size does not already
account for, and the smaller regulons rise: ATF3 to 65, REL to 97, HSF1 to 129. On the signed ULM
score, whose size dependence is weaker, HIF1A's residual rank is 20 of 601.

**The headline's denominator excludes the promiscuous regulons.** The sweep applies its size cap of 500
to the raw CollecTRI regulon before intersecting with the ranked list, so nine regulons sit outside the
592-member family a headline rank is read against. Five of those nine are the top five of the same
ranked list under CollecTRI-ULM: SP1 (score 11.51), JUN (10.80), RELA (9.98), NFKB (9.86) and AP1
(8.99), all above HIF1A's 8.30.

**The score is carried by jointly-owned genes.** Of HIF1A's 293 tested targets, 27 are claimed by no
other CollecTRI regulon and sum to 0.14% of its signed contribution. The 73 targets sitting in more
than 25 regulons carry 35%, and the mean target sits in 23.6 regulons. Its two largest contributors
are BHLHE40 (19 regulons) and CDKN1A (254).

NFKB1 decomposes the same way — 2 exclusive targets at 0.07%, 30% from targets in more than 25
regulons, mean 26.2 — and three of the ten largest contributors are the same genes for both factors
(NAMPT, SDC4, CCL4). On every axis measured the two are interchangeable: ULM 8.299 against 8.286, NES
2.329 against 2.348, size residual rank 20 against 17.

## That 0.14% is a cancellation

The 27 exclusively-claimed targets are 15% of HIF1A's signed total in *magnitude* — 72.6 of 474.6 units
— and net to 0.65 because 13 go up on the synovial-fluid side and 14 go down. Twelve of the 27 clear
the compartment's FDR cut, and those twelve split six each way, so the set carries evidence and
carries no direction. Counted on fold change the split reads 12 and 15:
TM9SF4 sits on a repressing edge, so it contributes positively while going down.

**16 of the 27 carry `default activation`: CollecTRI records no edge direction for them, so activation
is the direction the arithmetic assumes for those sixteen.** CollecTRI cites literature evidence for
the remaining 11.

What this licenses is narrow, and it is a statement about the regulon: the *direction* of HIF1A's
CollecTRI-ULM score comes from targets other regulons also claim. That is the same bound
`tf_activity_vs_regulon_size` reaches by another route.

Three figures carry it, in order: `tf_target_promiscuity` names the extremes of the strip and gives
both shares, `tf_selective_targets` names all 29 across the two factors, and
`tf_selective_targets_volcano` returns the same genes to the contrast, where the two-sidedness reads as
position.

## The canonical targets move in both directions

All eleven named HIF1A-selective targets carry a recorded activating edge. The glycolytic members go up
on the synovial-fluid side (PGK1 t = +9.76, LDHA +6.58, SLC2A1 +3.69, HK2 +3.16), BNIP3L (+3.18) and
VEGFA (+2.97) go up, BNIP3 (+1.29), CA9 (+1.25) and EGLN3 (+0.14) are flat, and PDK1 goes **down** at t
= −4.65 with a gene-level BH FDR of 0.0065. ADM is absent from the ranked list because the expression
filter dropped it.

## Signs matter across the network and leave HIF1A where it was

Forcing every edge positive reshuffles the ranking broadly: Spearman rho between the signed and
unsigned rank vectors is 0.708 over 601 factors, and 536 of them move by more than ten places, ATF3
(26% repressing edges among its 72 present targets) moving from rank 60 to rank 9. HIF1A carries 21
repressing edges among its 293 present targets, 26 among the 463 in the network, and the swap moves its
score from +8.30 to +9.51 while leaving its rank at 6. The committed unsigned-regulon NES pools
activating and repressing targets into one gene set, so that headline uses no sign at all.

## The symbol-vocabulary guard, and why it needed its own guard

The JIA count matrix carries pre-2019 HGNC symbols while CollecTRI carries current ones, so a renamed
target fails to join and leaves the regulon smaller with no error raised. All four named probes confirm
the vocabulary: the matrix holds `MB21D1`, `TMEM173`, `MARCH5` and `MRE11A`, and `CGAS`, `STING1`,
`MARCHF5` and `MRE11` are absent from it.

Of HIF1A's 463 network targets, 293 match, 89 are present in the count matrix and dropped by the
expression filter, 81 are absent from the count matrix, and 3 are silently lost to a rename that
`org.Hs.eg.db` can resolve (`MMUT` to `MUT`, `ATP5IF1` to `ATPIF1`, `TIGAR` to `C12orf5`). Recovering
them across the whole Treg network adds 394 edges over 124 renamed symbols and 184 factors, moves
HIF1A's regulon from 293 to 296 targets, and moves its ULM rank from 6 of 601 to 7 of 603.

A naive alias join would have been worse than the gap it closes. Many pre-2019 symbols were reassigned
as the official symbol of a *different* gene, so `PGF` resolves to the alias `PIGF`, `THPO` to `TPO`
and `KCNA5` to `HK2`, each of which would attach one gene's expression to another gene's regulon edge.
The recovery therefore rejects any candidate that is the official symbol of another Entrez id, which
discards 495 Treg edges. It also rejects 100 edges whose reference symbol reaches more than one
candidate present in the vocabulary and 6 whose reference symbol is ambiguous in `org.Hs.eg.db`, for
601 rejections of 995 Treg candidate rows. Every accepted and rejected resolution is published in
`tables/alias_recovery.csv`.

## Three ways these tables can be misread

**A high rank is evidence about the regulon.** Twelve of thirteen configurations put HIF1A in the top
twelve, which says the estimator and the network variant do not decide the answer. How much the answer
means is decided by `regulon_size_calibration.csv`: after size is accounted for, HIF1A sits 252nd of
592 on the statistic the headline was read off.

**Two regulons this similar cannot be separated by this contrast.** HIF1A and NFKB1 agree to three
significant figures on both statistics and share their largest contributors. A reading that singles
out one of them reads noise between them.

**A random-gene-set null is the weaker of the two nulls here.** HIF1A's observed ULM score of +8.30
sits far above the 95th percentile of 1,000 random regulons matched on size and repressing-edge
fraction (+1.62), and above the stricter draw that also reproduces its expression-decile composition
(+2.10), so its targets are more than an arbitrary bag of 293 genes. The informative comparison is
against other *real* regulons of the same size, which the size-conditional residual reports, and
there HIF1A is ordinary.

Both draws permute the annotation and leave the design fixed, and neither holds target promiscuity
fixed — the property the decomposition above shows carries the score. `03_results/19_regulon_nulls/` closes
both gaps and is where a reader should go before treating any p-value on this page as strong: against
regulons that preserve every target's in-degree the null's centre sits at +5.54, so HIF1A clears it by
about one unit, where the random-regulon draw gave six.

## tables/source_hash_manifest.csv

The CollecTRI human regulon table read across the compartment boundary is pinned at sha256
`4473c918…`, so a change to that network halts this analysis. Without the pin, every number downstream
of it would move in silence.

**How to read:** One row per cross-compartment source. `source_label` is the name this stage refers to
it by, `source_path` is repository-root-relative, and `sha256` is the digest of the bytes actually
read. The first run of the compute script writes the pin; every later run verifies against it and stops
on a mismatch. Verification is the only gate; this table asserts nothing about the network's content.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `verify_source_hash` | `unbiased_enrichment.tf_network.path = ../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` | `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/ranked_list_keycheck.csv

All three ranked lists are keyed by HGNC symbol with zero Ensembl-like identifiers, so the network
joins on the same vocabulary the reference sets use.

**How to read:** One row per sorted population. `n_ensembl_like` counts gene names matching
`^ENSG[0-9]{6,}` and `key` is the resulting verdict. An Ensembl-keyed ranked list intersects every
network at approximately zero, and both fgsea and decoupleR report that as an empty result, so the
failure is silent. The compute script therefore stops unless `key` reads `hgnc_symbol`. This is a
guard.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `read_ranked` + the key-check block | `tf_activity.populations = [treg, tcon, cd8]` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` |

## tables/symbol_vocabulary_probes.csv

All four named symbol pairs resolve the same way: the count matrix carries the pre-2019 symbol
(`MB21D1`, `TMEM173`, `MARCH5`, `MRE11A`) and none of the four current symbols appears anywhere in it,
so a lookup written against a current symbol returns nothing for these genes.

**How to read:** One row per probe pair. `matrix_symbol` is the pre-2019 spelling and `current_symbol`
the present HGNC symbol; the four boolean columns say whether each appears in the primary population's
ranked list and in the full pre-filter count-matrix vocabulary. A row with
`matrix_symbol_in_count_matrix = TRUE` and `current_symbol_in_count_matrix = FALSE` is a gene that any
reference set keyed on current symbols will drop from this dataset in silence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `probes` block | `tf_activity.symbol_vocabulary_probes = [MB21D1/CGAS, TMEM173/STING1, MARCH5/MARCHF5, MRE11A/MRE11]` | `03_results/03_pseudobulk/tables/ranked_treg.tsv`, `03_results/03_pseudobulk/tables/gene_symbols.csv` |

## tables/symbol_vocabulary_check.csv

Every focus regulon loses 27% to 41% of its network targets before any statistic is computed. For
HIF1A in Treg the 170 unmatched targets split into 89 dropped by the expression filter, 81 absent from
the count matrix entirely, and 3 lost to a resolvable rename.

**How to read:** One row per (population, factor). `n_targets_in_network` is the regulon size in
CollecTRI, `n_matched` the size the statistics were computed over, and three counts partition the
unmatched remainder by cause: `n_expression_filtered` are present in the pre-filter count matrix and
removed by `filterByExpr`, `n_absent_from_count_matrix` never appear in this dataset's annotation, and
`n_alias_recoverable` are renames `org.Hs.eg.db` resolves to a symbol the matrix does carry, named in
`alias_recoverable_symbols` in network-side spelling. The first two are facts about the dataset; the third
is a join failure, and the `alias_recovered` variant repairs it.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `vocab` block | `tf_activity.focus_tfs = [HIF1A, NFKB1, STAT3, CREB1, ATF3, REL, HSF1, EPAS1]` | `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/03_pseudobulk/tables/gene_symbols.csv` |

## tables/alias_recovery.csv

Resolving current network symbols against the pre-2019 vocabulary this matrix carries accepts 394 Treg
edges over 124 renamed symbols and 184 factors, and rejects 601. The largest rejected class, 495
edges, is candidates now owned by a different gene: `PGF` reaching `PIGF`, `THPO` reaching `TPO` and
`KCNA5` reaching `HK2` are three of them.

**How to read:** One row per network edge whose target failed to join directly. `reference_symbol` is
CollecTRI's spelling and `matrix_symbol` the candidate this dataset carries. `resolution` takes four
values: `accepted`, a retired symbol of the same Entrez id; `rejected_symbol_belongs_to_another_gene`,
where the candidate is the current official symbol of some other gene and accepting it would attach
one gene's expression to another gene's edge; `rejected_multiple_aliases_in_vocabulary`, where more
than one candidate is present in the ranked list; and
`rejected_reference_symbol_ambiguous_in_org_db`, where the network symbol reaches more than one Entrez
id. Only accepted rows enter the `alias_recovered` variant. `focus_tf` flags the eight factors carried
through every forensic. The map is `org.Hs.eg.db`'s `SYMBOL` and `ALIAS` tables.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `build_alias_map` | `tf_activity.network_variants includes alias_recovered` | `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `org.Hs.eg.db 3.22.0` |

## tables/network_variants.csv

The four network variants differ enough to be a real test axis. On the Treg universe the signed
network scores 601 factors over 25,135 edges with HIF1A at 293 targets; keeping only CollecTRI's
evidence-signed edges halves that to 138 targets and 388 factors; alias recovery raises it to 296
targets and 603 factors.

**How to read:** One row per (population, variant). `n_edges` counts edges surviving the intersection
with that population's ranked list, `n_tfs_ge_minsize` how many factors clear
`thresholds.gsea_min_size` and therefore appear in that configuration's ranking, `n_repressing_edges` zero
for the `unsigned` variant by construction. `hif1a_size` and `nfkb1_size` are given per row so a rank in
`hif1a_rank_cascade.csv` reads against the regulon it was computed over. Denominators differ between
variants, so ranks are comparable only alongside `n_tfs_scored`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `make_variants` + `variant_summary` | `tf_activity.network_variants = [signed, unsigned, literature_signed, alias_recovered]; tf_activity.default_sign_decision = "default activation"; thresholds.gsea_min_size = 5` | `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` |

## tables/tf_activity_all.csv

The full inferred-activity surface, 36 configurations over three sorted populations. The top of the
signed CollecTRI-ULM ranking for Treg is a crowd of promiscuous regulators — SP1 11.51, JUN 10.80, RELA
9.98, NFKB 9.86, AP1 8.99 — with HIF1A sixth at 8.30 and NFKB1 seventh at 8.29.

**How to read:** One row per (population, variant, method, factor). `score` is the estimator's own
statistic, comparable only within a configuration, since ULM and MLM report regression t-statistics on
different fits while consensus reports a mean of folded z-scores; scores are incomparable across the
`variant` axis too, because the regulons differ there. `padj` is BH within a configuration, `rank` is by
descending score, `pct_rank` normalises rank by `n_tfs_scored`, `regulon_size` is the targets present.
Use rank and `n_tfs_scored` together for any cross-configuration reading.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `run_one` (`decouple(statistics = c("ulm","mlm","wsum"), consensus_score = TRUE, minsize = 5)`) | `tf_activity.methods = [ulm, mlm, consensus]; tf_activity.consensus_statistics = [ulm, mlm, wsum]; thresholds.gsea_min_size = 5; thresholds.gsea_seed = 123` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/fgsea_family_size_cap.csv

Nine CollecTRI regulons exceed the sweep's raw-size cap of 500 and are therefore absent from the
592-member family a headline rank is read against. Five of the nine are the top five factors on the
same ranked list under CollecTRI-ULM — SP1, JUN, RELA, NFKB and AP1 — all scoring above HIF1A's 8.30.

**How to read:** One row per excluded regulon. `raw_targets` is the CollecTRI regulon size before
intersecting with the ranked list, the quantity the sweep's size filter acts on; `targets_present` is the
size after intersection. `in_fgsea_family` is `FALSE` for every row by construction, and `ulm_rank` places
the excluded factor in the signed CollecTRI-ULM ranking of the same contrast, where no size cap applies.
The exclusion is skewed in composition — these are the network's most promiscuous regulons — so a rank read
against the capped family is a rank among the survivors.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `CAPPED` block | `thresholds.gsea_max_size = 500; tf_activity.fgsea_database = TF_Targets` | `03_results/14_unbiased_enrichment/tables/gsea_all.csv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/hif1a_rank_cascade.csv

HIF1A's rank sits in the top twelve in twelve of the thirteen network-by-estimator configurations for
Treg. The multivariate estimator that moves it from rank 6 to rank 4 on the signed network moves NFKB1
from 7 to 138, REL from 8 to 108 and ATF3 from 60 to 97 on that same network.

**How to read:** One row per (population, factor, configuration) for the eight focus factors.
`configuration` is `variant / method`; the `unsigned_geneset / fgsea` row is the committed sweep value,
re-read here. `rank` is by descending activity within its configuration and
`n_tfs_scored` that configuration's denominator, which differs between variants, so read the pair
together and use `pct_rank` across them. ULM scores each regulon independently, MLM fits all regulons
jointly, so a factor losing rank under MLM has targets other regulons also claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `CASCADE` block | `tf_activity.focus_tfs = [HIF1A, NFKB1, STAT3, CREB1, ATF3, REL, HSF1, EPAS1]; tf_activity.network_variants = 4; tf_activity.methods = 3` | `03_results/18_tf_activity/tables/tf_activity_all.csv` substrate plus `03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## tables/target_decomposition.csv

Target by target, HIF1A's synovial-fluid-side contribution comes from genes many regulons share.
CDKN1A (t = +20.84), second only to BHLHE40, sits in 254 CollecTRI regulons, and of HIF1A's ten largest
contributors only S100A11 sits in fewer than ten.

**How to read:** One row per regulon target for HIF1A and NFKB1 on the primary population. `stat` is
the target's limma-voom moderated t, `mor` the recorded edge sign, and `contrib = sign(mor) * stat`
the signed contribution, so positive means the target moves with the synovial-fluid side once the
edge sign is applied. `n_regulons` counts every CollecTRI regulon containing that target,
`n_other_regulons` excludes the row's own factor, `selective` is `TRUE` when no other regulon claims
it, and `promiscuity_band` bins `n_regulons` at the configured breakpoints. `sign_decision` records
whether CollecTRI cites a PMID for the edge direction or applied `default activation`. `avg_expr` and
`padj_gene` come from the committed DE table. Contributions are arithmetic and carry no separate test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `DECOMP_TBL` block | `tf_activity.decompose_tfs = [HIF1A, NFKB1]; tf_activity.selective_max_regulons = 1; tf_activity.promiscuity_bands = [1, 2, 5, 10, 25]; tf_activity.primary_population = treg` | `03_results/03_pseudobulk/tables/ranked_treg.tsv`, `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/target_decomposition_summary.csv

Banding the same targets by how many regulons claim them shows the score concentrated in the shared
tail: HIF1A's 27 exclusively-claimed targets carry 0.14% of its signed total while its 73 targets in
more than 25 regulons carry 35%. NFKB1 splits 0.07% against 30%.

**How to read:** One row per (factor, promiscuity band). `n_targets` is the band's membership,
`sum_contrib` its total signed contribution, and `pct_of_total_contrib` that total as a percentage of the
factor's whole signed sum, so percentages add to 100 across bands within a factor. `mean_n_regulons` gives
the average promiscuity inside the band, which is what makes the top band interpretable: for HIF1A its
members sit in 67 regulons on average. Bands are cumulative upper bounds by name (`<=5` means at most five
regulons) and disjoint in membership, each target assigned to the tightest band it satisfies.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `decomp_summary` | `tf_activity.promiscuity_bands = [1, 2, 5, 10, 25]; tf_activity.decompose_tfs = [HIF1A, NFKB1]` | `03_results/18_tf_activity/tables/target_decomposition.csv` substrate |

## tables/regulon_size_calibration.csv

Ranking factors by the part of their activity that regulon size does not already account for moves
HIF1A from rank 4 to rank 252 of 592 on the unsigned-regulon NES and from rank 6 to rank 20 of 601 on
the signed ULM score. NFKB1, STAT3 and CREB1 move the same way, and the smaller regulons ATF3, REL and
HSF1 move up.

**How to read:** One row per factor scored on the primary population. `ulm_score` and `fgsea_nes` are
the two activity statistics with their own ranks. `*_size_expected` is the size-conditional expectation
from a loess fit of the statistic on log10 regulon size taken over the real regulons themselves;
`*_size_residual` is observed minus expected, so a positive residual is activity beyond what a regulon
of that size gets on this contrast. `*_size_residual_rank` ranks the residual descending, with
`n_tfs_*_ranked` as denominator. The two statistics differ in denominator because the unsigned-regulon
family excludes the nine regulons above the raw-size cap. A residual is a calibration device and
carries no p-value.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `fit_resid` + `rank_desc` | `tf_activity.primary_population = treg; thresholds.gsea_max_size = 500; loess span = 0.75, degree = 1` | `03_results/18_tf_activity/tables/tf_activity_all.csv` substrate plus `03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## tables/regulon_size_spearman.csv

Regulon size and inferred activity are strongly rank-correlated on this contrast in all three sorted
populations — NES against size 0.541 Treg, 0.582 Tcon, 0.431 CD8; ULM score against size 0.471, 0.482,
0.352 — and permuting the gene labels of the same Treg ranked list drops the ULM correlation to 0.084.

**How to read:** One row per population plus one row for the label-permuted Treg control. Each Spearman
rho is computed over every factor scored in that configuration, with denominators in `n_tfs_ulm` and
`n_tfs_fgsea`. The permuted row keeps the same ranked-list values and the same network and shuffles only
which gene carries which statistic, so it measures how much size-versus-activity structure survives
with no contrast signal present; the gap between the real and permuted rows is the part of the size
dependence that comes from the breadth of the synovial-fluid-side shift.
`spearman_ulm_abs_score_vs_size` repeats the correlation on the absolute score, which comes out lower
and places the dependence on signed activity.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `SPEAR` block (`cor(..., method = "spearman")` + `run_ulm` on the permuted vector) | `thresholds.gsea_seed = 123; thresholds.gsea_min_size = 5` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/14_unbiased_enrichment/tables/gsea_all.csv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/size_matched_null.csv

Every focus regulon beats a random regulon matched to it on size, on repressing-edge fraction and on
expression-decile composition. HIF1A's observed ULM score of +8.30 sits far above the 95th percentile of
1,000 such draws — +1.62 matching on size alone and +2.10 under the stricter expression-matched draw — so
these regulons are more than arbitrary bags of genes of the right size.

**How to read:** One row per (factor, statistic, null-match mode). `obs` is the observed value;
`null_mean`, `null_sd`, `null_q95` and `null_max` describe the draws. `pct_of_null` is the percentage
of draws below the observed value, `p_empirical` is `(draws at or above obs + 1) / (draws + 1)` so it
never returns exactly zero, and `z_vs_null` standardises the gap.

Two match modes: `size` draws uniformly from the ranked-list universe at the observed size, and
`size_and_expression` additionally reproduces the observed regulon's average-expression decile
composition, the stricter of the two because network targets are well-studied and highly expressed
genes. Both modes keep the observed regulon's mode-of-regulation composition by permuting its own
signs. This null asks whether a regulon beats an arbitrary gene set of matched size; the comparison
against real regulons of the same size is in `regulon_size_calibration.csv`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `null_scores` + `summarise_null` | `tf_activity.null_draws = 1000; tf_activity.null_expression_deciles = 10; thresholds.gsea_seed = 123; thresholds.gsea_min_size = 5; thresholds.gsea_max_size = 500` | `03_results/03_pseudobulk/tables/ranked_treg.tsv`, `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/signed_vs_unsigned.csv

CollecTRI's recorded per-edge signs change little for HIF1A: it carries 26 repressing edges out of 463
in the network, and forcing every edge positive moves its ULM score from +8.30 to +9.51 while leaving
its rank at 6.

**How to read:** One row per focus factor on the primary population. `score_signed` uses CollecTRI's
recorded mode of regulation, `score_unsigned` forces every edge to +1. `delta_score` and `delta_rank` are
signed minus unsigned, so a negative `delta_score` means the recorded signs lowered the score.
`n_repressing_edges` and `pct_repressing_edges` are counted over the targets present, which bounds how
much sign structure can matter for that factor. The committed unsigned-regulon NES pools activating and
repressing targets into one gene set, so it corresponds to the unsigned column.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `SIGNED_VS` block | `tf_activity.network_variants includes signed and unsigned; tf_activity.methods includes ulm` | `03_results/18_tf_activity/tables/tf_activity_all.csv` substrate, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/canonical_hif1a_targets.csv

Named one by one, the canonical HIF1A-selective targets move in both directions on this contrast. PGK1
(t = +9.76), LDHA (+6.58), SLC2A1 (+3.69), BNIP3L (+3.18), HK2 (+3.16) and VEGFA (+2.97) go up on the
synovial-fluid side; BNIP3 (+1.29), CA9 (+1.25) and EGLN3 (+0.14) are flat; PDK1 goes down at t = −4.65
with a gene-level BH FDR of 0.0065; ADM is absent because the expression filter dropped it.

**How to read:** One row per named target. `stat` is the moderated t on the SF-versus-PB contrast and
`direction_in_contrast` restates its sign in words. `in_ranked_list`, `in_count_matrix` and
`in_hif1a_regulon` say whether the gene was available to the statistic and whether CollecTRI assigns it
to HIF1A at all. `unmatched_cause` names why a gene is missing — `expression_filter` when the count matrix
held it and `filterByExpr` removed it, `absent_from_count_matrix` otherwise — and `alias_in_matrix` gives
the pre-2019 spelling where one exists.

`mor` is +1 in every row, so CollecTRI records none of these targets as repressed by HIF1A.
`log2FoldChange`, `padj_gene` and `avg_expr` come from the committed DE table. Listing a target that
failed to match is deliberate: a silent drop is the failure mode this stage exists to detect.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity.R` | `CANON_TBL` block | `tf_activity.canonical_hif1a_targets = [PDK1, BNIP3, BNIP3L, CA9, SLC2A1, VEGFA, EGLN3, ADM, HK2, LDHA, PGK1]; tf_activity.primary_population = treg` | `03_results/03_pseudobulk/tables/ranked_treg.tsv`, `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv` |

## tables/_overview/tf_rank_cascade.csv

The plotted cascade, eight factors across thirteen configurations for the primary population, carries
each rank with the denominator it was taken against.

**How to read:** The source table of `figures/_overview/tf_rank_cascade.png`, one row per (factor,
configuration). Columns are as in `hif1a_rank_cascade.csv`, restricted to the population the figure
draws and ordered by factor then by the configuration order the x axis uses. `pct_rank` is the
comparable quantity across configurations, because `n_tfs_scored` differs between network variants.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.network_variants = 4; tf_activity.methods = 3; figures.width_wide = 13` | `03_results/18_tf_activity/tables/hif1a_rank_cascade.csv` |

## tables/_overview/tf_target_promiscuity.csv

The plotted per-target decomposition for HIF1A and NFKB1, with each target's signed contribution beside
the number of CollecTRI regulons that contain it.

**How to read:** The source table of `figures/_overview/tf_target_promiscuity.png`, one row per
(factor, target). Columns are as in `target_decomposition.csv`, restricted to the two decomposed factors
and sorted by descending `contrib` within a factor, so the head of each block is the set of genes
carrying that factor's synovial-fluid-side score.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.decompose_tfs = [HIF1A, NFKB1]; tf_activity.selective_max_regulons = 1; figures.volcano_label_top = 10` | `03_results/18_tf_activity/tables/target_decomposition.csv`, `03_results/18_tf_activity/tables/target_decomposition_summary.csv` |

## tables/_overview/tf_activity_vs_regulon_size.csv

The plotted size calibration, every scored factor on both statistics with its size-conditional
expectation, joined to the size-and-expression-matched null for the eight focus factors.

**How to read:** The source table of `figures/_overview/tf_activity_vs_regulon_size.png`, one row per
(factor, statistic). `score` is the plotted y value and `size_expected` the loess curve at that factor's
size, so `score - size_expected` is the vertical distance from the curve. `focus_tf` marks the labelled
points. The null columns (`null_mean`, `null_sd`, `null_q95`, `obs`, `pct_of_null`, `p_empirical`,
`z_vs_null`) are populated for the focus factors alone, because the null was drawn for those regulons.
The two statistics sit on different scales and get separate y axes.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.null_draws = 1000; tf_activity.null_expression_deciles = 10; thresholds.gsea_max_size = 500` | `03_results/18_tf_activity/tables/regulon_size_calibration.csv`, `03_results/18_tf_activity/tables/regulon_size_spearman.csv`, `03_results/18_tf_activity/tables/size_matched_null.csv` |

## figures/_overview/tf_rank_cascade.png

HIF1A's inferred-activity rank on the sorted-Treg
synovial-fluid-versus-paired-blood contrast sits in the top 12 in 12
of the 13 network-by-estimator configurations, its remaining
placement rank 42 of 388 under literature-signed MLM. Ordered by the
span of ranks a factor traverses, HIF1A places 2 of 8 at 40 places,
and STAT3 at 23 places is the narrowest. The same axes move its
neighbours much further — NFKB1 from ULM rank 7 to 252 under unsigned
MLM, REL to 298 of 388 under literature-signed MLM — so the rank
instability the mouse anchor documented falls here on the NF-kB and
AP-1 members.

**How to read:** One line per factor, labelled at the right edge, each with its own
point shape so two lines of similar hue stay separable. The y axis is
rank by descending activity within a configuration, inverted on a log
scale, so rank 1 sits at the top. ULM scores each regulon on its own;
MLM fits every regulon jointly, so a factor whose targets are shared
loses rank there. The four variants: `signed` uses CollecTRI's
recorded per-edge mode of regulation, `unsigned` forces every edge
positive, `literature_signed` keeps only evidence-signed edges,
`alias_recovered` adds targets resolved to the pre-2019 symbol this
matrix carries. The leftmost column is the committed unsigned-regulon
fgsea rank. Denominators differ between configurations and sit in the
source table. Annotation tier: an inferred activity is a statistic
over target-gene expression.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.network_variants=[signed, unsigned, literature_signed, alias_recovered]; tf_activity.methods=[ulm, mlm, consensus]; thresholds.gsea_min_size=5` | `03_results/18_tf_activity/tables/hif1a_rank_cascade.csv` |

## figures/_overview/tf_target_promiscuity.png

The direction of HIF1A's CollecTRI-ULM score comes from targets many
other regulons also contain. The 27 of 293 targets HIF1A alone claims
hold 15% of its signed total in magnitude and net 0.14%, because 13
go up on the synovial-fluid side and 14 go down, while the 73 targets
in more than 25 regulons carry 35% net. NFKB1's 2 exclusive targets
split the same way (0.07% net), so joint ownership of the directional
high-t genes bounds both regulons equally.

**How to read:** One point per regulon target, faceted by factor. The x axis counts
how many CollecTRI regulons contain that target, on a log scale, so
points to the right are jointly owned and points at x = 1 belong to
this regulon alone. The y axis is the target's signed contribution,
its moderated t multiplied by the edge sign, so positive means the
target moves with the synovial-fluid side and the zero rule separates
the directions. Orange marks the exclusively-claimed targets. The 10
largest positive contributors per facet are named in black, and the 3
largest exclusively-claimed targets in each direction in orange,
placed above and below the x = 1 strip because 27 labels do not fit
inside it; `tf_selective_targets` names every one. The in-panel text
gives the exclusively-claimed share of the signed total twice, in
magnitude and net, with the up/down split behind the difference, then
the share from targets in more than 25 regulons. Annotation tier: a
contribution is arithmetic on the committed ranked list and carries
no separate test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.decompose_tfs=[HIF1A, NFKB1]; tf_activity.selective_max_regulons=1; tf_activity.primary_population=treg; figures.volcano_label_top=10` | `03_results/18_tf_activity/tables/target_decomposition.csv` |

## figures/_overview/tf_activity_vs_regulon_size.png

Inferred activity rises with regulon size across every factor tested:
Spearman rho = 0.47 between size and CollecTRI-ULM score over 601
factors, 0.54 between size and unsigned-regulon fgsea NES over 592
sets, falling to 0.08 when the gene labels are permuted. That places
the size dependence in the breadth of the synovial-fluid-side shift a
bigger regulon samples more thoroughly, and every large-regulon
factor in the headline table sits on that gradient.

**How to read:** One grey point per factor, faceted by statistic, x = the factor's
targets present in the ranked list on a log scale. The dashed
dark-grey curve is the size-conditional expectation fitted over the
real regulons themselves, so a point above it is more active than its
size alone accounts for. Coloured labelled points are the
headline-table factors, each with its own shape. The open triangle
below each is the 95th percentile of random regulons matched to that
factor on size, on repressing-edge fraction and on average-expression
decile composition, and the stalk spans that percentile to the
observed value, so a short stalk means a factor barely beats a
matched bag of genes. In-panel text gives the size-versus-activity
Spearman correlation, and for the left facet the same correlation
after the gene labels are permuted. The facets use free y axes
because the statistics differ in scale, so position relative to the
curve and the triangle is the comparable quantity. The
unsigned-regulon facet omits the regulons above the sweep's size cap.
Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.null_draws=1000; tf_activity.null_expression_deciles=10; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500` | `03_results/18_tf_activity/tables/regulon_size_calibration.csv + regulon_size_spearman.csv + size_matched_null.csv` |

## figures/_overview/hif1a_rank_cascade_linear.png

HIF1A's rank stays between 2 and 12 in 12 of the 13 configurations
and reaches 42 of 388 in the thirteenth, literature-signed MLM. On
the linear rank axis the mouse anchor uses for its Hif1a cascade,
that traverse is nearly flat, where the murine one runs rank 1 to 12
to 142 and back to 8. The two panels are comparable in shape alone:
the ranked lists differ in length (13,999 genes here) and the factors
scored differ between configurations, both of which sit in the source
table.

**How to read:** One factor, one line, on the configuration axis and order
`tf_rank_cascade` uses, so a column means the same in both. The y
axis is HIF1A's rank by descending activity among the factors scored
in that configuration, linear and inverted, so rank 1 is at the top
and a step's height is the size of the rank move. Labels give the
rank, the factors scored and the score behind it. Point colour is
that score; the four estimators share no scale, so colour compares
within an estimator alone. Annotation tier: an inferred activity is a
statistic over target-gene expression.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.primary_population=treg; tf_activity.network_variants=[signed, unsigned, literature_signed, alias_recovered]; tf_activity.methods=[ulm, mlm, consensus]` | `03_results/18_tf_activity/tables/hif1a_rank_cascade.csv + ranked_list_keycheck.csv` |

## figures/_overview/tf_selective_targets.png

Named one by one, the targets HIF1A alone claims carry magnitude
without direction: 13 of the 27 go up on the synovial-fluid side and
14 go down, holding 15% of the regulon's signed total in magnitude
and netting 0.14%. Glycolytic members fall on both sides (PGAM1
+5.68, GBE1 +3.93 up; PFKL -2.55, TKTL1 -4.25 down). 16 of the 27
carry no recorded evidence for the edge direction, so activation is
assumed for those 16. NFKB1's 2 such targets split one each way (GCA
+2.80, BST1 -2.47).

**How to read:** One row per target that no other CollecTRI regulon contains, faceted
by factor, ordered by signed contribution. The bar runs from zero to
the target's signed contribution, its moderated t multiplied by the
edge sign, so length is magnitude and side is direction; colour
restates the direction. Row pitch is equal in both panels, so a bar
length means the same thing in each. A filled point means CollecTRI
records literature evidence for the edge's direction; an open point
means the direction was assumed activating by default, which holds
for 16 of HIF1A's 27. A dashed bar marks a repressing edge, which
flips the contribution's sign away from the gene's own direction;
TM9SF4 is the only one here, so its positive contribution comes from
a gene that goes down in synovial fluid. In-panel text gives the
set's size, its up/down split, and its share of the factor's signed
total in magnitude and net, which differ because the set nearly
cancels. Annotation tier: a contribution is arithmetic on the
committed ranked list and carries no separate test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_activity_viz.R` | `save_overview` | `tf_activity.decompose_tfs=[HIF1A, NFKB1]; tf_activity.selective_max_regulons=1; tf_activity.primary_population=treg; tf_activity.default_sign_decision=default activation` | `03_results/18_tf_activity/tables/target_decomposition.csv` |

## figures/_overview/tf_selective_targets_volcano.png

Placed on the contrast they were scored on, the 27 targets HIF1A
alone claims sit on both sides: 12 go up in synovial fluid and 15 go
down, 12 clear FDR 0.05, and those 12 split 6 up against 6 down, so
the set carries evidence and carries no direction. That is what makes
its 0.14% share of HIF1A's signed contribution a cancellation: in
magnitude the same 27 targets are 15% of that total. NFKB1's 2
exclusively-claimed targets split one each way, and the FDR cut is
cleared by neither.

**How to read:** The standard volcano of the committed donor-pseudobulk contrast: x is
log2 fold change, synovial fluid over paired blood; y is raw p on a
-log10 scale; colour is the significance category, decided on FDR
while the axis keeps raw p for resolution. The dashed horizontal rule
is the raw p that realises the FDR cut, the vertical rules the
fold-change cut. Only the targets no other CollecTRI regulon claims
are ringed and named, circles for HIF1A and triangles for NFKB1, the
ring unfilled so the point inside keeps its category colour. The
spread of the named genes across both halves is the point, so nothing
else is labelled: count how many sit either side of zero and how many
clear the rules. The split here is by fold change; by signed
contribution it reads 13 up and 14 down, since TM9SF4 sits on a
repressing edge and contributes positively while going down.
Annotation tier: this restates a committed DE table and a set
membership, so nothing here is a new test.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_tf_selective_volcano_viz.R` | `save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; tf_activity.decompose_tfs=[HIF1A, NFKB1]; tf_activity.selective_max_regulons=1; tf_activity.primary_population=treg` | `03_results/03_pseudobulk/tables/de_SFvsPB_treg.csv, 03_results/18_tf_activity/tables/target_decomposition.csv` |

