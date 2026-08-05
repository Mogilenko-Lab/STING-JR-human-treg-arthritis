# 13_arm_decomposition: artifact captions

_**Abbreviations:** HGNC = HUGO Gene Nomenclature Committee symbol, HSR = heat-shock response, UPR = unfolded-protein response, IFN = interferon, SF = synovial fluid, PB = peripheral blood._

The mouse 39 °C anchor hands this compartment three up arms in human projection space, and the
larger of them is scored, purged and read against curated lenses elsewhere in these results.
None of that tells me what the arms are **made of**. That is a separate and much plainer
question, and this is where I answer it: which curated, versioned, anchor-independent programs
*contain* the genes of each mouse-derived up arm, gene by gene, and how much of each arm no
lens contains at all.

## Membership is not enrichment

Containment and enrichment are different measurements with different failure modes, which is
why this sits apart from the scoring results. A lens containing a gene is
arithmetic over two committed text files: it is exact, it has no p-value, and it cannot be
underpowered. A lens *enriching* in a ranked list is a statistical claim about a contrast,
answered on donor-level pseudobulk. Reading the first as evidence for the second is an
error, so nothing here carries an NES, an FDR, a direction, a
confidence interval or an effect size, and no row reaches `effect_sizes_treg_arthritis.csv` or
any `03_results/master/` accumulator. A band on the figure means "these genes are in that list",
and it means nothing more.

The corollary matters as much: a lens containing **few** genes of an arm reports a small
containment count, which is a fact about composition. I report the small containments
with their genes named rather than folding them into an "other" bucket.

## What the four arms are

Every arm is named by how it was derived, and each carries the mouse anchor's own gate. The
sizes are read from the anchor's frozen `manifest.csv` at run time, so a change to that
contract stops this analysis instead of silently changing a denominator.

| Arm | Mouse contrast | Gate | Genes | Frozen file |
|---|---|---|---|---|
| `WT_heat_up` | WT: 39 vs 37 °C | `fdr_logfc` | 199 | `signatures/WT_heat/WT_heat_up.txt` |
| `KO_heat_up` | cGAS-KO: 39 vs 37 °C | `fdr_logfc` | 218 | `signatures/KO_heat/KO_heat_up.txt` |
| `Interaction_up` | heat × genotype interaction | `fdr_logfc` | 7 | `signatures/Interaction/Interaction_up.txt` |
| `Interaction_up_fdrOnly` | the same interaction, relaxed gate | `fdr_only` | 18 | `signatures/Interaction/Interaction_fdrOnly_up.txt` |

Down arms are out of scope here by decision, not by oversight: no down arm is read, tallied or
drawn anywhere in these artifacts.

## What the nine lenses are

Nothing was curated for this analysis. Each lens is read verbatim from a frozen file, and none
of them knows anything about the mouse anchor, so an intersection measures something rather
than restating a result. The first seven are byte-identical to the ones used for the
`11_heat_decomposition` tables, which is what lets the two sets of tables join on `program`.

| Program key | Curated set | n | Source |
|---|---|---|---|
| `hypoxia` | `HALLMARK_HYPOXIA` | 200 | frozen MSigDB Hallmark |
| `nfkb_tnfa` | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 | frozen MSigDB Hallmark |
| `inflammatory` | `HALLMARK_INFLAMMATORY_RESPONSE` | 200 | frozen MSigDB Hallmark |
| `t_activation` | `HALLMARK_IL2_STAT5_SIGNALING` | 199 | frozen MSigDB Hallmark |
| `ifn_type_i` | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 | frozen MSigDB Hallmark |
| `upr_er` | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 | frozen MSigDB Hallmark |
| `hsr_curated` | `HSR_core` | 56 | frozen curated HSR core, Reactome/GO-derived |
| `sting_specific_published` | `de_Cevins_sting_specific_up` | 21 | published de Cevins et al. 2023 Table S6 |
| `ifn_generic_axis` | `ifn_only_up` | 200 | the STING positive-control compartment's generic type-I IFN axis |
| `unassigned` | the remainder, contained by none of the above | n/a | n/a |

The two axes from the STING positive-control compartment are the reason the panel widened, and
what they are is load-bearing for how their counts read. `sting_specific_published` is the
**published** 21-gene de Cevins Table S6 IFN-independent STING-activation signature: the genes
most specific to the SAVI disease-associated monocyte cluster after every type-I IFN transcript
and every IFN-β-inducible gene was removed. It is a published gene list whose own IFN-β
time-course validation in that compartment is an objective failure at three donors, directionally separating and underpowered, so a positive read on it is
consistent with STING pathway activity and is never proof of it.

`ifn_generic_axis` is the 200-gene up half of the generic type-I interferon axis, the genes
most induced by IFN-β at 24 h in three healthy donors. Overlap with it reflects a generic
type-I interferon response, correlative with but not diagnostic of STING activation. Reading
either axis as a mechanism
label would be exactly the overclaim both were frozen to prevent.

## Three ways these tables can be misread

**The bands are not a partition of genes.** The lenses overlap, so a gene can be contained by
several of them. `arm_program_gene.csv` therefore emits one row per (arm, program, gene), so a
gene in three lenses appears three times, and it carries the multiplicity on every row so the
sharing is never invisible. Of the 199 `WT_heat_up` genes, 67 are contained by at least one
lens, and those 67 carry 100 memberships between them, with as many as four lenses on a single
gene; adding the per-lens counts therefore over-counts by 33 and shrinks the 132-gene
remainder, which is the largest single part of the arm. Forcing a priority-ordered disjoint
assignment would silently decide which program gets credit for a shared gene, so the sharing is
published instead, per gene in `arm_program_multiplicity.csv`.

**Two accountings, and only one of them sums.** `n_intersect` and `frac_of_arm` are the plain
duplicated counts, which do not sum to the arm. `weight_fractional = 1/n_programs_for_gene`
gives every gene one unit of mass split evenly across the lenses that contain it, so summing
that column over an arm returns the arm size exactly. Both are published because either alone
misleads: the counts imply a partition that does not exist, and the fractional shares hide how
many genes a lens actually contains. The figure draws widths from the fractional accounting and
prints the duplicated counts inside the bands, which is why the widths and the numbers disagree.

**The arms are not independent.** The mouse anchor's contract records the three contrasts as
linearly dependent by construction: WT_heat = KO_heat + Interaction. That algebra concerns the
model coefficients and does not carry to the thresholded gene lists, so the arms are neither
independent nor set sums of one another, and the dependence is read by counting shared
genes.

`WT_heat_up` and `KO_heat_up` share 182 genes, 182 of 199 and 182 of 218, so their
near-identical composition is expected structure and the two rows carry one observation
between them.
`Interaction_up` shares **zero** genes with either, because a gene can pass the interaction
gate while failing both main-effect gates; and `Interaction_up_fdrOnly` contains all 7
`Interaction_up` genes among its 18, so those two rows are one contrast read at two gates.
Every one of these counts is derivable from the `gene` column of `arm_program_gene.csv` and is
printed by `02_analysis/scripts/13_arm_decomposition.py` on each run.

## What the membership returns

The two large arms have the same shape and it is mostly unclaimed. Nine curated lenses contain
67 of the 199 `WT_heat_up` genes and 68 of the 218 `KO_heat_up` genes, leaving remainders of 132
and 150, in both cases the largest single part. What the lenses do claim is dominated by
inflammatory gene content: 35 `WT_heat_up` genes sit in `HALLMARK_TNFA_SIGNALING_VIA_NFKB` and
21 in `HALLMARK_INFLAMMATORY_RESPONSE`, against 18 in `HALLMARK_HYPOXIA` and 14 in
`HALLMARK_IL2_STAT5_SIGNALING`.

**Why this stage says 132 unassigned where `11_heat_decomposition` says 137.** Both numbers are
correct: they are the same arm measured against different lens
panels. `11_heat_decomposition` reads `WT_heat_up` against **seven** lenses and leaves 137 genes
unclaimed. This stage adds two more, the 21 published IFN-independent STING genes and the
200-gene generic type-I interferon axis, and those two claim exactly **5 further genes**
(`IL1RN`, `IL2`, `OAS3`, `SAMD4A`, `TGM1`), giving 132. Restricting `arm_program_gene.csv` to the
original seven programs reproduces 137 exactly, and `13_arm_decomposition.py` asserts that
against the committed `11_heat_decomposition` table on every run rather than leaving it to a
reader to check. A remainder is only meaningful against the panel it was measured on: expect it
to shrink as lenses are added, and read the lens panel alongside every remainder quoted anywhere
in this project. The arm size is 199 either way.

Three containments are small enough to be worth naming outright, since a reader scanning the
figure will not see them. The curated HSR core contains 2 of the 199 `WT_heat_up` genes
(`HSPA1A`, `HSPH1`); `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` contains 1 (`ATF3`); and
`HALLMARK_INTERFERON_ALPHA_RESPONSE` contains 1 (`PROCR`). Of the 21 published
IFN-independent STING genes, 2 fall in `WT_heat_up` (`PLAUR`, `PTGS2`) and 1 in `KO_heat_up`
(`PLAUR`), and the generic type-I interferon axis contains 6 and 4 respectively. So the gene
content of the two large arms is largely inflammatory, is minimally proteostatic by these
curated definitions, and is not STING-specific by composition.

The two Interaction arms invert that shape entirely, and their thinness is the first thing to
note about them. In the 7-gene `Interaction_up`, `HALLMARK_INTERFERON_ALPHA_RESPONSE` contains
4 (`IRF7`, `MX1`, `RTP4`, `TRIM5`) and the generic type-I interferon axis 2, while
`HALLMARK_INFLAMMATORY_RESPONSE` contains 2 and every hypoxia, activation, proteostasis and
published-STING lens contains none. At the relaxed `fdr_only` gate the pattern holds and
sharpens: 10 of 18 in Hallmark type-I interferon, 6 in the generic axis, still 0 in the
published IFN-independent STING set.

Read as composition and nothing else, that says the genes the mouse interaction contrast
selected are predominantly
interferon-annotated, and that the annotation available for them is the generic interferon kind
rather than the IFN-independent STING kind. It says nothing about whether any of them separates
synovial fluid from paired blood in human data; at 7 and 18 genes these arms are thin sets, and
that question belongs to the donor-pseudobulk panels.

## A second thing a containment count depends on: the vocabulary

Everything above is counted over the frozen lists in full, which is the right universe for a
composition question. It is not the universe an enrichment statistic lives in. A ranked list
carries only the genes that were detected and survived `filterByExpr`, and this compartment's
count matrix is frozen to a CellRanger hg19 HGNC vintage while the reference sets ship current
symbols, so a set meeting that ranking loses members two ways at once — to power and to
nomenclature.

Both losses are large enough to change a number a reader would quote. `HALLMARK_HYPOXIA` is 200
genes as curated; against the Treg ranked list 139 match by exact symbol and 4 more only after
their current symbols are resolved into the vintage the matrix carries (`CAVIN1` as `PTRF`,
`CAVIN3` as `PRKCDBP`, `ERO1A` as `ERO1L`, `NOCT` as `CCRN4L`), giving 143. `WT_heat_up` is 202
as curated; against the same ranking 119 match by exact symbol and one more only by alias
(`DYNLT2B` as `TCTEX1D2`), giving 120. The 18-gene containment
above is therefore 12 on that ranking, and the six genes that fall out — `ADM`, `ADORA2B`,
`CCN1`, `EGFR`, `F3`, `TGM2` — are absent for reasons the tables separate rather than pool.

`figures/_overview/arm_hypoxia_euler.png` draws both universes side by side for exactly this
reason: 18 and 12 are both right, and either quoted without its vocabulary is wrong. Any
hypoxia set size taken from elsewhere in this project should be checked for the same thing,
because a table written before the alias map was built carries the exact-match count and reads
like the whole set.

## tables/arm_program_gene.csv

One row per (arm, program, gene) shows that the 199 `WT_heat_up` genes generate 232 rows, 100
memberships across nine curated lenses plus 132 unclaimed genes, so the arm's curated
composition is dominated by genes no lens names.

**How to read:** The alluvial substrate: `arm` → `program` → `gene`, one row per membership. A
gene contained by k lenses appears k times, each row carrying `n_programs_for_gene = k` and
`weight_fractional = 1/k`; summing `weight_fractional` within an arm returns the arm size
exactly, while counting rows over-counts by the amount of sharing. A gene no lens contains gets
exactly one row with `program = unassigned`, `curated_set = (none)`,
`n_programs_for_gene = 0` and `weight_fractional = 1`, so the remainder is a first-class member
of the same table rather than an absence from it. `gate` repeats the mouse anchor's gate for
that arm (`fdr_logfc`, or `fdr_only` for the relaxed Interaction variant). Rows are ordered by
arm, then by the declared lens order, then by symbol. Membership only: containment of a gene
by a lens, with no NES, direction or effect size anywhere in the file. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition.py` | `arm_program_gene` | `arms = 4 mouse-derived up arms at gates fdr_logfc + fdr_only; lenses = 9 frozen curated sets; measurement=membership_not_enrichment; evidence_tier=secondary_annotation` | `../mouse_anchor/03_results/human_projection/signatures/{WT_heat/WT_heat_up.txt,KO_heat/KO_heat_up.txt,Interaction/Interaction_up.txt,Interaction/Interaction_fdrOnly_up.txt}`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt`, `../sting_positive_control/03_results/06_reference_axis/signatures/{sting_specific_up.txt,ifn_only_up.txt}` |

## tables/arm_program_summary.csv

Curated lenses contain 67 of the 199 `WT_heat_up` genes and 68 of the 218 `KO_heat_up` genes,
with TNFA/NF-kB the largest single containment at 35 and 33, while the thin Interaction arms are
predominantly interferon-annotated (4 of 7 and 10 of 18 in Hallmark type-I interferon) and carry
no published-STING gene at all.

**How to read:** One row per (arm, program), plus one `unassigned` row per arm. `n_intersect` is
how many of that arm's genes the curated set contains and `frac_of_arm` its share of `n_arm`;
`n_curated_set` is the lens size, so a small `n_intersect` against a large lens is a statement
about the arm and not about the lens. `genes` is the semicolon-delimited, alphabetically sorted
membership, and it is empty where the lens contains nothing of that arm, which is the frozen
state for that pair. Rows within an arm **overlap**, so `n_intersect` does
not sum to `n_arm`; use `weight_fractional` in `arm_program_gene.csv` for an accounting that
does. Membership only: no NES, no p-value, no effect size. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition.py` | `arm_program_summary` | `arms = 4 mouse-derived up arms at gates fdr_logfc + fdr_only; lenses = 9 frozen curated sets; evidence_tier=secondary_annotation` | `../mouse_anchor/03_results/human_projection/signatures/{WT_heat/WT_heat_up.txt,KO_heat/KO_heat_up.txt,Interaction/Interaction_up.txt,Interaction/Interaction_fdrOnly_up.txt}`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt`, `../sting_positive_control/03_results/06_reference_axis/signatures/{sting_specific_up.txt,ifn_only_up.txt}` |

## tables/arm_program_multiplicity.csv

26 of the 67 claimed `WT_heat_up` genes are contained by two or more curated lenses: 21 by two,
3 by three, and `ATF3` and `PLAUR` by four each, which is why the per-lens counts
cannot be summed.

**How to read:** One row per (arm, gene): every gene of every arm, whether or not a lens
contains it. `n_programs` counts the lenses containing that gene and `programs` names them,
semicolon-delimited in the declared lens order; a gene no lens contains reads `n_programs = 0`
and `programs = unassigned`, so the remainder is visible here as rows rather than as silence.
This is the audit trail behind the overlapping assignment: filter to `n_programs >= 2` to see
exactly which genes two bands of the figure share, and how independent any two containment
counts really are. A plain re-expression of `arm_program_gene.csv` with no statistic in it.
Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition.py` | `arm_program_multiplicity` | `arms = 4 mouse-derived up arms; lenses = 9 frozen curated sets; evidence_tier=secondary_annotation` | `../mouse_anchor/03_results/human_projection/signatures/{WT_heat/WT_heat_up.txt,KO_heat/KO_heat_up.txt,Interaction/Interaction_up.txt,Interaction/Interaction_fdrOnly_up.txt}`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt`, `../sting_positive_control/03_results/06_reference_axis/signatures/{sting_specific_up.txt,ifn_only_up.txt}` |

## tables/_overview/arm_program_composition.csv

The plotted table pairs each band's duplicated gene count with its fractional width, which is
where the two accountings visibly diverge: 35 `WT_heat_up` genes in TNFA/NF-kB render as a
12.1% band rather than a 17.6% one, because 20 of those 35 are also contained by another lens
and so contribute a fraction of a gene each.

**How to read:** One row per plotted band, ordered by arm and then by `band_order` (lenses by
total fractional share across all arms, widest first, remainder always last). `n_intersect` and
`frac_of_arm` are the duplicated accounting copied from `arm_program_summary.csv`, the numbers
printed inside the bands. `weight_fractional_sum` and `frac_of_arm_fractional` are the fractional
accounting summed from `arm_program_gene.csv`, the band widths, which total exactly 1.0 per
arm. The `arm_*` columns repeat that arm's not-a-partition counts on every row so a caption can
quote them without a second file: `arm_n_claimed`, `arm_n_unassigned`, `arm_n_claims_total`,
`arm_n_excess_claims` (how much summing the per-lens counts over-counts) and
`arm_max_lenses_per_gene`. `is_partition` is `False` on every row by construction, and
`measurement = membership_not_enrichment` states the tier in the data rather than only in prose.
Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition_viz.py` | `composition_table` | `accounting=fractional (1/n_programs_for_gene) for widths, duplicated for printed counts; evidence_tier=secondary_annotation` | `03_results/13_arm_decomposition/tables/arm_program_summary.csv`, `03_results/13_arm_decomposition/tables/arm_program_gene.csv`, `03_results/13_arm_decomposition/tables/arm_program_multiplicity.csv` |

## figures/_overview/arm_program_composition.png

9 curated anchor-independent lenses contain 67 of the 202 WT_heat_up
genes and 68 of the 221 KO_heat_up genes, leaving remainders of 135
and 153 genes as the largest single part of each large arm.
Inflammatory gene content dominates what the lenses do claim: 35
TNFA/NF-kB and 21 inflammatory-response genes in WT_heat_up, against 2
in the curated HSR core and 2 of the 21 published IFN-independent
STING genes. The thin Interaction arms invert that shape: Hallmark
type-I interferon contains 4 of the 7 genes at the fdr_logfc gate and
10 of the 19 at the relaxed fdr_only gate.

**How to read:** Each band counts how many of one arm's genes a frozen curated lens
contains. That is set arithmetic over committed gene lists, so no NES,
FDR, direction or effect size appears here or in these tables. One row
per arm, named by how it was derived and labelled with its gene count
and the mouse anchor's gate (fdr_logfc, with fdr_only the relaxed
Interaction variant frozen as Interaction_fdrOnly_up.txt). Band width
is the fractional share: a gene in k lenses gives 1/k to each, so
widths total 1.0. The number in a band is the duplicated count of that
arm's genes the lens contains, so numbers and widths measure different
things. A band too narrow for a digit carries its count to the right
of the bar, so every lens keeps its number. Grey is the remainder, the
genes no lens contains, left unnamed on purpose. The lenses overlap,
so the bands are a containment tally: 67 of the 202 WT_heat_up genes
are contained by at least one lens, and those 67 carry 100
memberships, with up to 4 lenses on one gene, so the printed counts
exceed the claimed genes by 33. Per-gene multiplicity is in
arm_program_multiplicity.csv. The four rows share structure by
construction: the mouse contrasts are linearly dependent as model
coefficients (WT_heat = KO_heat + Interaction), and the two
Interaction rows are one contrast at two gates, so agreement between
rows is expected. That algebra holds for the coefficients and stops at
the thresholded lists: WT_heat_up and KO_heat_up share 185 genes,
Interaction_up shares 0 with either, and Interaction_up_fdrOnly holds
all 7 Interaction_up genes among its 19. Annotation tier, firewalled
from the donor-pseudobulk claim spine; no effect-size row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition_viz.py` | `plot_composition` | `colors.okabe_ito + colors.diverging.up; label_floor=0.045; accounting=fractional (1/n_programs_for_gene); evidence_tier=secondary_annotation` | `03_results/13_arm_decomposition/tables/arm_program_summary.csv, 03_results/13_arm_decomposition/tables/arm_program_gene.csv, 03_results/13_arm_decomposition/tables/arm_program_multiplicity.csv` |

## figures/_overview/arm_hypoxia_euler.png

Curated hypoxia accounts for a small minority of the mouse 39
°C-derived up arm in either vocabulary, and the vocabulary decides how
small: the frozen lists share 18 genes of the arm's 202 against 182
hypoxia genes the arm does not carry, while restricting both to the
Treg donor-pseudobulk ranked list leaves only 120 of the arm's 202
genes testable and drops the shared count to 12 — so the arm's hypoxia
content is 9% of the curated arm but 10% of the part of it this
contrast can actually test, and 4 of the 143 testable hypoxia genes
are visible only because alias resolution recovered them.

**How to read:** The same two gene lists, read in two vocabularies. Every area equals
its gene count, solved numerically; the largest residual across both
panels is 9.9e-14 genes. A two-set Euler is exactly solvable for every
valid configuration, because the shared area falls continuously from
the smaller set's size to zero as the circles part, so every area here
is exact; the configurations with no exact solution begin at three
sets. Both panels share one area-per-gene scale and one bounding box,
so the right panel is smaller because it holds fewer genes. Orange is
the mouse WT iTreg 39-versus-37 °C up arm in human projection; blue is
frozen MSigDB Hallmark hypoxia, curated independently of the anchor,
so the overlap is an independent measurement. Each region carries its
count, the shared one inside the lens with its name above on a grey
leader because the lens is too narrow for both. Left is the frozen
lists in full, 202 and 200 genes, the universe the composition bar and
the committed membership tables report. Right keeps only what the Treg
donor-pseudobulk ranked list carries, the universe an enrichment
statistic on that ranking is computed over: the arm falls to 120 and
hypoxia to 143. So the panels report 18 and 12 shared genes and both
are correct, each within its own vocabulary; quoting either without
its vocabulary is the misreading this figure exists to prevent. The 6
shared genes the left panel carries and the right panel drops are ADM,
ADORA2B, CCN1, EGFR, F3, TGM2. Absence has three causes and the source
table splits them: a symbol in the count matrix and outside the
ranking was dropped by filterByExpr, a symbol in the CellRanger
reference and outside the matrix was never detected in sorted T cells,
and a symbol outside the reference is a vocabulary miss. Alias
resolution runs first and only ever adds, so 139 hypoxia genes match
exactly and 4 more only once their current symbols resolve into the
hg19-vintage vocabulary this matrix carries (CAVIN1->PTRF,
CAVIN3->PRKCDBP, ERO1A->ERO1L, NOCT->CCRN4L) — which is what lifts the
testable size from 139 to 143. The arm recovers 1 of its own the same
way (DYNLT2B->TCTEX1D2), which is what lifts it from 119 to 120. This
is membership: the figure and its table carry no NES, FDR, direction
or effect size, and no row reaches effect_sizes_treg_arthritis.csv or
any 03_results/master/ accumulator. A small overlap bounds how much of
the arm is hypoxia gene content. Whether temperature and hypoxia are
separable in this niche is undecidable from cross-sectional human
data. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition_viz.py` | `plot_arm_hypoxia_euler` | `arm = WT_heat_up; lens = HALLMARK_HYPOXIA (gene_sets.project_frozen); vocabulary = symbol_alias.ranked_list at population=treg; alias pairs from symbol_alias.map_path, accepted only; colours = colors.okabe_ito.orange + the hypoxia band hue of this stage's program palette; fill_alpha=0.45; evidence_tier=secondary_annotation` | `03_results/13_arm_decomposition/tables/arm_program_summary.csv, 03_results/13_arm_decomposition/tables/arm_program_gene.csv, 00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt, 00_data/references/symbol_alias/symbol_alias_map.csv, 03_results/03_pseudobulk/tables/ranked_treg.tsv, 03_results/03_pseudobulk/tables/gene_symbols.csv, 03_results/00_build/tables/reference_feature_symbols.csv` |

## tables/_overview/arm_hypoxia_euler.csv

The plotted regions with their gene names, and the ledger behind the
two vocabularies: 139 of the 200 frozen hypoxia genes match the Treg
ranked list by exact symbol and 4 more only after alias resolution, so
what looks like a 139-gene set is a 143-gene one.

**How to read:** One row per (`universe` x `region`), six rows. `region` is `arm_only`,
`shared` or `lens_only` and `n_genes` is the count the figure draws
that area to; `genes` names them, semicolon-delimited and sorted, so
any region can be checked gene by gene. `universe` is `frozen_sets`
for the lists as curated and `treg_ranked_list` for the restriction,
and `vocabulary` records the file the restriction was made against.
The `arm_*` and `lens_*` columns repeat that universe's ledger on
every row: `n_nominal` is the curated size, `n_exact_match` how many
symbols match the ranked list verbatim, `n_via_alias` how many more
are recovered by resolving a current symbol into this matrix's hg19
vintage (named in `*_alias_pairs_applied`), and the three `n_absent_*`
columns split what is left: `expression_filtered` for symbols the
count matrix carries but filterByExpr dropped, `undetected` for
symbols the CellRanger reference carries but the matrix does not,
`absent_from_reference` for a vocabulary miss outright. Those three
are reported separately because collapsing them reads a power fact and
a nomenclature fact as the same biological absence. On `frozen_sets`
rows the ledger columns are trivial by construction, since that
universe applies no restriction. `circle_radius_*`, `centre_distance`
and `shared_area_solved` are the geometry the figure drew, and
`shared_area_residual_genes` is how far the drawn shared area misses
the count in gene units — the proof behind `is_area_proportional`.
This is membership: the file carries no NES, p-value or effect size.
Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition_viz.py` | `euler_table` | `rows = 2 vocabularies x 3 regions; arm = WT_heat_up; lens = HALLMARK_HYPOXIA; vocabulary = symbol_alias.ranked_list at population=treg` | `03_results/13_arm_decomposition/tables/arm_program_summary.csv, 03_results/13_arm_decomposition/tables/arm_program_gene.csv, 00_data/references/msigdb_hallmark/HALLMARK_HYPOXIA.txt, 00_data/references/symbol_alias/symbol_alias_map.csv, 03_results/03_pseudobulk/tables/ranked_treg.tsv, 03_results/03_pseudobulk/tables/gene_symbols.csv, 03_results/00_build/tables/reference_feature_symbols.csv` |
