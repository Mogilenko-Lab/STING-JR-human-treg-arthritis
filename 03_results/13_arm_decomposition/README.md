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

## tables/source_hash_manifest.csv

All fourteen inputs are recorded by SHA-256, and the published de Cevins STING gene set matches
the hash this compartment already committed for it, so the STING containment counts here and
the STING tally in the `11_heat_decomposition` tables read the same 21 genes.

**How to read:** One row per input file. `source_label` names the role (`arm_*`, `lens_*`, or
the projection manifest), `source_path` is the path relative to the super-repo root, and
`sha256` is the digest read on this run. `pin_status` is the column that matters and it
distinguishes two different strengths of guarantee:
`verified_against_11_heat_decomposition` means the digest was checked against a hash already
committed in this compartment, so a changed source is a hard stop that halts the run; `recorded`
means the digest is a provenance record of this run only, which a reviewer can diff against a
later run but which nothing enforced here. Reading a `recorded` row as a pin would overstate
what it does.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition.py` | `write_source_hashes` | `pinned lens = sting_specific_published (verified against the committed savi_sting_specific_up hash); all other inputs recorded` | `03_results/11_heat_decomposition/tables/source_hash_manifest.csv`, `../mouse_anchor/03_results/human_projection/{manifest.csv,signatures/**}`, `00_data/references/{msigdb_hallmark,temp_hsr_lens}/*.txt`, `../sting_positive_control/03_results/06_reference_axis/signatures/{sting_specific_up.txt,ifn_only_up.txt}` |

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

Nine curated anchor-independent lenses contain 67 of the 199
WT_heat_up genes and 68 of the 218 KO_heat_up genes, leaving
remainders of 132 and 150 genes as the largest single part of each
large arm, while what the lenses do claim is dominated by inflammatory
gene content (35 TNFA/NF-kB and 21 inflammatory-response genes in
WT_heat_up) against 2 in the curated HSR core and 2 of the 21
published IFN-independent STING genes. The thin Interaction arms
invert that shape: Hallmark type-I interferon contains 4 of the 7
genes at the fdr_logfc gate and 10 of the 18 at the relaxed fdr_only
gate.

**How to read:** Each band counts how many of one arm's genes a frozen curated lens
contains. That is set arithmetic over committed gene lists, so no NES,
FDR, direction or effect size appears here or anywhere in these
tables. One row per arm, named by how it was derived and labelled with
its gene count and the mouse anchor's gate (fdr_logfc, with fdr_only
the relaxed Interaction variant frozen as Interaction_fdrOnly_up.txt).
Band width is the fractional share: a gene in k lenses gives 1/k to
each, so widths total 1.0. The number in a band is the duplicated
count of that arm's genes the lens contains, so numbers and widths
disagree. A band too narrow for a digit carries its count to the right
of the bar, so every lens keeps its number. Grey is the remainder, the
genes no lens contains, left unnamed on purpose. The lenses overlap,
so the bands are a containment tally: 67 of the 199 WT_heat_up genes
are contained by at least one lens, and those 67 carry 100
memberships, with up to 4 lenses on one gene, so the printed counts
exceed the claimed genes by 33. Per-gene multiplicity is in
arm_program_multiplicity.csv. The four rows share structure by
construction: the mouse contrasts are linearly dependent as model
coefficients (WT_heat = KO_heat + Interaction), and the two
Interaction rows are one contrast at two gates, so agreement between
rows is expected. That algebra holds for the coefficients and stops at
the thresholded lists: WT_heat_up and KO_heat_up share 182 genes,
Interaction_up shares none with either, and Interaction_up_fdrOnly
holds all 7 Interaction_up genes among its 18. Annotation tier,
firewalled from the donor-pseudobulk claim spine; no effect-size row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_arm_decomposition_viz.py` | `plot_composition` | `colors.okabe_ito + colors.diverging.up; label_floor=0.045; accounting=fractional (1/n_programs_for_gene); evidence_tier=secondary_annotation` | `03_results/13_arm_decomposition/tables/arm_program_summary.csv, 03_results/13_arm_decomposition/tables/arm_program_gene.csv, 03_results/13_arm_decomposition/tables/arm_program_multiplicity.csv` |
