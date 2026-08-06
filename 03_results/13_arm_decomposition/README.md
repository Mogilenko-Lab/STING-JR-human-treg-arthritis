# 13_arm_decomposition — What the mouse up arms are made of

The mouse 39 °C anchor hands this compartment four up arms in human projection. Elsewhere those
arms are scored, purged and read against curated lenses. This stage asks the question underneath
all of them: **which curated programs contain the genes of each arm, and how much of each arm no
program names.**

The answer is arithmetic over committed gene lists. Nine curated, versioned, anchor-independent
lenses are intersected with each arm gene by gene. Two figures draw the result — a composition
bar per arm and an area-proportional Euler of the arm against curated hypoxia — and three tables
carry the membership at gene, program and arm level.

**Containment and enrichment are separate measurements.** A lens containing a gene is exact and
carries no p-value. A lens enriching in a ranked list is a statistical claim, answered on
donor-level pseudobulk elsewhere. Nothing here carries an NES, an FDR, a direction or an effect
size, and no row reaches `../master/`.

## The four arms and the nine lenses

Arm sizes are read from the anchor's frozen `manifest.csv` at run time, so a change to that
contract stops the stage.

| Arm | Mouse contrast | Gate | Genes |
|---|---|---|---|
| `WT_heat_up` | WT, 39 vs 37 °C | `fdr_logfc` | 202 |
| `KO_heat_up` | cGAS-KO, 39 vs 37 °C | `fdr_logfc` | 221 |
| `Interaction_up` | heat × genotype | `fdr_logfc` | 7 |
| `Interaction_up_fdrOnly` | the same contrast, relaxed gate | `fdr_only` | 19 |

Each lens is read verbatim from a frozen file and knows nothing of the mouse anchor, so an
intersection is an independent measurement. The first seven are byte-identical to the lenses
[`../11_heat_decomposition/`](../11_heat_decomposition/) uses, which lets the two stages join on
`program`.

| Program key | Curated set | n |
|---|---|---|
| `hypoxia` | `HALLMARK_HYPOXIA` | 200 |
| `nfkb_tnfa` | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 |
| `inflammatory` | `HALLMARK_INFLAMMATORY_RESPONSE` | 200 |
| `t_activation` | `HALLMARK_IL2_STAT5_SIGNALING` | 199 |
| `ifn_type_i` | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 |
| `upr_er` | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 |
| `hsr_curated` | `HSR_core`, Reactome/GO-derived | 56 |
| `sting_specific_published` | the published interferon-independent STING signature (de Cevins et al. 2023, Table S6) | 21 |
| `ifn_generic_axis` | a generic type-I interferon axis, the genes most induced by IFN-β at 24 h in healthy donors | 200 |

Both interferon-family lenses report gene content. A positive read on the published STING set is
consistent with STING pathway activity and stops short of proving it, and overlap with the
generic axis reflects a type-I interferon response of any origin.

## What the membership returns

**The two large arms have the same shape, and most of it is unclaimed.** Nine lenses contain 67
of the 202 `WT_heat_up` genes and 68 of the 221 `KO_heat_up` genes, leaving remainders of 135
and 153 — the largest single part of each arm (`tables/arm_program_summary.csv`).

**The claimed genes are largely inflammatory.** 35 `WT_heat_up` genes sit in
`HALLMARK_TNFA_SIGNALING_VIA_NFKB` and 21 in `HALLMARK_INFLAMMATORY_RESPONSE`, against 18 in
`HALLMARK_HYPOXIA` and 14 in `HALLMARK_IL2_STAT5_SIGNALING`. Three containments are small enough
to name outright: the curated HSR core holds 2 genes, the unfolded-protein response 1, and
Hallmark type-I interferon 1. The published STING signature contributes 2 to `WT_heat_up` and 1
to `KO_heat_up`, and the generic interferon axis 6 and 4.

**The two Interaction arms invert that shape.** In the 7-gene arm, Hallmark type-I interferon
contains 4 genes and the generic axis 2, while hypoxia, activation, proteostasis and the
published STING set contain none. The relaxed gate sharpens it: 10 of 19 in Hallmark type-I
interferon, 6 in the generic axis. By composition, the genes the interaction contrast selected
carry generic interferon annotation. At 7 and 19 genes these are thin sets, and whether they
separate synovial fluid from blood is a question for the donor-pseudobulk panels.

## Three properties to hold while reading

**The bands overlap, so they form no partition.** 67 of the 202 `WT_heat_up` genes are claimed at
all, and those 67 carry 100 memberships, with ATF3 and PLAUR claimed by four lenses each. Adding
the per-lens counts over-counts by 33 and shrinks the 135-gene remainder. The sharing is
published per gene in `tables/arm_program_multiplicity.csv`.

**Two accountings, and one of them sums.** `n_intersect` is the duplicated count.
`weight_fractional = 1/n_programs_for_gene` gives each gene one unit of mass split across the
lenses containing it, and summing it over an arm returns the arm size exactly. The figure draws
widths from the fractional accounting and prints the duplicated counts inside the bands.

**The arms share structure by construction.** The mouse contrasts are linearly dependent as
model coefficients. That algebra holds for the coefficients and stops at the thresholded lists:
`WT_heat_up` and `KO_heat_up` share 185 genes, `Interaction_up` shares 0 with either, and
`Interaction_up_fdrOnly` holds all 7 `Interaction_up` genes among its 19. Agreement between
those rows is expected structure.

## Why this stage counts 135 unassigned where the neighbouring stage counts 137

Both are the same arm measured against different lens panels.
[`../11_heat_decomposition/`](../11_heat_decomposition/) reads `WT_heat_up` against seven
lenses. This stage adds the two interferon-family lenses, and those claim exactly five further
genes — `IL1RN`, `IL2`, `OAS3`, `SAMD4A`, `TGM1`. Restricting `arm_program_gene.csv` to the
original seven programs reproduces the neighbouring count exactly, and the compute script
asserts that against the committed table on every run. A remainder is meaningful against the
panel it was measured on, so quote every remainder together with its panel.

## The vocabulary decides the count

The counts above are taken over the frozen lists in full, which is the right universe for a
composition question. An enrichment statistic lives in a different universe: a ranked list
carries only genes that were detected and survived `filterByExpr`, and this compartment's count
matrix is frozen to a CellRanger hg19 symbol vintage while the reference sets ship current
symbols.

Both losses are large enough to change a quoted number. Against the Treg ranked list,
`HALLMARK_HYPOXIA` matches 139 of its 200 genes by exact symbol and 4 more by alias
(`CAVIN1`→`PTRF`, `CAVIN3`→`PRKCDBP`, `ERO1A`→`ERO1L`, `NOCT`→`CCRN4L`), giving 143.
`WT_heat_up` matches 119 and one more by alias (`DYNLT2B`→`TCTEX1D2`), giving 120. The 18-gene
containment becomes 12 (`tables/_overview/arm_hypoxia_euler.csv`).

---

## Figures

### `figures/_overview/arm_program_composition.png`

**How much of each mouse up arm the curated lenses claim.**
One horizontal bar per arm, labelled with its gene count and the anchor's gate. Band width is
the fractional share, so a gene in k lenses gives 1/k to each and the widths total 1.0. The
number printed in a band is the duplicated count of that arm's genes the lens contains, so
widths and numbers measure different things. A band too narrow for a digit carries its count to
the right of the bar. Grey is the remainder.

Nine lenses claim 67 of the 202 `WT_heat_up` genes and 68 of the 221 `KO_heat_up` genes, and
inflammatory content dominates what they claim. The thin Interaction arms invert that shape into
type-I interferon. Read each band on its own: the 67 claimed genes carry 100 memberships, so the
printed counts exceed the claimed genes by 33.
*Source* `tables/_overview/arm_program_composition.csv` ·
`02_analysis/scripts/13_arm_decomposition_viz.py`.

### `figures/_overview/arm_hypoxia_euler.png`

**The same two gene lists, drawn in two vocabularies.**
Two area-proportional Euler panels sharing one area-per-gene scale and one bounding box. Orange
is the mouse 39 °C-derived up arm in human projection, blue is frozen `HALLMARK_HYPOXIA`, and
each region carries its gene count. Left is the frozen lists in full, 202 and 200 genes. Right
keeps only what the Treg donor-pseudobulk ranked list carries, 120 and 143 genes. Every area is
exact — the largest residual across both panels is 9.9e-14 genes.

Curated hypoxia accounts for a small minority of the arm in either vocabulary: 18 shared genes
of 202 on the left, 12 of 120 on the right. Both counts are correct within their own vocabulary,
and quoting either without its vocabulary is the misreading this figure prevents. The six shared
genes the right panel drops are ADM, ADORA2B, CCN1, EGFR, F3 and TGM2. Absence has three causes
and the source table splits them: dropped by the expression filter, undetected in sorted T
cells, or absent from the reference vocabulary outright.
*Source* `tables/_overview/arm_hypoxia_euler.csv` ·
`02_analysis/scripts/13_arm_decomposition_viz.py`.

---

## Tables

### `tables/arm_program_gene.csv` — the membership itself

One row per (arm, program, gene). A gene contained by k lenses appears k times, each row
carrying `n_programs_for_gene = k` and `weight_fractional = 1/k`. A gene no lens contains gets
exactly one row reading `program = unassigned`, so the remainder sits in the table as rows of its
own. `gate` repeats the anchor's gate for that arm.

The 202 `WT_heat_up` genes generate 232 rows — 100 memberships across nine lenses plus 135
unclaimed genes.

### `tables/arm_program_summary.csv` — the per-lens tally

One row per (arm, program), plus one `unassigned` row per arm. `n_intersect` is how many of that
arm's genes the lens contains and `frac_of_arm` its share of `n_arm`. `n_curated_set` is the
lens size, so a small `n_intersect` against a large lens is a statement about the arm.
`genes` is the semicolon-delimited membership.

Rows within an arm overlap, so summing `n_intersect` exceeds `n_arm` by the amount of that
sharing. Use `weight_fractional` in `arm_program_gene.csv` for an accounting that totals the arm
exactly.

### `tables/arm_program_multiplicity.csv` — the sharing, per gene

One row per (arm, gene): every gene of every arm, whether or not a lens contains it.
`n_programs` counts the lenses containing it and `programs` names them. A gene no lens contains
reads `n_programs = 0`.

26 of the 67 claimed `WT_heat_up` genes are contained by two or more lenses — 21 by two, 3 by
three, and ATF3 and PLAUR by four. Filter to `n_programs >= 2` to see exactly which genes two
bands of the figure share.

### `tables/_overview/arm_program_composition.csv`

The figure's same-stem source, one row per plotted band, ordered by arm then by `band_order`.
`n_intersect` and `frac_of_arm` are the duplicated accounting printed inside the bands.
`weight_fractional_sum` and `frac_of_arm_fractional` are the fractional accounting that sets the
widths and totals 1.0 per arm. The `arm_*` columns repeat that arm's counts on every row —
`arm_n_claimed`, `arm_n_unassigned`, `arm_n_claims_total`, `arm_n_excess_claims`,
`arm_max_lenses_per_gene` — so a caption can quote them from one file. `is_partition` reads
False on every row by construction.

One line shows where the two accountings diverge: 35 `WT_heat_up` genes in TNFα/NF-κB render as
an 11.9% band where the duplicated accounting gives 17.3%, because 20 of those 35 sit in another
lens too.

### `tables/_overview/arm_hypoxia_euler.csv`

Six rows, one per (`universe` × `region`). `region` is `arm_only`, `shared` or `lens_only`,
`n_genes` is the count the figure draws that area to, and `genes` names them.
`universe` is `frozen_sets` for the lists as curated and `treg_ranked_list` for the restriction,
with `vocabulary` recording the file the restriction was made against.

The `arm_*` and `lens_*` columns carry the three-way symbol ledger: `n_exact_match`,
`n_via_alias` with the pairs applied, and three `n_absent_*` columns splitting
`expression_filtered`, `undetected` and `absent_from_reference`. Those three stay separate
because collapsing them reads a power fact and a nomenclature fact as one biological absence.
`circle_radius_*`, `centre_distance` and `shared_area_residual_genes` are the geometry, and the
residual is the proof behind `is_area_proportional`.
