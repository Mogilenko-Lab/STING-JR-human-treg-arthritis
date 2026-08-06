# 11_heat_decomposition — Where each part of the arm sits in the same rankings

The mouse 39 °C up arm enriches toward synovial fluid in every population and survives a
hypoxia-gene purge. It is also 199 genes doing many different things. This stage splits the
projected signature into subcomponents and asks where each one sits in the same donor-pseudobulk
rankings.

**Every testable up-arm part enriches toward synovial fluid**, so the shift is broad rather than
localised. The 137-gene remainder that no curated program claims stays strongly enriched (+2.21
Treg, +2.27 Tcon, +2.10 CD8). The TNFα/NF-κB part is the most CD4-selective (+2.24, +2.32, +1.23),
and the curated IL2-STAT5 activation proxy is the weakest in Treg (+1.32 at FDR 0.22).

**Two nulls carry as much weight as the positives.** The curated heat-shock core contributes 2 of
the 199 up genes and type-I interferon contributes 1, both far under the size floor, so those
curated gene contents explain very little of `WT_heat_up`. The down arm tells the same story: 83
of its 94 genes belong to no named program, and nothing in it separates the two tissues.

The cGAS/STING tally falls the same way. Of the published 21-gene interferon-independent STING
signature, PLAUR and PTGS2 sit in the up arm and none in the down arm, and PLAUR is itself one of
the 18 hypoxia-purged genes. Two of 21 is far below the size floor, so it stays a tally.

The whole stage is annotation tier and writes no effect-size row.

## Where the parts come from

The parts are defined by intersection with curated, versioned, **anchor-independent** public gene
sets: the frozen curated HSR core plus six MSigDB Hallmark programs. The `WT_heat_up`
leading-edge taxonomy is deliberately unused — it covers only the 66 genes that are the union of
the three populations' leading edges, so scoring subsets of it would score genes selected because
they had already enriched.

| Presumption | Curated set | n |
|---|---|---|
| curated HSR core (Reactome/GO) | `HSR_core` | 56 |
| unfolded-protein response | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 |
| hypoxia | `HALLMARK_HYPOXIA` | 200 |
| TNFα / NF-κB signalling | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 |
| type-I interferon | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 |
| inflammatory response | `HALLMARK_INFLAMMATORY_RESPONSE` | 200 |
| IL2-STAT5 activation | `HALLMARK_IL2_STAT5_SIGNALING` | 199 |
| no named program | the residual | 137 up / 83 down |

The six Hallmark programs are used **whole**, with no taxonomy refinement, where the HSR union
was refined from 176 down to 56. That asymmetry is deliberate: for a purge or a claim test the
unrefined set is the conservative choice, because a larger curated set claims more of the mouse
signature and so understates what is left over.

The frozen lists live under `00_data/references/`, which is untracked. The reproducer scripts are
tracked, so any clone regenerates byte-identical lists, and a size drift in the installed
`msigdbr` is a hard stop.

## Three properties of this decomposition to hold

**The parts overlap, and the size of the overlap shows they are no partition.** 62 of the 199 up
genes are claimed by a curated set at all, and those 62 carry 92 claims between them, because 25
belong to two or three sets at once. So the bars and the NES rows do not sum to the arm: adding
the named parts double-counts 30 claims and shrinks the 137-gene remainder, which is the largest
single part. Forcing a priority-ordered disjoint partition would silently decide which program
gets credit for a shared gene, so the full per-gene membership is published instead.

**Small parts are reported with their size.** A part whose intersection with a ranked list falls
under `gsea_min_size = 5` gets no score and is reported as untestable with its size and its
reason, on the coverage figure's face as well as in the tables. Silent truncation would read as
full coverage.

**The genes no curated set claims are their own part.** `unassigned` is a first-class part with
its own score, and it is the largest.

---

## Figures

### `figures/_overview/heatdecomp_arm_coverage.png`

**How much of each arm the curated presumptions claim.**
One bar per mouse arm and curated presumption; length is how many of that arm's genes the curated
set contains. Warm brown gives the 202-gene up arm and cool blue the 96-gene down arm. The
right-hand text gives the count, then the testability: parts reaching 5 genes in the ranked lists
are tested, smaller parts are marked under the floor, and a part with no gene in that arm says so.

Curated sets claim 62 of the 202 up genes and 11 of the 96 down genes, so **the largest part of
the projected signature — 140 up genes — belongs to no named program**, and the curated heat-shock
core contributes 2.

**Read each bar on its own.** 25 of the claimed up-arm genes sit in up to three sets, so adding
the named bars double-counts 30 claims and shrinks the remainder. That count prints on the face,
per arm in `decomposition_assignment_multiplicity.csv` and per gene in
`decomposition_gene_assignment.csv`. The remainder is reported as a remainder: unnamed, and
supporting no mechanism. The published interferon-independent STING signature contributes PLAUR
and PTGS2 here, tallied in `sting_axis_overlap.csv`.

This is arithmetic over committed files. The panel carries no enrichment statistic on its face.
*Source* `tables/_overview/heatdecomp_arm_coverage.csv` ·
`02_analysis/scripts/11_heat_decomposition_viz.py`.

### `figures/_overview/heatdecomp_runsum_{up,down}_<part>.png` — six panels

**Where one part of the signature sits along each population's ranking.**
Two stacked panels sharing a fractional-rank x axis. Top, the weighted running enrichment score
as each population's ranked list is walked from synovial-up (left) to blood-up (right), so a
positive left-shifted excursion is synovial enrichment. Bottom, where that part's genes sit in
each ranking, in matching colour. Legend labels carry the testable gene count, the NES and the
FDR, and no other glyph marks significance. **The y range is shared across the whole
decomposition family**, so curve heights compare between these six figures.

| Panel | The part | What it returns |
|---|---|---|
| `up_unassigned` | The 140 up-arm genes no presumption claims | The strongest part in CD8 (+2.10) and strongly enriched in Treg (+2.21) and Tcon (+2.27) — the shift spreads across the arm. |
| `up_nfkb_tnfa` | The 35 TNFα/NF-κB up genes | Strong in Treg (+2.24) and Tcon (+2.32), weak in CD8 (+1.23, FDR 0.22) — the most CD4-selective part. |
| `up_hypoxia` | The 18 hypoxia-overlap up genes | Enriched in all three (+1.81 to +2.07), so the part carries a shift of its own. Whether the **whole** set reduces to it is the separate question the deletion panel in [`../09_heat_hypoxia/`](../09_heat_hypoxia/) answers. |
| `up_inflammatory` | The 21 inflammatory-response up genes | +1.48 to +2.11, tracking the whole up arm, so their separation is the broad synovial shift itself. |
| `up_t_activation` | The 14 IL2-STAT5 up genes | The weakest testable part in Treg (+1.32, FDR 0.22) while reaching +1.89 in Tcon, so the Treg shift rests on more than a curated activation program. |
| `down_unassigned` | The 85 down-arm genes no presumption claims | +0.97 Treg, +1.41 Tcon, −1.12 CD8, none significant — this remainder separates the tissues in neither direction. The whole 96-gene down arm does reach significance in Tcon, at the up arm's sign, and the whole-set panels carry that. |

**The parts overlap, so read each score on its own.** Adding them, or ranking them as shares of
the whole, would double-count genes.
*Source* `tables/_overview/heatdecomp_runsum_*.csv` and the `runsum_interactive_*` traces ·
`02_analysis/scripts/11_heat_decomposition_viz.py`.

---

## Tables

### `tables/decomposition_overlap.csv`

One row per mouse arm and presumption, plus an `unassigned` row per arm. `n_intersect` is how
many of that arm's genes the curated set contains, `frac_of_mouse_arm` its share of the 199 or 94,
and `frac_of_curated_set` how much of the public set the arm covers. `genes` is the
semicolon-delimited membership. Rows overlap, so `n_intersect` does not sum to the arm.

Curated sets claim 62 of 199 up and 11 of 94 down, with TNFα/NF-κB the largest single claim at 35
and the curated HSR core at 2.

### `tables/decomposition_gene_assignment.csv` · `tables/decomposition_assignment_multiplicity.csv`

The audit trail behind the overlap. The first holds one row per signature gene, with
`subcomponents` the semicolon-delimited list of every presumption claiming it and
`n_subcomponents` the count; a gene claimed by none reads `unassigned`. 25 of the 62 claimed up
genes and 2 of the 11 claimed down genes belong to more than one presumption, with ATF3, CDKN1A,
F3, PLAUR and SERPINE1 claimed by three.

The second aggregates that per arm. `n_claimed` counts genes claimed at least once and
`n_unassigned` the rest; `n_claimed_once` and `n_claimed_multiply` split the claimed genes, with
`max_subcomponents_per_gene` the worst case. `n_claims_total` is the sum over genes of how many
sets claim each, and `n_excess_claims` is that total minus the number of claimed genes — exactly
the amount by which summing the per-set bars over-counts. `is_partition` reads False whenever
`n_claimed_multiply` is non-zero, which is the whole point: the coverage figure states the
constraint as a measured count.

### `tables/decomposition_nes.csv`

One row per population, mouse arm and part — **every requested part, whether or not it could be
scored**. `n_genes` is the part's size, `set_size_in_ranked` its intersection with that
population's ranked list, and `testable` is False when that intersection falls under
`gsea_min_size`, with `untestable_reason` naming which condition failed. Positive `nes` means
enrichment toward synovial-up genes and `padj` is BH across the parts scored within one
population.

Every testable up-arm part enriches toward synovial fluid (NES +1.23 to +2.32) and no down-arm
part reaches significance; 27 of the 48 population-by-part cells are untestable and carry their
reason.

### `tables/sting_axis_overlap.csv`

One row per mouse arm. `n_intersect` and `genes_intersect` are the shared genes,
`n_intersect_also_in_hypoxia` counts how many the hypoxia purge also removes, and
`testable_as_gsea_arm` records that the overlap sits under the size floor. That last column is why
this is a gene tally and why no STING enrichment score appears in this stage.

The published interferon-independent STING signature contributes PLAUR and PTGS2 to the up arm
and nothing to the down arm, and PLAUR is itself a hypoxia gene.

### The remaining files

| File | What it holds |
|---|---|
| `tables/_signatures_decomp/<part>_{up,down}.txt` | The sixteen sub-signature lists — the decomposition itself. Plain sorted HGNC symbols, arm and presumption carried by the filename; an empty file means no gene of that arm belongs to that set. Their sizes are the shape of the finding: 137 up unclaimed against 2 in the curated HSR core, with three down-arm parts empty. |
| `tables/decomp_gsea_{treg,tcon,cd8}.csv` · `.rds` | Each population's thirteen non-empty sub-signatures scored in one run, so the parts within a population share a single BH correction and their FDRs compare directly. A part under the floor appears with NA statistics and its intersection size; read testability from `decomposition_nes.csv`. |
| `tables/runsum_interactive_decomp_gsea_<population>_<part>_{up,down}.csv` | One row per ranked gene per part, including parts too small to score, so a shape can be inspected even where no score is reported. |
| `tables/_overview/heatdecomp_arm_coverage.csv` | One row per plotted bar: `n_intersect` the bar length, `frac_of_mouse_arm` its share, `set_size_in_ranked_min` / `_max` the testable range across the three populations, and `n_populations_testable` of `n_populations`. |
| `tables/_overview/heatdecomp_runsum_{up,down}_<part>.csv` | One row per population per figure: the NES and FDR printed in the legend, the testable size, the part's full size, and the ranked-list length. |
| `tables/source_hash_manifest.csv` | The SHA-256 pin on the published STING signature this stage reads for the tally. |
