# 11_heat_decomposition — artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = BH-adjusted p-value._

The mouse 39 °C up-arm enriches toward synovial fluid in every sorted population and survives a hypoxia-gene purge, but it is 199 genes doing many different things. Here I split the projected signature into subcomponents and ask where each one sits in the same donor-pseudobulk ranked lists. The whole stage is annotation tier: no row reaches `effect_sizes_treg_arthritis.csv` or any `03_results/master/` accumulator.

The parts are defined by intersection with curated, versioned, anchor-independent public gene sets — the frozen curated HSR core (Reactome/GO) plus six MSigDB Hallmark programs. The `WT_heat_up` leading-edge taxonomy is deliberately unused: it covers only the 66 genes that are the union of the three populations' leading edges, so scoring subsets of it would score genes selected because they had already enriched.

Parts overlap, and the size of the overlap is what shows they are not a partition. 62 of the 199 up genes are claimed by a curated set at all, and those 62 carry 92 claims between them because 25 of them belong to two or three sets at once. So the bars and the NES rows do not sum to the arm: adding the named parts double-counts 30 claims and shrinks the 137-gene remainder, which is the largest single part. Forcing a priority-ordered disjoint partition would silently decide which program gets credit for a shared gene, so the full per-gene membership is published instead, the per-arm multiplicity is tabulated in `decomposition_assignment_multiplicity.csv`, and both are printed on the coverage figure's face. The genes no curated set claims are reported as their own `unassigned` part.

Several parts are small. A part whose intersection with a ranked list falls under `gsea_min_size` = 5 gets no NES and is reported as untestable with its size and its reason, on the face of the coverage figure as well as in the tables. Silent truncation would read as full coverage.

### Where the curated sets came from

The curated HSR core (Reactome/GO) is anchor-independent and MSigDB-derived: 56 genes, the `hsf1_core_hsr` plus `co_chaperone` categories of a per-gene taxonomy over the 176-gene union of three human MSigDB v2026.1.Hs sets (`REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS` 101, `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` 82, `GOBP_RESPONSE_TO_HEAT` 104). It shares two genes with `WT_heat_up`, so intersecting the two measures something rather than restating it.

The six Hallmark programs come from the offline `msigdbr` 26.1.0 package, frozen with validated sizes by `02_analysis/scripts/freeze_hallmark_sets.R`, and are used WHOLE — no taxonomy refinement, where the HSR union was refined from 176 down to 56. That asymmetry deserves naming: for a purge or a claim test the unrefined set is the conservative choice, because a larger curated set claims more of the mouse signature and so understates what is left over.

| Presumption | Curated set | n | Frozen by |
|---|---|---|---|
| curated HSR core (Reactome/GO) | `HSR_core` (taxonomy `hsf1_core_hsr` + `co_chaperone`) | 56 | `02_analysis/scripts/freeze_hsr_lens.R` |
| unfolded-protein response | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| hypoxia | `HALLMARK_HYPOXIA` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| TNFA / NF-kB signalling | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| type-I interferon | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| inflammatory response | `HALLMARK_INFLAMMATORY_RESPONSE` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| IL2-STAT5 activation | `HALLMARK_IL2_STAT5_SIGNALING` | 199 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| no named program | the residual, claimed by none of the above | 137 up / 83 down | — |

The frozen lists themselves live under `00_data/references/`, which is not committed. The reproducer scripts are, so any clone regenerates byte-identical lists, and a size drift in the installed `msigdbr` is a hard stop rather than a silent shift.

### What the decomposition returns

Every testable up-arm part enriches toward synovial fluid, so the shift is broad rather than localised. The 137-gene remainder that no curated program claims remains strongly enriched (+2.21 Treg, +2.27 Tcon, +2.10 CD8). The TNFA/NF-kB part is the most CD4-selective (+2.24 Treg, +2.32 Tcon, +1.23 CD8), and the curated IL2-STAT5 activation proxy is the weakest in Treg (+1.32 at FDR 0.22).

Two nulls carry as much weight as the positives. The curated HSR core (Reactome/GO) contributes 2 of the 199 up genes and type-I interferon contributes 1, both far under the size floor, so these curated gene contents explain very little of `WT_heat_up`. The down arm tells the same story: 83 of its 94 genes belong to no named program, and nothing in it separates synovial fluid from blood.

The cGAS/STING tally falls the same way. Of the published 21-gene interferon-independent STING signature, PLAUR and PTGS2 sit in the up arm and none in the down arm, and PLAUR is itself one of the 18 hypoxia-purged genes. Two of 21 is far below the size floor, so it stays a tally rather than an arm and does not support reading `WT_heat_up` as STING-specific.

## tables/decomposition_overlap.csv

Curated public gene sets claim 62 of the 199 `WT_heat_up` genes and 11 of the 94 `WT_heat_down` genes, with TNFA/NF-kB the largest single claim at 35 and the curated HSR core (Reactome/GO) at 2.

**How to read:** One row per mouse arm and presumption, plus an `unassigned` row per arm. `n_intersect` is how many of that arm's genes the curated set contains, `frac_of_mouse_arm` its share of the 199 or 94, and `frac_of_curated_set` how much of the public set the mouse arm covers. `genes` is the semicolon-delimited membership. Rows overlap because presumptions overlap, so `n_intersect` does not sum to the arm. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `overlap_tallies` | `signature=WT_heat; evidence_tier=secondary_annotation` | `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt` |

## tables/decomposition_gene_assignment.csv

25 of the 62 claimed up genes and 2 of the 11 claimed down genes belong to more than one curated presumption, with ATF3, CDKN1A, F3, PLAUR and SERPINE1 claimed by three, which is why the parts overlap instead of partitioning.

**How to read:** One row per mouse-signature gene. `subcomponents` is the semicolon-delimited list of every presumption claiming it and `n_subcomponents` the count; a gene claimed by none reads `unassigned`. This is the audit trail behind the overlapping decomposition — it shows exactly which genes two parts share, so a reader can judge how independent any two NES values are. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `gene_assignment` | `signature=WT_heat` | `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt` |

## tables/decomposition_assignment_multiplicity.csv

Neither arm's curated assignment is a partition: on the up arm 62 of 199 genes are claimed and carry 92 claims, so 30 claims are duplicates of a gene already counted, and on the down arm 11 claimed genes carry 13 claims.

**How to read:** One row per mouse arm. `n_claimed` counts genes claimed by at least one curated set and `n_unassigned` the rest; `n_claimed_once` and `n_claimed_multiply` split the claimed genes by how many sets claim them, with `max_subcomponents_per_gene` giving the worst case. `n_claims_total` is the sum over genes of how many sets claim each, and `n_excess_claims` is that total minus the number of claimed genes — exactly the amount by which summing the per-set bars over-counts. `is_partition` is `False` whenever `n_claimed_multiply` is non-zero, which is the whole point of the table: it exists so the coverage figure can state the constraint as a measured count rather than as a caution, and so a reader who sums the bars can see how far wrong that goes. Annotation tier; a plain aggregation of `decomposition_gene_assignment.csv` with no statistic in it.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `assignment_multiplicity` | `signature=WT_heat` | `03_results/11_heat_decomposition/tables/decomposition_gene_assignment.csv` (in-memory from `gene_assignment`) |

## tables/decomposition_nes.csv

Every testable up-arm part enriches toward synovial fluid (NES +1.23 to +2.32) while no down-arm part reaches significance, and 27 of the 48 population-by-part cells are untestable and carry their reason.

**How to read:** One row per population, mouse arm and part — every requested part, whether or not fgsea could score it. `n_genes` is the part's size, `set_size_in_ranked` its intersection with that population's ranked list, and `testable` is False when that intersection falls under `gsea_min_size`, with `untestable_reason` naming which condition failed. Positive `nes` means enrichment toward SF-up genes and `padj` is BH across the parts scored within one population. Annotation tier, never pooled with the pseudobulk spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `decomposition_nes` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/11_heat_decomposition/tables/_signatures_decomp/*.txt` |

## tables/sting_axis_overlap.csv

The published interferon-independent STING signature contributes PLAUR and PTGS2 to the mouse 39 C-derived up arm and nothing to the down arm, and PLAUR is itself a `HALLMARK_HYPOXIA` gene, so the overlap does not support reading `WT_heat_up` as STING-specific.

**How to read:** One row per mouse arm. `n_intersect` and `genes_intersect` are the shared genes, `n_intersect_also_in_hypoxia` counts how many of them the hypoxia purge also removes, and `testable_as_gsea_arm` records that the overlap sits under `gsea_min_size`. That last column is the reason this is a gene tally and no STING NES appears anywhere in the stage. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `sting_axis_overlap` | `gsea_min_size=5; evidence_tier=secondary_annotation` | `../sting_positive_control/03_results/06_reference_axis/signatures/sting_specific_up.txt`, `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt` |

## tables/_signatures_decomp/{subcomponent}_{up,down}.txt

The sixteen sub-signature lists are the decomposition itself, and their sizes are the shape of the finding: 137 up genes unclaimed against 2 in the curated HSR core (Reactome/GO), with three down-arm parts empty.

**How to read:** Plain newline-delimited HGNC symbols, one per line, sorted; the mouse arm and the presumption are both carried by the filename. An empty file means no gene of that arm belongs to that curated set. These are inputs rather than results — the exact gene universe handed to fgsea, regenerated from the frozen references on every run. Diff `unassigned_up.txt` against the full arm to see what the presumptions collectively claim.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `prepare_signature_dirs` | `signature=WT_heat` | `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt` |

## tables/decomp_gsea_{treg,tcon,cd8}.csv

Each population's thirteen non-empty sub-signatures are scored in one fgsea run, so the parts within a population share a single BH correction and their FDRs are directly comparable.

**How to read:** One file per sorted population, one row per sub-signature named `<part>_<arm>`. Positive `nes` is enrichment toward SF-up genes, `set_size` is the part's intersection with that ranked list, and `core_enrichment` is the slash-separated leading edge. A part under the size floor appears with `NA` statistics and its intersection size only; read testability from `decomposition_nes.csv`, which spells out the reason. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `run_fgsea` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/11_heat_decomposition/tables/_signatures_decomp/*.txt` |

## tables/decomp_gsea_{treg,tcon,cd8}.rds

The clusterProfiler `gseaResult` objects preserve the same decomposition runs the CSVs summarize, so any part's running sum can be reconstructed exactly if a later display needs it.

**How to read:** One RDS per population, written by the fgsea helper. These are compute substrates rather than separate statistics; interpret their NES and FDR through the sibling `decomp_gsea_{treg,tcon,cd8}.csv` and the summarized `decomposition_nes.csv`. Positive enrichment still means SF-up. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/helpers/fgsea_prerank.R` | `(top-level)` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/11_heat_decomposition/tables/_signatures_decomp/*.txt` |

## tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_{subcomponent}_{up,down}.csv

The running-sum substrates place each part of the mouse signature along each population's ranked list, including the parts too small to score, so a shape can be inspected even where no NES is reported.

**How to read:** One row per ranked gene. `running_es` is the weighted enrichment trace, `hit` marks the part's genes and `leading_edge` those inside the enrichment peak. Positive, left-shifted peaks correspond to SF-up enrichment. A file exists for every non-empty part, including parts under the size floor, whose trace is a legitimate curve carrying no NES. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/helpers/fgsea_prerank.R` | `(top-level)` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/11_heat_decomposition/tables/_signatures_decomp/*.txt` |

## tables/_overview/heatdecomp_arm_coverage.csv

The plotted coverage table pairs each presumption's claim on a mouse arm with how large that claim stays once intersected with the ranked lists, which is where the small parts lose testability.

**How to read:** One row per plotted bar. `n_intersect` is the bar length, `frac_of_mouse_arm` its share of the arm, `set_size_in_ranked_min`/`_max` the range of testable sizes across the three populations, and `n_populations_testable` of `n_populations` how many cleared `gsea_min_size`. Rows are ordered up arm first, then by claim size. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `coverage_table` | `gsea_min_size=5; figures.top_n=20` | `03_results/11_heat_decomposition/tables/decomposition_overlap.csv`, `03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## tables/_overview/heatdecomp_runsum_{up,down}_{subcomponent}.csv

The annotated numbers behind each running-sum figure: one part's NES, FDR, testable size and ranked-list length in every population.

**How to read:** One row per population. `nes` and `padj` are the values printed in that figure's legend, `set_size_in_ranked` the genes actually walked, `n_genes` the part's full size, and `n_ranked_genes` the length of the ranked list, which is why traces end at slightly different x. The traces themselves are the `runsum_interactive_*` tables. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `subcomponent_traces` | `evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/decomposition_nes.csv`, `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_*.csv` |

## figures/_overview/heatdecomp_arm_coverage.png

Curated public gene sets claim only 62 of the 199 mouse 39 C-derived
up genes and 11 of the 94 down genes, so the largest part of the
projected signature — 137 up genes — belongs to no named program, and
the curated HSR core (Reactome/GO) contributes 2 genes. The bars are
not a partition: 25 of those 62 claimed up genes belong to two or
three curated sets at once, so the 62 carry 92 claims and summing the
named bars double-counts 30 of them.

**How to read:** ANSWERS what the projected set is made of, by membership over frozen
versioned gene lists — arithmetic over committed files, not an effect
estimate, and no NES on the face. One bar per mouse arm and curated
presumption; length is how many of that arm's genes the curated set
contains. Warm brown = the 199-gene up arm, cool blue = the 94-gene
down arm. The right-hand text gives the count, then the testability:
parts reaching 5 genes in the ranked lists are tested, smaller parts
are marked under the floor, and a part with no gene in that arm says
so. **Do not sum the bars.** The assignment is not a partition — 25 of
the 62 claimed up-arm genes sit in two or three sets, so adding the
named bars double-counts 30 claims and shrinks the 137-gene remainder,
the largest single part. That count is on the face, per arm in
decomposition_assignment_multiplicity.csv, and per gene in
decomposition_gene_assignment.csv. The remainder is reported as a
remainder: it is not named, and is evidence for no mechanism. The
published 21-gene interferon-independent STING signature contributes
only PLAUR and PTGS2 here, tallied in sting_axis_overlap.csv.
Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_coverage` | `gsea_min_size=5; figures.top_n=20; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/decomposition_overlap.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv, 03_results/11_heat_decomposition/tables/decomposition_assignment_multiplicity.csv, 03_results/11_heat_decomposition/tables/sting_axis_overlap.csv` |

## figures/_overview/heatdecomp_runsum_up_unassigned.png

The 137 up-arm genes that no curated presumption claims give the
strongest synovial-fluid enrichment of any part in CD8 (+2.10) and
remain strongly enriched in Treg (+2.21) and Tcon (+2.27), so the
shift is not carried by any single named program.

**How to read:** This part is the residual: the up-arm genes belonging to none of the
curated presumptions. Top panel: the weighted running enrichment score
as each population's ranked list is walked from synovial-fluid-up
(left) to blood-up (right); a positive, left-shifted excursion is
synovial-fluid enrichment and a negative trace the opposite. Bottom
panel: where this part's genes sit in each ranking, in matching
colour. Legend labels carry the testable gene count, the NES and the
FDR, and no other glyph marks significance. The y-range is shared
across the whole decomposition family, so curve heights compare
between figures. CORROBORATES; does not answer the niche question —
annotation tier, firewalled from the confirmatory WT_heat spine, no
effect-size row. The parts overlap, so their NES may not be added or
ranked as shares of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_unassigned_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_nfkb_tnfa.png

The 35 TNFA/NF-kB up-arm genes enrich toward synovial fluid strongly
in Treg (+2.24) and Tcon (+2.32) and only weakly in CD8 (+1.23, FDR
0.22), making the inflammatory-signalling part the most CD4-selective
of the decomposition.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_TNFA_SIGNALING_VIA_NFKB. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across the whole decomposition family, so curve heights compare
between figures. CORROBORATES; does not answer the niche question —
annotation tier, firewalled from the confirmatory WT_heat spine, no
effect-size row. The parts overlap, so their NES may not be added or
ranked as shares of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_nfkb_tnfa_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_hypoxia.png

The 18 hypoxia-overlap up-arm genes enrich toward synovial fluid in
all three populations (+1.81 to +2.07), so this part carries a shift
of its own — which is a separate question from whether the whole set's
enrichment is reducible to it, and that one is answered by the
deletion panel rather than here.

**How to read:** This part is the up-arm genes that also sit in HALLMARK_HYPOXIA, the
same 18 the whole-signature purge removes. Top panel: the weighted
running enrichment score as each population's ranked list is walked
from synovial-fluid-up (left) to blood-up (right); a positive, left-
shifted excursion is synovial-fluid enrichment and a negative trace
the opposite. Bottom panel: where this part's genes sit in each
ranking, in matching colour. Legend labels carry the testable gene
count, the NES and the FDR, and no other glyph marks significance. The
y-range is shared across the whole decomposition family, so curve
heights compare between figures. CORROBORATES; does not answer the
niche question — annotation tier, firewalled from the confirmatory
WT_heat spine, no effect-size row. The parts overlap, so their NES may
not be added or ranked as shares of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_hypoxia_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_inflammatory.png

The 21 inflammatory-response up-arm genes track the whole up-arm
(+1.48 to +2.11), adding no separation of their own beyond the broad
synovial-fluid shift.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_INFLAMMATORY_RESPONSE. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across the whole decomposition family, so curve heights compare
between figures. CORROBORATES; does not answer the niche question —
annotation tier, firewalled from the confirmatory WT_heat spine, no
effect-size row. The parts overlap, so their NES may not be added or
ranked as shares of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_inflammatory_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_t_activation.png

The 14 IL2-STAT5 activation up-arm genes are the weakest testable part
in Treg (+1.32, FDR 0.22) while reaching +1.89 in Tcon, so a curated
T-cell activation program does not account for the Treg shift.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_IL2_STAT5_SIGNALING, read as a curated proxy for T-cell
activation. Top panel: the weighted running enrichment score as each
population's ranked list is walked from synovial-fluid-up (left) to
blood-up (right); a positive, left-shifted excursion is synovial-fluid
enrichment and a negative trace the opposite. Bottom panel: where this
part's genes sit in each ranking, in matching colour. Legend labels
carry the testable gene count, the NES and the FDR, and no other glyph
marks significance. The y-range is shared across the whole
decomposition family, so curve heights compare between figures.
CORROBORATES; does not answer the niche question — annotation tier,
firewalled from the confirmatory WT_heat spine, no effect-size row.
The parts overlap, so their NES may not be added or ranked as shares
of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_t_activation_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_down_unassigned.png

The 83 down-arm genes no presumption claims sit nowhere in particular
— NES +0.97 in Treg, +1.41 in Tcon and -1.12 in CD8, none of them
significant — so this remainder does not separate synovial fluid from
blood in either direction. Read that as a statement about the
remainder and not about the arm — the whole 94-gene down arm does
reach significance in Tcon, at the same sign as the up arm, and that
result belongs to the whole-set panels rather than to this one.

**How to read:** This part is the residual of the mouse down arm: the genes belonging
to none of the curated presumptions. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across the whole decomposition family, so curve heights compare
between figures. CORROBORATES; does not answer the niche question —
annotation tier, firewalled from the confirmatory WT_heat spine, no
effect-size row. The parts overlap, so their NES may not be added or
ranked as shares of the whole.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_unassigned_down.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |
