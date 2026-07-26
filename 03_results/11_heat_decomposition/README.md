# 11_heat_decomposition — artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = BH-adjusted p-value._

The mouse 39 °C up-arm enriches toward synovial fluid in every sorted population and survives a hypoxia-gene purge, but it is 199 genes doing many different things. Here I split the projected signature into subcomponents and ask where each one sits in the same donor-pseudobulk ranked lists. The whole stage is annotation tier: no row reaches `effect_sizes_treg_arthritis.csv` or any `03_results/master/` accumulator.

The parts are defined by intersection with curated, versioned, anchor-independent public gene sets — the frozen HSR core plus five MSigDB Hallmark programs. The `WT_heat_up` leading-edge taxonomy is deliberately unused: it covers only the 66 genes that are the union of the three populations' leading edges, so scoring subsets of it would score genes selected because they had already enriched.

Parts overlap. A gene in two curated programs appears in both, 25 of the 199 up genes are multiply claimed, and the NES rows therefore do not sum to the arm. Forcing a priority-ordered disjoint partition would silently decide which program gets credit for a shared gene, so the full per-gene membership is published instead. The genes no curated set claims are reported as their own `unassigned` part.

Several parts are small. A part whose intersection with a ranked list falls under `gsea_min_size` = 5 gets no NES and is reported as untestable with its size and its reason, on the face of the coverage figure as well as in the tables. Silent truncation would read as full coverage.

### Where the curated sets came from

The HSR core is anchor-independent and MSigDB-derived: 56 genes, the `hsf1_core_hsr` plus `co_chaperone` categories of a per-gene taxonomy over the 176-gene union of three human MSigDB v2026.1.Hs sets (`REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS` 101, `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` 82, `GOBP_RESPONSE_TO_HEAT` 104). It shares two genes with `WT_heat_up`, so intersecting the two measures something rather than restating it.

The five Hallmark programs come from the offline `msigdbr` 26.1.0 package, frozen with validated sizes by `02_analysis/scripts/freeze_hallmark_sets.R`, and are used WHOLE — no taxonomy refinement, where the HSR union was refined from 176 down to 56. That asymmetry deserves naming: for a purge or a claim test the unrefined set is the conservative choice, because a larger curated set claims more of the mouse signature and so understates what is left over.

| Presumption | Curated set | n | Frozen by |
|---|---|---|---|
| HSF1 thermal core | `HSR_core` (taxonomy `hsf1_core_hsr` + `co_chaperone`) | 56 | `02_analysis/scripts/freeze_hsr_lens.R` |
| unfolded-protein response | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 113 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| hypoxia | `HALLMARK_HYPOXIA` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| TNFA / NF-kB signalling | `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| type-I interferon | `HALLMARK_INTERFERON_ALPHA_RESPONSE` | 97 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| inflammatory response | `HALLMARK_INFLAMMATORY_RESPONSE` | 200 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| IL2-STAT5 activation | `HALLMARK_IL2_STAT5_SIGNALING` | 199 | `02_analysis/scripts/freeze_hallmark_sets.R` |
| no named program | the residual, claimed by none of the above | 137 up / 83 down | — |

The frozen lists themselves live under `00_data/references/`, which is not committed. The reproducer scripts are, so any clone regenerates byte-identical lists, and a size drift in the installed `msigdbr` is a hard stop rather than a silent shift.

### What the decomposition returns

Every testable up-arm part enriches toward synovial fluid, so the shift is broad rather than localised. The strongest part in Treg and CD8 is the 137-gene remainder that no curated program claims (+2.15 and +2.02). The TNFA/NF-kB part is the most CD4-selective (+2.06 Treg, +2.20 Tcon, +1.43 CD8), and the curated IL2-STAT5 activation proxy is the weakest in Treg (+1.12 at FDR 0.39).

Two nulls carry as much weight as the positives. The canonical HSF1 thermal core contributes 2 of the 199 up genes and type-I interferon contributes 1, both far under the size floor, so in this projection the mouse thermal program is neither a heat-shock-transcript program nor an interferon program. The down arm tells the same story: 83 of its 94 genes belong to no named program, and nothing in it separates synovial fluid from blood.

The cGAS/STING axis falls the same way. Of the published 21-gene interferon-independent STING signature, PLAUR and PTGS2 sit in the up arm and none in the down arm, and PLAUR is itself one of the 18 hypoxia-purged genes. Two of 21 is far below the size floor, so it stays a tally rather than an arm: the mouse thermal program is essentially not a STING program.

## tables/decomposition_overlap.csv

Curated public gene sets claim 62 of the 199 mouse up genes and 11 of the 94 down genes, with TNFA/NF-kB the largest single claim at 35 and the HSF1 thermal core at 2.

**How to read:** One row per mouse arm and presumption, plus an `unassigned` row per arm. `n_intersect` is how many of that arm's genes the curated set contains, `frac_of_mouse_arm` its share of the 199 or 94, and `frac_of_curated_set` how much of the public set the mouse arm covers. `genes` is the semicolon-delimited membership. Rows overlap because presumptions overlap, so `n_intersect` does not sum to the arm. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `overlap_tallies` | `signature=WT_heat; evidence_tier=secondary_annotation` | `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt` |

## tables/decomposition_gene_assignment.csv

25 of the 199 up genes and 2 of the 94 down genes belong to more than one curated presumption, with ATF3, CDKN1A, F3, PLAUR and SERPINE1 claimed by three, which is why the parts overlap instead of partitioning.

**How to read:** One row per mouse-signature gene. `subcomponents` is the semicolon-delimited list of every presumption claiming it and `n_subcomponents` the count; a gene claimed by none reads `unassigned`. This is the audit trail behind the overlapping decomposition — it shows exactly which genes two parts share, so a reader can judge how independent any two NES values are. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `gene_assignment` | `signature=WT_heat` | `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt`, `00_data/references/msigdb_hallmark/HALLMARK_*.txt`, `00_data/references/temp_hsr_lens/HSR_core.txt` |

## tables/decomposition_nes.csv

Every testable up-arm part enriches toward synovial fluid (NES +1.12 to +2.15) while no down-arm part reaches significance, and 27 of the 48 population-by-part cells are untestable and carry their reason.

**How to read:** One row per population, mouse arm and part — every requested part, whether or not fgsea could score it. `n_genes` is the part's size, `set_size_in_ranked` its intersection with that population's ranked list, and `testable` is False when that intersection falls under `gsea_min_size`, with `untestable_reason` naming which condition failed. Positive `nes` means enrichment toward SF-up genes and `padj` is BH across the parts scored within one population. Annotation tier, never pooled with the pseudobulk spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `decomposition_nes` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/11_heat_decomposition/tables/_signatures_decomp/*.txt` |

## tables/sting_axis_overlap.csv

The published interferon-independent STING signature contributes PLAUR and PTGS2 to the mouse up arm and nothing to the down arm, and PLAUR is itself a `HALLMARK_HYPOXIA` gene, so the mouse thermal program is essentially not a STING program.

**How to read:** One row per mouse arm. `n_intersect` and `genes_intersect` are the shared genes, `n_intersect_also_in_hypoxia` counts how many of them the hypoxia purge also removes, and `testable_as_gsea_arm` records that the overlap sits under `gsea_min_size`. That last column is the reason this is a gene tally and no STING NES appears anywhere in the stage. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition.py` | `sting_axis_overlap` | `gsea_min_size=5; evidence_tier=secondary_annotation` | `../sting_positive_control/03_results/06_reference_axis/signatures/sting_specific_up.txt`, `03_results/09_heat_hypoxia/tables/_signatures_full/WT_heat_{up,down}.txt` |

## tables/_signatures_decomp/{subcomponent}_{up,down}.txt

The sixteen sub-signature lists are the decomposition itself, and their sizes are the shape of the finding: 137 up genes unclaimed against 2 in the HSF1 thermal core, with three down-arm parts empty.

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

Curated public gene sets claim only 62 of the 199 mouse heat up genes
and 11 of the 94 down genes, so the largest part of the projected
signature — 137 up genes — belongs to no named program, and the
canonical HSF1 thermal core contributes 2 genes.

**How to read:** One bar per mouse arm and curated presumption; length is how many of
that arm's genes the curated set contains. Warm brown = the 199-gene
up arm, cool blue = the 94-gene down arm. The right-hand text gives
the count, then the testability: parts reaching 5 genes in the ranked
lists are tested, smaller parts are marked as under the floor, and a
part with no gene in that arm says so. Presumptions overlap, so bars
share genes and do not sum to the arm; per-gene memberships are in
decomposition_gene_assignment.csv. The published 21-gene interferon-
independent STING signature contributes only PLAUR and PTGS2 here,
tallied in sting_axis_overlap.csv. Annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_coverage` | `gsea_min_size=5; figures.top_n=20; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/decomposition_overlap.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv, 03_results/11_heat_decomposition/tables/sting_axis_overlap.csv` |

## figures/_overview/heatdecomp_runsum_up_unassigned.png

The 137 up-arm genes that no curated presumption claims give the
strongest synovial-fluid enrichment of any part in Treg (+2.15) and
CD8 (+2.02), so the shift is not carried by any single named program.

**How to read:** This part is the residual: the up-arm genes belonging to none of the
curated presumptions. Top panel: the weighted running enrichment score
as each population's ranked list is walked from synovial-fluid-up
(left) to blood-up (right); a positive, left-shifted excursion is
synovial-fluid enrichment and a negative trace the opposite. Bottom
panel: where this part's genes sit in each ranking, in matching
colour. Legend labels carry the testable gene count, the NES and the
FDR, and no other glyph marks significance. The y-range is shared
across every figure of this decomposition family, so curve heights
compare between figures. Annotation tier, firewalled from the
confirmatory WT_heat effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_unassigned_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_nfkb_tnfa.png

The 35 TNFA/NF-kB up-arm genes enrich toward synovial fluid strongly
in Treg (+2.06) and Tcon (+2.20) and only weakly in CD8 (+1.43, FDR
0.087), making the inflammatory-signalling part the most CD4-selective
of the decomposition.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_TNFA_SIGNALING_VIA_NFKB. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across every figure of this decomposition family, so curve
heights compare between figures. Annotation tier, firewalled from the
confirmatory WT_heat effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_nfkb_tnfa_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_hypoxia.png

The 18 hypoxia-overlap up-arm genes enrich toward synovial fluid in
all three populations (+1.82 to +2.07), confirming that the hypoxic
co-exposure is real even though removing these genes barely moves the
whole-signature NES.

**How to read:** This part is the up-arm genes that also sit in HALLMARK_HYPOXIA, the
same 18 the whole-signature purge removes. Top panel: the weighted
running enrichment score as each population's ranked list is walked
from synovial-fluid-up (left) to blood-up (right); a positive, left-
shifted excursion is synovial-fluid enrichment and a negative trace
the opposite. Bottom panel: where this part's genes sit in each
ranking, in matching colour. Legend labels carry the testable gene
count, the NES and the FDR, and no other glyph marks significance. The
y-range is shared across every figure of this decomposition family, so
curve heights compare between figures. Annotation tier, firewalled
from the confirmatory WT_heat effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_hypoxia_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_inflammatory.png

The 21 inflammatory-response up-arm genes track the whole up-arm
(+1.52 to +1.96), adding no separation of their own beyond the broad
synovial-fluid shift.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_INFLAMMATORY_RESPONSE. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across every figure of this decomposition family, so curve
heights compare between figures. Annotation tier, firewalled from the
confirmatory WT_heat effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_inflammatory_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_up_t_activation.png

The 14 IL2-STAT5 activation up-arm genes are the weakest testable part
in Treg (+1.12, FDR 0.39) while reaching +1.80 in Tcon, so a curated
T-cell activation program does not account for the Treg shift.

**How to read:** This part is the up-arm genes that also sit in
HALLMARK_IL2_STAT5_SIGNALING, read as a curated proxy for T-cell
activation. Top panel: the weighted running enrichment score as each
population's ranked list is walked from synovial-fluid-up (left) to
blood-up (right); a positive, left-shifted excursion is synovial-fluid
enrichment and a negative trace the opposite. Bottom panel: where this
part's genes sit in each ranking, in matching colour. Legend labels
carry the testable gene count, the NES and the FDR, and no other glyph
marks significance. The y-range is shared across every figure of this
decomposition family, so curve heights compare between figures.
Annotation tier, firewalled from the confirmatory WT_heat effect-size
spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_t_activation_up.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |

## figures/_overview/heatdecomp_runsum_down_unassigned.png

The 83 down-arm genes no presumption claims sit nowhere in particular
— NES +0.95 in Treg, +1.27 in Tcon and -1.04 in CD8, none of them
significant — so the mouse down arm does not separate synovial fluid
from blood in either direction.

**How to read:** This part is the residual of the mouse down arm: the genes belonging
to none of the curated presumptions. Top panel: the weighted running
enrichment score as each population's ranked list is walked from
synovial-fluid-up (left) to blood-up (right); a positive, left-shifted
excursion is synovial-fluid enrichment and a negative trace the
opposite. Bottom panel: where this part's genes sit in each ranking,
in matching colour. Legend labels carry the testable gene count, the
NES and the FDR, and no other glyph marks significance. The y-range is
shared across every figure of this decomposition family, so curve
heights compare between figures. Annotation tier, firewalled from the
confirmatory WT_heat effect-size spine.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/11_heat_decomposition_viz.py` | `plot_subcomponent_runsum` | `figures.running_sum_heights=[2.4, 0.7]; thresholds.gsea_fdr=0.05; gsea_min_size=5; evidence_tier=secondary_annotation` | `03_results/11_heat_decomposition/tables/runsum_interactive_decomp_gsea_{treg,tcon,cd8}_unassigned_down.csv, 03_results/11_heat_decomposition/tables/decomposition_nes.csv` |
