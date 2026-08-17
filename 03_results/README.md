# Results — JIA sorted synovial and blood T cells

GSE160097 is FACS-sorted Treg, Tcon and CD8 from the synovial fluid and the paired peripheral
blood of the same juvenile idiopathic arthritis patients: seven donors, six paired in each
population after QC. That pairing is what lets this compartment ask the niche question.

**Abbreviations.** SF = synovial fluid (inflamed joint). PB = peripheral blood. Treg =
CD4⁺CD127ˡᵒCD25⁺ regulatory. Tcon = CD4⁺CD25⁻ conventional. CD8 = CD8⁺CD45RO⁺ memory.
NES = normalised enrichment score.

## The question, and what the tree does with it

The mouse anchor hands over three frozen up arms in human symbols. This compartment asks whether
the mouse 39 °C-derived signature separates the inflamed synovial niche from paired blood within
a frozen cell state, and whether any part of that separation resists being reduced to the other
stresses the same niche imposes.

The arc runs in one direction. Build and QC the object, freeze the sort labels, aggregate to
donor-level pseudobulk and fit the paired contrast, then score the mouse arms against those
rankings. The rest of the tree presses on the answer: purging hypoxia genes, decomposing the arm
against curated lenses, running the whole unbiased sweep for calibration, and localising every
score on two embeddings.

## The two tracks, and how they rank

**The confirmatory spine is the only track that supports a claim.** Donor-level pseudobulk on raw
counts within frozen sort labels, limma-voom, then pre-ranked enrichment on the moderated t. Each
donor casts one vote. [`03_pseudobulk`](03_pseudobulk/), [`05_scoring`](05_scoring/),
[`09_heat_hypoxia`](09_heat_hypoxia/) and [`14_unbiased_enrichment`](14_unbiased_enrichment/)
carry it, and [`master`](master/) accumulates its effect sizes.

**Everything per-cell is annotation.** A per-cell score localises a program on a map. It pools
thousands of cells from donors of unequal yield, so a tissue difference read off a colour or a
violin is pseudoreplicated. Read those figures for where a score sits. Read the confirmatory ones
for whether it separates.

## The reading order

| Stage | Question it answers |
|---|---|
| [`00_build`](00_build/) | What the forty libraries contain, and which strata the design holds. |
| [`01_qc`](01_qc/) | Which cells survive, and what the high-mitochondrial pocket is. |
| [`02_annotation`](02_annotation/) | Whether the FACS sort labels track the transcriptome. |
| [`03_pseudobulk`](03_pseudobulk/) | The donor-paired synovial-versus-blood contrast, per population. |
| [`05_scoring`](05_scoring/) | Whether the mouse 39 °C arms enrich in that contrast. |
| [`07_embedding`](07_embedding/) | The population-of-interest harvest design, previewed. |
| [`08_harvest_readout`](08_harvest_readout/) | Per-cell hypoxia and proteostasis readouts on one scale. |
| [`09_heat_hypoxia`](09_heat_hypoxia/) | Whether the arm's enrichment reduces to its hypoxia gene content. |
| [`10_hsr_lens`](10_hsr_lens/) | What an anchor-independent proteostasis lens returns. |
| [`11_heat_decomposition`](11_heat_decomposition/) | Where each part of the arm sits in the same rankings. |
| [`12_treg_localisation`](12_treg_localisation/) | Per-cell score distributions inside the Treg gate. |
| [`13_arm_decomposition`](13_arm_decomposition/) | What the arms are made of, by membership against nine curated lenses. |
| [`14_unbiased_enrichment`](14_unbiased_enrichment/) | The whole sweep — what this contrast contains, with no set privileged. |
| [`15_coresh_search`](15_coresh_search/) | Where else in public human data the synovial up-arm co-varies. |
| [`16_narrative_scoring`](16_narrative_scoring/) | Every program localised on the full-object map. |
| [`17_treg_reembedding`](17_treg_reembedding/) | The same programs on a Treg-only, batch-corrected map. |
| [`18_tf_activity`](18_tf_activity/) | What an inferred HIF1A activity is, and what bounds it. |
| [`19_regulon_nulls`](19_regulon_nulls/) | Two nulls that hold more of the network structure fixed. |
| [`interactive`](interactive/) | Per-cell substrates for the explorer views. |
| [`master`](master/) | The confirmatory effect sizes, with intervals. |
| [`objects`](objects/) | Recomputable checkpoints. |

## What the compartment established

**The mouse 39 °C up arm separates the niche, and it does so in every sorted population.** NES
2.59 in Treg on 120 of 202 arm genes, 2.68 in Tcon on 131, 2.06 in CD8 on 114, every pooled FDR
below 1e-4 (`14_unbiased_enrichment/tables/gsea_all.csv`). The Treg score sits between the Tcon
and CD8 scores, so the separation is pan-T, with Treg one of the three populations carrying it.

**Its up and down arms move together.** The down arm reaches NES 1.432 at FDR 0.0354 on 64 genes
in Tcon, at the same sign as the up arm, and carries no direction in Treg or CD8
(`05_scoring/tables/gsea_pseudobulk_*.csv`). Both arms carry information, and the pattern is a
shared non-directional shift. Recapitulating the mouse contrast would take the two arms apart.

**By composition the arm is largely inflammatory.** Nine curated anchor-independent lenses
contain 67 of the 202 `WT_heat_up` genes, leaving 135 unclaimed as the largest single part. What
the lenses do claim runs 35 TNFα/NF-κB genes and 21 inflammatory-response genes against 2 in the
curated heat-shock core and 2 of the 21 published interferon-independent STING genes
(`13_arm_decomposition/tables/arm_program_summary.csv`).

**Its enrichment survives deleting its hypoxia gene content.** Removing the 18 HALLMARK_HYPOXIA
overlap genes costs 0.126 to 0.164 NES and leaves all three populations significant
(`09_heat_hypoxia/tables/gene_purge_nes_comparison.csv`). That is a statement about gene content
and nothing more.

**The niche moves an enormous amount at once.** 1,443 of 11,514 tests reach pooled FDR 0.05 in
Treg, 2,043 of 11,752 in Tcon, 939 of 11,532 in CD8, and the largest single effects run toward
blood — `WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS` at NES −3.45. Among the programs moving toward the
joint in sorted CD4 cells the arm sits near the top of the distribution. That calibration is what
[`14_unbiased_enrichment`](14_unbiased_enrichment/) supplies.

**Temperature and hypoxia are jointly imposed by the inflamed joint** and stay entangled in
cross-sectional human data. Nothing here separates them, and no artifact in this tree asserts
that one is a confound of the other.

## Signature provenance

Every gene set scored in this tree comes from outside the JIA data. This is the index: the
accession, how the list was derived, and the paper reference where one exists. Each stage README
repeats the entries it uses, so a caption stands on its own.

| Set | Origin | How the list was derived |
|---|---|---|
| `WT_heat_up/down`, `KO_heat_up/down`, `Interaction_up`, `Interaction_fdrOnly_up` | **GSE329522**, this project's own mouse anchor. No paper reference recorded. | Bulk RNA-seq of induced regulatory T cells differentiated from primary murine splenic CD4⁺ T cells, in a 2×2 design of genotype (WT, cGAS-KO) × temperature (37 °C, 39 °C), 5 biological replicates per group over 20 libraries. Projected to human orthologs with pinned offline babelgene and read from `../mouse_anchor/03_results/human_projection/signatures/`. |
| `sting_specific_up`, scored as `sting_specific_published` (21 genes) | SAVI PBMC **GSE226598**. de Cevins et al. 2023, *Cell Reports Medicine*, PMID 38118407, PMC10772457, Table S6 (supplement `mmc7.xlsx`, sheet "SAVI signature"). | The published interferon-independent STING-activation signature: the genes most specific to the SAVI disease-associated monocyte cluster after every type-I interferon transcript and every IFN-β-inducible gene is removed, so it is interferon-independent by construction. SAVI is monogenic and PBMC-derived, which makes this a positive-control reference: overlap with it is consistent with STING activation. |
| `ifn_only_up` / `ifn_only_down`, scored as `ifn_generic_axis` (200 up, 200 down) | **GSE226572**, an interferon-β time course from the same study family, de Cevins et al. 2023, PMID 38118407. The list itself is derived in this project. | `IFNb_vs_0h` donor-pseudobulk differential expression over 3 healthy donors, paired within donor, called at FDR below 0.05 and absolute log2 fold change at least 1.0. None of the 21 `sting_specific` genes appears in the up set. |
| `eTreg_up` / `eTreg_down` / `score_eTreg` | **GSE161426**, 26 bulk RNA-seq samples of sorted CD4 populations. Mijnheer / Lutter et al. 2021, *Nature Communications*, PMID 33976194, doi 10.1038/s41467-021-22975-7. | Derived here as a synovial-fluid Treg against peripheral-blood Treg contrast on the deposited log2-normalised matrix `GSE161426_Gene_expression_table_log2.xlsx` (32,584 genes × 26 samples). GEO carries a matrix alone for this series, so the list is computed in this compartment. |
| `HSR_core` (56) · `HSR_sensitivity` (176) | **MSigDB v2026.1.Hs**, retrieved offline through **msigdbr 26.1.0**. | Union of `REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS`, `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` and `GOBP_RESPONSE_TO_HEAT`, mapping to 176 genes. `HSR_core` is the taxonomy categories `hsf1_core_hsr` and `co_chaperone` inside that union. |
| The six `HALLMARK_*` sets | **MSigDB Hallmark collection H**, *Homo sapiens*, **v2026.1.Hs**, through **msigdbr 26.1.0**. | Frozen one symbol per line under `00_data/references/msigdb_hallmark/`, each with a validated expected size: HYPOXIA 200, UNFOLDED_PROTEIN_RESPONSE 113, TNFA_SIGNALING_VIA_NFKB 200, INTERFERON_ALPHA_RESPONSE 97, INFLAMMATORY_RESPONSE 200, IL2_STAT5_SIGNALING 199. |
| CollecTRI regulons | The CollecTRI transcription-factor-to-target collection. No paper reference recorded. | Human table read under a SHA-256 pin from `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, built locally because `decoupleR::get_collectri()` fails against OmnipathR 3.18.4. |
| PROGENy footprints | The PROGENy pathway-footprint model, human, 500 targets per pathway. No paper reference recorded. | Scored with decoupleR MLM. |
| CoReSh compendium | A public GEO compendium distributed through Synapse **syn66227307**. No paper reference recorded. | Consumed read-only from the shared reference cache, human (`hsa`) half. |

The two SAVI-derived axes are frozen in a separate compartment of this project and reached across
a relative path. Their accessions, derivations and citations above are what a reader needs, and
the path is an internal pointer.

## Layout

Every stage holds `figures/` and `tables/`, with `_overview/` for cross-population artifacts and
`by_contrast/<population>/` where a stage runs per population. A figure's source table is its
same-stem neighbour, so `figures/_overview/foo.png` draws from `tables/_overview/foo.csv` and the
viz script computes nothing.

`.png` and `.pdf` files are ignored repo-wide. The tables and these READMEs are the tracked
record, and every figure regenerates from the script its caption names.
