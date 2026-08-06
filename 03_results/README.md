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

## Layout

Every stage holds `figures/` and `tables/`, with `_overview/` for cross-population artifacts and
`by_contrast/<population>/` where a stage runs per population. A figure's source table is its
same-stem neighbour, so `figures/_overview/foo.png` draws from `tables/_overview/foo.csv` and the
viz script computes nothing.

`.png` and `.pdf` files are ignored repo-wide. The tables and these READMEs are the tracked
record, and every figure regenerates from the script its caption names.
