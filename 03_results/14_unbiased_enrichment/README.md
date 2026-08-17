# 14_unbiased_enrichment — The whole sweep, and where the mouse arm sits in it

Every targeted test in this compartment hands the JIA synovial-fluid-versus-paired-blood ranking
one named signature and asks whether it enriches. It does. A targeted test answers for the set it
was handed and stays silent on whether that answer is remarkable. This stage asks the unbiased
counterpart of the same three frozen rankings: **what does this niche contrast contain across
curated databases, with no set privileged?**

The two questions are informative together. If a hundred inflammatory programs move as strongly
as the mouse-derived arm, the arm's enrichment describes the niche. If almost nothing else moves
that strongly, the arm is distinctive. This stage says which.

**What was computed.** Pre-ranked fgsea over 15 gene-set collections in each of the three
sorted populations — 11,514 tests in Treg, 11,752 in Tcon, 11,532 in CD8 — on the same
donor-level pseudobulk rankings the confirmatory answer rests on, with two multiplicity
corrections published side by side. Then decoupleR MLM on the human PROGENy model, run both on
the contrast statistics and paired across donors. PROGENy needs no gene-set list, so it inherits
no curation decision.

**What was drawn.** Sixteen overview panels carry the reading — the calibration scatter, the top
sets, the two PROGENy panels, the arm and program NES-by-cell-state dots, the named sets placed
against the whole sweep, and nine running sums. Beneath them sits a browse surface of 141 panels
over every population × collection cell.

**Tier.** Enrichment statistics, correlative throughout. A set enriching means its gene content
moves with the synovial-fluid side of this ranking. No row here becomes an effect size or
reaches [`../master/`](../master/).

## What the sweep returns

**The niche contrast moves an enormous amount at once.** 1,443 of the 11,514 Treg tests reach
pooled FDR < 0.05 — 1,379 toward synovial fluid and 64 toward blood. Tcon is larger (2,043 of
11,752) and CD8 smaller (939 of 11,532) (`tables/gsea_all.csv`). This is what an inflamed tissue
niche against paired circulating cells looks like, and it is the context every targeted result
in this compartment reads against.

**The mouse-derived up arm is strong, and its rank depends on the population.** `WT_heat_up`
reaches NES +2.59 at pooled FDR 3.7e-12 in Treg on 120 testable genes. Among the 1,379
synovial-side significant sets it ranks **5th**, behind `eTreg_up` (+3.72),
`REACTOME_INTERLEUKIN_10_SIGNALING` (+2.64), `GOCC_COPII_COATED_ER_TO_GOLGI_TRANSPORT_VESICLE`
(+2.61) and `KO_heat_up` (+2.60). In Tcon it ranks 3rd of 1,908 at +2.68. In CD8 it is mid-pack,
73rd of 799 at +2.06, behind MHC class II, antigen presentation and interferon sets.

So the arm sits at the top of the synovial-side distribution in sorted CD4 cells and inside it
in CD8. A great deal co-enriches with it in every population. That is compatible with what this
compartment established by other routes: the synovial enrichment is shared across sorted
populations, and by curated composition the up arm is largely inflammatory. The sweep adds the
part a targeted test omits — what the arm is at the top of.

**`KO_heat_up` tracks it row for row** (+2.60 Treg, +2.66 Tcon, +2.06 CD8), which is the
arithmetic to expect from arms sharing 185 genes, and a further reason to read the WT arm as
carrying no cGAS-specific content. **`Interaction_up` reaches nothing anywhere** (+1.47 / +1.41
/ +1.54, pooled FDR 0.24 / 0.27 / 0.18): at 6 testable genes its null is a statement about size.

**The largest single shifts in this contrast are downward.** The top Treg results by pooled FDR
are translation and ribosome sets moving toward blood — `WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS` at
−3.45, `REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION` at −3.39, `KEGG_RIBOSOME` at −3.35. Only 64
of the 1,443 pooled-significant Treg sets go that way, so this is a small number of very large
effects. A targeted test of up arms sees none of it, which is the clearest argument for running
the sweep.

**Hypoxia and interferon both move, and the four oxygen-named sets disagree.**
`HALLMARK_HYPOXIA` reaches +2.27 at pooled FDR 1.4e-08 in Treg on 143 testable genes, while
`GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS` reaches +1.20 and clears nothing
(`tables/_overview/program_nes_by_cell_state.csv`). The generic type-I interferon axis enriches
robustly (+2.31 Treg, +2.33 Tcon, +2.60 CD8), while the published 21-gene interferon-independent
STING signature reaches pooled FDR 0.104 in Treg on 13 genes and clears the threshold in Tcon
alone (+1.95, 0.009, 16 genes). An interferon-like reading of this contrast survives. A
STING-specific reading stays unsupported.

**The curated proteostasis core changes sign between populations.** `HSR_core` reads +1.49 in
Treg, −1.33 in Tcon and −1.13 in CD8, clearing FDR 0.05 in none. That is a Treg-versus-others
sign difference at trend level.

## Two corrections, both published

`padj` is Benjamini-Hochberg within one collection, the correction a single-collection run
reports. `padj_pooled` corrects across every test asked of that population's ranking, the honest
correction for a sweep that interrogates one ranking fifteen times.

Pooling redistributes as much as it tightens. In Treg the total barely moves while the small
hand-picked collections lose and GO_BP gains, because Benjamini-Hochberg divides by rank as well
as multiplying by family size (`tables/gsea_pooled_summary_by_db.csv`). The effect that matters
is the first: a deliberately chosen collection of fifty stops benefiting from having been tested
alone. That correction is available only to a sweep.

**One set, one hypothesis.** Six MSigDB Hallmark sets are held twice, once live and once as this
compartment's frozen re-pin. Both copies are scored, their gene content is verified identical
set by set, and only the canonical MSigDB copy enters the pooled family. The agreement between
copies doubles as a drift check on the frozen files, and the whole record sits in
`tables/geneset_alias_map.csv`.

## Two gates the sweep runs before it starts

Count matrices here are keyed by Ensembl id while every reference set matches on HGNC symbol. A
leaked Ensembl-keyed ranking intersects all fifteen collections at approximately zero, and fgsea
returns empty rows with no error raised, so the failure reads as a biological null.
`tables/ranked_list_keycheck.csv` measures the key vocabulary and stops the run above a 0.5
Ensembl-like fraction, and `tables/geneset_overlap.csv` publishes the per-collection gene
overlap as positive evidence that the join joined.

The second gate is comparability. `WT_heat_up` runs as an ordinary member of the sweep, and its
NES must land on the published value off the same ranked list. It does, to within fgsea's own
stochastic normalisation and at identical effective set sizes, in all three populations
(`tables/wt_heat_up_reproduction.csv`). The script stops before touching the large collections
if that fails, and it was never tuned to agree.

## Signature provenance

The fifteen collections and the PROGENy model, with the accession or resource version each was
retrieved at. `tables/geneset_manifest.csv` carries the same record per collection.

| Collection | Origin | How the sets were obtained |
|---|---|---|
| Hallmark, KEGG, Reactome, WikiPathways, GO_BP, GO_MF, GO_CC | **MSigDB 2026.1.Hs**, *Homo sapiens*, retrieved offline through **msigdbr 26.1.0**. | Collections `H`, `C2/CP:KEGG_LEGACY` (the release serves KEGG under the legacy name), `C2/CP:REACTOME`, `C2/CP:WIKIPATHWAYS`, `C5/GO:BP`, `C5/GO:MF`, `C5/GO:CC`. |
| `project_frozen` | The same MSigDB release, re-pinned in this compartment. | The six Hallmark sets frozen under `00_data/references/msigdb_hallmark/` plus `HSR_core`, held as a drift check against the live copies. |
| `HSR_lens` | **MSigDB v2026.1.Hs** through **msigdbr 26.1.0**. | `HSR_sensitivity` is the union of `REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS`, `REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE` and `GOBP_RESPONSE_TO_HEAT` (176 genes). `HSR_core` (56) is its `hsf1_core_hsr` and `co_chaperone` taxonomy categories. |
| `TCR_activation` | A frozen literature-grounded human T-cell activation panel of 66 symbols. No paper reference recorded. | Curated as a human CSV spanning TCR-proximal signalling, early costimulation, immediate-early transcription factors and activation effector genes. |
| `mouse_projection` | **GSE329522**, this project's own mouse anchor. No paper reference recorded. | Bulk RNA-seq of induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype (WT, cGAS-KO) × temperature (37 °C, 39 °C), 5 biological replicates per group over 20 libraries. The thresholded up arms, projected to human orthologs with pinned offline babelgene. |
| `sting_axes` | `sting_specific_up`: SAVI PBMC **GSE226598**, de Cevins et al. 2023, *Cell Reports Medicine*, PMID 38118407, PMC10772457, Table S6 (`mmc7.xlsx`, sheet "SAVI signature"). `ifn_only_up`: **GSE226572**, an interferon-β time course from the same study family. | The published 21-gene set is the genes most specific to the SAVI disease-associated monocyte cluster after every type-I interferon transcript and every IFN-β-inducible gene is removed. The 200-gene generic axis is derived in this project from `IFNb_vs_0h` donor pseudobulk over 3 healthy donors, paired within donor, at FDR below 0.05 and absolute log2 fold change at least 1.0. |
| `eTreg_lens` | **GSE161426**, 26 bulk RNA-seq samples of sorted CD4 populations. Mijnheer / Lutter et al. 2021, *Nature Communications*, PMID 33976194, doi 10.1038/s41467-021-22975-7. | Derived here as a synovial-fluid Treg against peripheral-blood Treg contrast on the deposited log2-normalised matrix `GSE161426_Gene_expression_table_log2.xlsx`. That is what makes it another cohort's synovial-versus-blood contrast. |
| `MitoPathways` | MitoCarta 3.0, human build. No paper reference recorded. | Read from the RNAseq-toolkit reference build under `01_modules/`. |
| `TF_Targets` | CollecTRI transcription-factor-to-target regulons. No paper reference recorded. | The human table read under a SHA-256 pin from `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, each regulon flattened to one unsigned gene set. |
| PROGENy | The PROGENy pathway-footprint model, human, 500 targets per pathway. No paper reference recorded. | Scored with decoupleR MLM, so no gene-set list enters it. |

The two SAVI-derived axes are frozen in a separate compartment of this project and read across a
relative path. The accessions, derivations and citation above are what their rows in the sweep
rest on.

---

## Figures

### `figures/_overview/arm_nes_by_cell_state.png`

**The mouse up arms scored in each sorted cell state, at the confirmatory tier.**
One dot per arm and cell state. Rows are the three arms, and inside a row the three cell states
are offset vertically and coloured. x, normalised enrichment score for synovial fluid over
paired blood. A filled dot clears pooled FDR 0.05 and an open dot sits above it. The annotation
column gives how many of the arm's genes reached that population's ranked list, against how many
the frozen arm holds, then the adjusted p.

`WT_heat_up` rises toward synovial fluid in all three states — NES 2.5918 in Treg on 120 of 202
genes, 2.6823 in Tcon on 131, 2.0625 in CD8 on 114, every one below pooled FDR 1e-4. The Treg
score sits between the Tcon and CD8 scores, so the separation reads as pan-T. `KO_heat_up`
tracks it row for row. The 7-gene interaction arm reaches 1.4050 to 1.5435 on 6 testable genes
and stays above pooled FDR 0.05 in every state. Read every score with its gene count: resolution
scales with it, and these arms span 6 to 147 testable genes.
*Source* `tables/_overview/arm_nes_by_cell_state.csv` ·
`02_analysis/scripts/14_unbiased_enrichment_viz.py`.

### `figures/_overview/program_nes_by_cell_state.png`

**Every set the per-cell maps colour by, on one donor-level geometry.**
The same dot geometry as the arm panel, with rows grouped into three blocks: oxygen and HIF
response, cGAS-STING and interferon, then proteostasis, inflammation and the effector-Treg
reference. x, NES for synovial fluid over paired blood. Fill marks pooled FDR < 0.05.

The four hypoxia-named sets disagree with each other, spanning `HALLMARK_HYPOXIA` at pooled FDR
1.4e-08 down to `GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS` at 0.362 in Treg, so a hypoxia reading
depends on which curated set is asked. The curated proteostasis core changes sign between the
sorted populations (+1.49 Treg, −1.33 Tcon, −1.13 CD8) while clearing FDR in none of them. The
generic interferon axis clears comfortably in all three, and the published STING signature clears
in Tcon alone.
*Source* `tables/_overview/program_nes_by_cell_state.csv` ·
`02_analysis/scripts/14_program_nes_by_cell_state_viz.py`.

### `figures/_overview/progeny_paired_forest.png`

**Signalling footprints, one donor at a time.**
One row per PROGENy pathway, ordered by the Treg paired difference. Inside a row the three
populations are offset vertically — Treg above, Tcon centre, CD8 below — with colour repeating
the key. x, the mean within-donor difference in footprint activity, synovial fluid minus that
donor's own blood, and the bar through it is the 95% interval of the paired t-test. A filled
marker reaches FDR < 0.05.

Ten of the fourteen footprints separate synovial Treg from paired blood over six paired donors.
EGFR is the largest at +10.9 (FDR 2.7e-04), followed by JAK-STAT +9.0, NFκB +8.3, TGFβ +7.9 and
Hypoxia +7.4 (FDR 1.2e-03), with WNT at −3.5 and TNFα at −2.0 on the blood side
(`tables/progeny_sf_vs_pb.csv`). So the joint carries a measurable hypoxia footprint alongside
the inflammatory ones.

Two limits bind the reading. Activity is computed on expression centred within a population, so
a difference compares pathways inside one population and carries no meaning across populations.
Six pairs is also a small n, so an open marker leaves the question open. A footprint is inferred
from target-gene expression, and pathway activity itself is untested.
*Source* `tables/_overview/progeny_paired_forest.csv` ·
`02_analysis/scripts/14_progeny_paired_forest_viz.R`.

### `figures/_overview/progeny_activity_panel.png`

**The same fourteen footprints read off the contrast statistic.**
One row per pathway, ordered by its Treg score, with a grey line joining the three populations
so the spread within a row is the between-population difference. x, the MLM activity score on
the donor-pseudobulk moderated-t statistics. Colour is population. A solid point reaches FDR <
0.05 on this test, and a black ring marks one that also reaches it in the independent
donor-paired test, so ringed and solid is corroborated twice.

In Treg the largest footprint is JAK-STAT (+9.40, FDR 8.6e-20), then EGFR, TGFβ and Hypoxia
(+5.40, 2.4e-07), with WNT the only negative. Eight of the fourteen are significant on both
tests. The pattern repeats in all three populations, so it is a property of the niche.

The inflammatory and low-oxygen readouts rise together in one cross-sectional contrast, and that
is what these data show. Which of the two precedes the other stays untested here.
*Source* `tables/_overview/progeny_activity_panel.csv` ·
`02_analysis/scripts/14_unbiased_enrichment_viz.R`.

### `figures/_overview/named_sets_in_sweep.png`

**Seventeen named sets placed against every set tested.**
Columns are the three populations on one shared row axis. Each of the 17 upper rows is one named
set, coloured by comparison thread. Below the dashed separator the bottom row is every set
tested in that population, one grey point each, which is the distribution a marker is read
against. x, NES clamped to ±3.5. A filled marker reaches pooled FDR < 0.05, and marker area is
genes reaching the ranked list. Grey text gives pooled FDR and rank within that population's
whole sweep.

`HALLMARK_HYPOXIA` ranks 74 of 11,514 by pooled FDR in Treg. The best-placed of the six
cGAS-STING sets ranks 2,021, and `eTreg_up` — another cohort's synovial-versus-blood contrast —
ranks 1, which marks what a set built to separate exactly these two tissues reaches here.

Set size tracks the outcome closely: 43.6% of Treg sets of 130 to 150 genes reach pooled
significance against 6.1% of sets of 10 to 22 genes, the band five of the six cGAS-STING sets
fall in. The selection of rows, the reason for each, and the two excluded substring matches are
committed in `tables/sweep_named_sets.csv`.
*Source* `tables/_overview/named_sets_in_sweep.csv` ·
`02_analysis/scripts/14_sweep_named_sets_viz.R`.

### `figures/_overview/pooled_overview_by_population.png`

**How much of each collection moves, and where the mouse arms fall inside it.**
Columns are the three populations. Below the dashed line each row is a collection, ordered by
how many of its sets reach significance. Above it each row is one mouse-derived up arm. A small
point is a set at pooled FDR < 0.05, warm brown concentrating on the synovial-fluid side and
blue on paired blood. x, exact NES clamped to ±3.5. Yellow diamonds are the arms, filled when
significant, and grey text gives each collection's count or each arm's NES and FDR.

Read it for calibration. Far right in a dense row is strong and ordinary at once.
*Source* `tables/_overview/pooled_overview_by_population.csv` ·
`02_analysis/scripts/14_unbiased_enrichment_viz.R`.

### `figures/_overview/treg_top_sets.png`

**The twenty largest movers in the Treg contrast.**
One row per gene set, the top ten in each direction by absolute NES among pooled-significant
sets, with the subtitle stating how many the cap leaves out. Identifiers are wrapped in full,
each with its collection in brackets. Right of zero the set's genes concentrate on the
synovial-fluid side. Point size is genes reaching the ranked list, so a large NES on a small
point rests on few genes. The grey number is pooled FDR, and a black ring marks a mouse-derived
arm.

`WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS` reaches −3.45 toward blood against `eTreg_up` at +3.72
toward synovial fluid, so the inflammatory gain this compartment reads comes alongside a loss of
translation and ribosomal programs at least as large.
*Source* `tables/_overview/treg_top_sets.csv` ·
`02_analysis/scripts/14_unbiased_enrichment_viz.R`.

### `figures/_overview/runsum_<set>.png` — nine panels

**Where one named set sits along all three rankings.**
Three stacked panels sharing a fractional-rank x axis. Top, the weighted running enrichment
score as each population's list is walked from synovial-up to blood-up. Middle, where that set's
genes sit in each ranking. Bottom, the ranked moderated t each curve was computed on. Legend
labels carry the NES, the pooled FDR and the testable size. The y range is shared across the
family, so curve heights compare between these nine figures.

The nine sets are the three mouse up arms, the two frozen axes, the curated proteostasis core,
`HALLMARK_HYPOXIA`, the effector-Treg reference, and one translation set for the blood-side
comparison.
*Source* `tables/_overview/runsum_<set>.csv` ·
`02_analysis/scripts/14_unbiased_enrichment_viz.R`.

### `figures/by_contrast/<population>/<COLLECTION>/` — the browse surface

**141 panels over every population × collection cell.**
Layout is `figures/by_contrast/<population>/<COLLECTION>/{dotplot,facet,barplot,running_sum}.png`,
with the rows behind each panel in the mirrored path under `tables/by_contrast/`. Population
directories are Treg, Tcon and CD8, and the contrast inside each is the same donor-paired
synovial-fluid-versus-blood comparison.

Four panel types share their glyphs. **dotplot**: x is GeneRatio, point size is −log10 pooled
adjusted p, fill is NES with orange positive and blue negative, and a black outline marks pooled
FDR < 0.05. **facet** splits the same dotplot into an NES > 0 block and an NES < 0 block.
**barplot** draws NES bars from zero, ordered by NES. **running_sum** stacks the enrichment
curve, the gene-hit ticks and the ranked moderated t, with the score clamped to [−1, 1] so
curves compare between collections.

Start with the three Hallmark dotplots, one per population: fifty named programs on a top-20
axis, with hypoxia and both interferon-response sets among them. The Tcon panel is where the
mouse arm scores highest and where the published STING signature reaches pooled FDR 0.009. The
CD8 panel is where the generic interferon axis carries the single highest signed NES of the
whole sweep.

Two reading rules. The panel types rank by different metrics — adjusted p for dotplot and facet,
|NES| for barplot and running_sum — so read an absence against the rule named in that panel's own
subtitle. Five collections are also too small for the full battery and carry fewer panels. That
omission is a redundancy judgement recorded per collection, and every set of every collection
appears in `gsea_all.csv`.
*Source* `tables/by_contrast/<population>/<COLLECTION>/*.csv` ·
`02_analysis/scripts/14c_gsea_battery_viz.R`.

---

## Tables

### `tables/gsea_all.csv` — the primary table

One row per population and gene set, each identifier appearing exactly once per population, so a
row count is a usable rank denominator. `nes` is positive when the set's genes concentrate on
the synovial-fluid side. `padj` is Benjamini-Hochberg within one collection and `padj_pooled`
across all `n_tests_pooled` tests asked of that ranking. A row significant under `padj` alone
depends on the rest of the sweep going uncounted.

Keep the two per-collection denominators apart: `n_sets_scored_in_db` is what `padj` was
corrected over, and `n_tests_in_db` is what the collection contributes to the pooled family.

### `tables/gsea_pooled_summary_by_db.csv` — the collection-level counts

One row per population and collection. `sig_per_database` counts sets significant under the
within-collection correction and `sig_pooled` under the family-wide one, each against its own
denominator. `sig_lost_to_pooling` goes negative where pooling is looser, which is a real
property of Benjamini-Hochberg. `min_pvalue` and `n_at_min_pvalue` expose a p-value floor with
many ties.

### `tables/gsea_<population>_<database>.csv` — the leading edges

45 files, one per population × collection. Columns match `gsea_all.csv`, with two additions:
the full slash-delimited `core_enrichment`, and the within-collection `padj` corrected over that
file's own rows, which is what a standalone single-collection run reports. `is_pooled_alias` is
TRUE for a set scored here and pooled under another collection. Use these when a leading edge is
needed.

Fewer sets appear than the manifest offers, because a set whose effective size against that
ranking falls outside 5 to 500 is declined. That is a coverage fact about the data.

### The two PROGENy tables and their donor substrate

**`tables/progeny_activity.csv`** — one row per pathway and population, from decoupleR MLM on
the donor-pseudobulk moderated-t statistics. **The `nes` column here holds an MLM activity
statistic**, named for schema compatibility with the fgsea outputs it sits beside. Its
magnitudes run much larger than any normalised enrichment score, so never rank or threshold this
file against a `gsea_*.csv` file. `set_size` and `core_enrichment` describe the model, because
MLM has no leading edge.

**`tables/progeny_donor_activity.csv`** — one row per pseudobulk sample and pathway, 546 rows,
the substrate that lets the direction be tested paired. `activity` is the MLM score against that
population's own centred log-CPM matrix, and `pvalue` is decoupleR's per-sample value.
Comparable within a population, never across.

**`tables/progeny_sf_vs_pb.csv`** — one row per population and pathway. The test is a paired t
on per-donor activities, matched donor for donor, with `n_paired_donors` at six in every
population. `mean_difference` is synovial fluid minus paired blood with a 95% interval, and
`padj` is BH across the fourteen pathways. Read it against `progeny_activity.csv`: agreement
between the contrast-statistic score and the six-donor paired test is what the forest's black
rings encode.

### The five provenance and gate tables

| File | What it holds |
|---|---|
| `tables/geneset_manifest.csv` | One row per collection: sets in the source, sets surviving the nominal 5-to-500 bounds, sets pooled elsewhere as duplicates, and the resolved source with its msigdbr and MSigDB releases. MSigDB 2026.1.Hs serves KEGG as `CP:KEGG_LEGACY`, and the loader records the substitution. |
| `tables/ranked_list_keycheck.csv` | One row per population. `frac_keys_ensembl_like` is the share of ranking keys matching an `ENSG…` pattern, and the script stops above 0.5 with the diagnosis attached. `n_ranked` gives 13,999 / 14,411 / 14,014 symbols. |
| `tables/geneset_overlap.csv` | One row per population and collection: how many genes of the collection's whole union appear in that ranking, as a share of each side. A GO_BP overlap under 2,000 stops the run. |
| `tables/geneset_alias_map.csv` | One row per population and duplicated set. `kept_copy` and `dropped_copy` record which entered the pooled family, `gene_content_identical` is verified set by set, and `abs_nes_difference` is the drift check — a value above 0.01 stops the run. |
| `tables/wt_heat_up_reproduction.csv` | One row per population. `nes_published` against `nes_this_stage` off the same ranking, with equal effective set sizes required for `reproduced` to read TRUE. |

### The running-sum substrates

`tables/runsum_interactive_index.csv` names every emitted substrate table with its collection,
set, effective size, NES and both corrections. Read `always_emitted` as the selection audit:
TRUE means a mouse arm, present unconditionally, and FALSE means a curated set picked on pooled
p-value, split evenly between directions so a synovial-side comparator exists.

`tables/runsum_interactive_<population>_<set>.csv` carries the gene-by-gene walk, one row per
ranked gene: `stat`, `running_es` recomputed off the fitted object's own gene list and exponent,
`hit` and `leading_edge`. These are substrate for an interactive comparison.

### `tables/_overview/<figure stem>.csv`

One file per overview figure, written in the same call that writes the figure, so every panel
ships with the numbers behind it. Each is deliberately narrowed to what is drawn — read the
stage tables above for anything wider.
