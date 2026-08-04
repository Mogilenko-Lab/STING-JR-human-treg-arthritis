# 14_unbiased_enrichment: artifact captions

_**Abbreviations:** SF = synovial fluid, PB = peripheral blood, NES = normalized enrichment score, FDR = Benjamini-Hochberg adjusted p-value, MLM = multivariate linear model._

Everything this compartment has published about the inflamed juvenile-arthritis joint so far asks a **targeted** question of the same contrast: I hand the sorted JIA synovial-fluid-versus-paired-blood ranking one named, mouse-derived signature and ask whether it enriches. It does. But a targeted test can only answer about the set it was handed, and it cannot tell a reader whether that answer is remarkable. So here I ask the unbiased counterpart of the same three frozen ranked lists: **what does this niche contrast actually contain, across curated databases, with no set privileged?**

That question is the honest complement to the targeted one, and the two are only informative together. If a hundred inflammatory programs move as strongly as the mouse-derived arm, the arm is one of many co-enriching programs and its enrichment says something about the niche rather than something specific about the arm. If almost nothing else moves that strongly, the arm is distinctive. This stage exists to say which of those is the case, and nothing more.

Every number here is a ranked-list enrichment statistic or a signalling-footprint activity score on the very same donor-level pseudobulk contrast the confirmatory answer already rests on. No row is written to `03_results/master/` and no row becomes an effect size. Language stays correlative throughout: a set enriching means its gene content moves with the synovial-fluid side of this ranking, not that the program the set is named for is present, and certainly not that it is driving anything.

### What was scored, and with what parameters

Two methods, chosen because they fail differently. Pre-ranked fgsea needs a gene-set list and therefore inherits every curation and size decision baked into that list; decoupleR MLM on PROGENy needs no list at all and so cannot inherit them.

The fgsea run parameters are not new: they are the ones the targeted signature test in `05_scoring` already used, read from `02_analysis/config/analysis_config.yaml` — `clusterProfiler::GSEA(by = "fgsea")` with `exponent = 1`, `eps = 0` (exact multilevel p-values), `nPermSimple = 100000`, `pvalueCutoff = 1` so every set is reported rather than pre-filtered, `pAdjustMethod = "BH"`, `seed = TRUE` at `gsea_seed = 123`, and effective set-size bounds `gsea_min_size = 5` to `gsea_max_size = 500`. Using this compartment's own frozen `gsea_min_size = 5` rather than the mouse anchor's 15 is deliberate: 5 is the value the published contrast was scored at, and it keeps the small frozen axes (`sting_specific_up`, 12 to 15 genes in these rankings) testable instead of silently dropping them.

| Collection | Source | Sets tested |
|---|---|---|
| `GO_BP` | MSigDB C5 GO:BP, msigdbr 26.1.0 / MSigDB 2026.1.Hs | 7,276 |
| `GO_MF` | MSigDB C5 GO:MF | 1,825 |
| `Reactome` | MSigDB C2 CP:REACTOME | 1,818 |
| `GO_CC` | MSigDB C5 GO:CC | 1,026 |
| `WikiPathways` | MSigDB C2 CP:WIKIPATHWAYS | 924 |
| `TF_Targets` | CollecTRI human regulons, unsigned, one set per transcription factor | 763 |
| `KEGG` | MSigDB C2 CP:KEGG_LEGACY | 186 |
| `Hallmark` | MSigDB H | 50 |
| `project_frozen` | this compartment's six frozen MSigDB Hallmark lists plus the frozen curated `HSR_core` | 7 |
| `mouse_projection` | the mouse-derived human-projected UP arms: `WT_heat_up`, `KO_heat_up`, `Interaction_up` | 3 |
| `sting_axes` | the frozen SAVI reference axes `sting_specific_up` and `ifn_only_up` | 2 |

Three things in that table a reviewer should check rather than take.

**KEGG is the legacy collection.** MSigDB 2026.1.Hs renamed `CP:KEGG` to `CP:KEGG_LEGACY` and added a much larger `CP:KEGG_MEDICUS`. The loader retries with the legacy name so the collection matches the mouse anchor's, and `geneset_manifest.csv` records the substitution rather than hiding it.

**The TF collection is not the mouse anchor's.** The mouse sweep uses MSigDB C3 TFT:GTRD. On the human side this repository already carries a curated CollecTRI regulon table, so that is used instead, pooling each factor's activating and repressing targets into one unsigned set. The signed use of the same network is decoupleR ULM — a different method with a different readout, and mixing the two would be quiet.

**Down arms are out of scope by design**, so `mouse_projection` carries three sets and not six.

### The reproduction gate, and the silent failure it guards

The count matrices in this compartment are keyed by Ensembl id while every reference gene set matches on HGNC symbol. A ranked list that leaked Ensembl ids would intersect all eleven collections at approximately zero, and fgsea reports that as empty or NA rows rather than as an error — it looks exactly like a biological null. Two checks stand in front of that. The ranked-list keys are tested before anything runs (`ranked_list_keycheck.csv`: 0% of keys are Ensembl-like in all three populations, 13,999 / 14,411 / 14,014 symbols), and the per-collection gene overlap is published as a first-class table rather than assumed (`geneset_overlap.csv`: the GO_BP union meets 9,352 of the Treg ranking's symbols, 9,648 of Tcon's, 9,411 of CD8's).

The harder check is whether this sweep is on the same footing as the targeted result it exists to calibrate. `WT_heat_up` is run as an ordinary member of the sweep, and its NES must land on the published value off the same ranked list. It does, to within fgsea's own stochastic normalisation, with identical effective set sizes:

| Population | Published NES | This sweep | Difference | Effective set size |
|---|---|---|---|---|
| Treg | +2.591524 | +2.591633 | 1.1 × 10⁻⁴ | 119 = 119 |
| Tcon | +2.680878 | +2.680206 | 6.7 × 10⁻⁴ | 130 = 130 |
| CD8 | +2.070985 | +2.071197 | 2.1 × 10⁻⁴ | 113 = 113 |

The script stops before touching the large collections if this fails, and it was never tuned to agree.

A second, unplanned consistency check falls out of the design: six sets are scored twice in every population, once from `msigdbr` and once from this compartment's independently frozen copy. In Treg `HALLMARK_IL2_STAT5_SIGNALING` gives +2.5891 and +2.5874 on the same 167 effective genes, and across all eighteen duplicate pairs the largest disagreement is 5.0 × 10⁻³ in NES — inside fgsea's stochastic normalisation. The frozen lists have not drifted. That check is worth keeping and is recorded in `geneset_alias_map.csv`; what it must not do is let one hypothesis into the multiplicity correction twice, which is the subject of *One set, one hypothesis* below.

### What the sweep returns

**The niche contrast moves an enormous amount at once.** Of 11,236 tests asked of the Treg ranking, 1,490 reach FDR < 0.05 after pooling — 1,423 toward synovial fluid and 67 toward blood. Tcon is larger still (2,165 of 11,459) and CD8 smaller (1,027 of 11,242). Fifty-six per cent of MSigDB Hallmark is significant in Treg. This is what an inflamed tissue niche against paired circulating cells looks like, and it is the context every targeted result in this compartment should be read against.

**The mouse-derived up arm is strong, and how ordinary it is depends on the compartment.** `WT_heat_up` reaches NES +2.59 at pooled FDR 3.7 × 10⁻¹² in Treg, and 21 other pooled-significant sets carry a larger absolute NES. Ranking on absolute NES is the fair comparison — a signed ranking would flatter the arm by ignoring everything moving the other way — but the raw count of 21 needs one qualification before it can be read, because **19 of those 21 are on the blood side**. Conditioned on direction, the picture separates by compartment:

| Population | Larger \|NES\| among significant sets | of which blood-side | Rank among synovial-fluid-side significant sets |
|---|---|---|---|
| Treg | 21 | 19 | **3rd of 1,423** |
| Tcon | 30 | 29 | **2nd of 2,029** |
| CD8 | 141 | 68 | 74th of 897 |

So in the two CD4 compartments the arm is close to the strongest thing moving toward synovial fluid — only `REACTOME_INTERLEUKIN_10_SIGNALING` (+2.64) and `GOCC_COPII_COATED_ER_TO_GOLGI_TRANSPORT_VESICLE` beat it in Treg, and only `HALLMARK_INTERFERON_GAMMA_RESPONSE` in Tcon. In CD8 it is genuinely mid-pack, with 73 synovial-fluid-side sets above it that are dominated by MHC class II and antigen presentation (`REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION`, the `RFX5` / `RFXAP` / `RFXANK` regulons) and by interferon.

Both readings therefore hold, and neither on its own is honest. The arm is not a *uniquely* distinctive feature — a great deal co-enriches with it, and it does not outrank the largest shifts in the contrast, which are downward. But among programs moving toward the inflamed joint in sorted CD4 cells it is at the top of the distribution rather than inside it. That is compatible with what this compartment established by a different route — its synovial enrichment is not cell-subset-selective, and by curated composition its up arm is largely inflammatory — while adding the part a targeted test could not supply: what it is at the top *of*.

**`KO_heat_up` tracks it almost exactly** (+2.56 Treg, +2.65 Tcon, +2.08 CD8, ranks 5 / 3 / 63 among synovial-fluid-side significant sets), which is the arithmetic one should expect from arms that are linearly dependent by construction, and is a further reason not to read the WT arm as carrying anything cGAS-specific. **`Interaction_up` reaches nothing anywhere** (+1.47 / +1.40 / +1.54, pooled FDR 0.24 / 0.26 / 0.17): at six genes in the ranked list it is under any useful power floor, and its non-significance is a statement about its size, not about its biology.

**The largest single shifts in this contrast are downward, and no mouse arm speaks to them.** Ranked by pooled FDR, the top ten Treg results are all translation and ribosome sets moving toward blood — `WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS` at NES −3.46 (4.6 × 10⁻³²), `REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION` at −3.39, `KEGG_RIBOSOME` at −3.35, `GOBP_CYTOPLASMIC_TRANSLATION` at −3.10. These magnitudes exceed anything on the synovial-fluid side. Only 67 of the 1,490 pooled-significant Treg sets are down-going, so this is a small number of very large effects rather than a broad loss. It is also entirely invisible to a targeted test of up arms, which is the clearest argument for running this stage at all.

**Two frozen axes behave as their own compartment predicts.** The generic type-I interferon axis `ifn_only_up` enriches robustly (+2.32 Treg at 1.7 × 10⁻⁶, +2.34 Tcon, +2.61 CD8) while the 21-gene STING-specific axis does not clear the pooled threshold in Treg (+1.52 at 0.195, on 12 ranked genes) and reaches it only in Tcon (+1.87 at 0.018, on 15). A reading in which this niche contrast is interferon-like but not STING-specific survives; a STING-specific reading does not. And the curated `HSR_core` proteostasis lens sits at +1.49 in Treg at pooled FDR 0.146, negative in Tcon (−1.34) and CD8 (−1.15) — a Treg-versus-others sign difference at trend level, which is where it stood before and where it stays.

### What pooling across the family does, and why both corrections are published

`padj` is Benjamini-Hochberg within one collection, which is what a single-collection run reports and what makes a row comparable to a published per-collection number. `padj_pooled` corrects across every test asked of that population's ranked list, which is the honest correction for a sweep that interrogates one ranking eleven times.

Pooling is **not** a uniform tightening, and the tables show that plainly rather than implying otherwise. In Treg the total barely moves (1,492 significant per-database against 1,490 pooled) but it redistributes: the small hand-picked collections lose (Hallmark 31 → 28, KEGG 41 → 37, Reactome 264 → 249, `TF_Targets` 187 → 152) while GO_BP gains (663 → 694). Benjamini-Hochberg divides by rank as well as multiplying by family size, so a set whose rank rises faster than the family grows comes out marginally tighter — 52% of rows do.

The effect that matters is the first one: a small, deliberately chosen collection no longer benefits from having been tested on its own. That is exactly the correction a targeted test cannot apply to itself.

### One set, one hypothesis — the pooled family is deduplicated

Two collections here legitimately hold the same gene set. `project_frozen` re-pins six MSigDB Hallmark sets to files so the decomposition and purge stages have a size-validated asset that cannot move under an `msigdbr` upgrade, and the `Hallmark` collection then fetches those same six live. The gene content is identical — verified set by set, not assumed — so letting both copies into one pooled Benjamini-Hochberg family put six exact duplicate hypotheses in it and made every population's rank denominator six too high. The correction is now applied to a family that contains each hypothesis once: **11,236 Treg, 11,459 Tcon, 11,242 CD8**, each six lower than the figure this stage previously reported.

Both copies are still scored, because their agreement is a real check on whether the frozen files have drifted from the package they were frozen from, and both still appear in their own per-database table so each of those reads as a standalone single-collection run. Only the canonical copy enters the pooled family, and `geneset_alias_map.csv` records which copy was kept, which was dropped, and both copies' statistics side by side. Canonicality goes to MSigDB, since it is the home of a `HALLMARK_` identifier and a frozen file here is the re-pin of it.

The resolution is structural rather than a list of set names: the loader detects any identifier present in more than one collection, refuses to proceed if the two copies' genes differ, and refuses to demote a collection that is not one of the three file-backed ones. Adding a seventh frozen Hallmark set to the configuration is therefore deduplicated without editing code.

**What this did and did not move.** Nothing about the underlying tests changed: `nes`, `pvalue` and per-collection `padj` are bit-for-bit identical to the previous run for every surviving row, and the Spearman rank correlation of `padj_pooled` old against new is 1.000 in all three populations. Pooled q-values shift in the third or fourth significant figure, at most by 6.8%.

No set became significant. Twelve borderline rows between q = 0.0498 and q = 0.0499 crossed to just above 0.05, none of them a set this compartment's reading names. One headline number moved: `WT_heat_up` in CD8 goes from 75th to 74th by signed NES, because a removed duplicate (`HALLMARK_INTERFERON_ALPHA_RESPONSE`, +2.25) had been counted above it. This is a multiplicity-family correction, not a change of result.

One reporting consequence is visible in the overview figure: the `project_frozen` row now reads **0 of 1**, because six of that collection's seven sets are reported under `Hallmark` and only `HSR_core` remains in the pooled family — where it does not reach significance. The row is thin but the numerator and denominator are on the same basis, which they were not before.

### The gene-set-free read

decoupleR MLM on the human PROGENy model (14 pathways, continuous weights, top 500 targets each) is run two ways on the same contrast: on the donor-pseudobulk moderated-`t` statistics, mirroring how the mouse anchor runs it, and on the per-donor pseudobulk activities so the direction can be tested paired across the six donors who span both arms.

In Treg the largest footprint is JAK-STAT (+9.40, FDR 8.6 × 10⁻²⁰), then EGFR (+8.25), TGFβ (+5.60), Hypoxia (+5.40, 2.4 × 10⁻⁷), PI3K (+5.22), VEGF (+3.71) and NF-κB (+3.49), with WNT the only negative (−2.88). Eight of the fourteen are significant on both tests. The pattern is the same in all three sorted populations, which again says this is a property of the niche rather than of a cell subset.

The Hypoxia footprint rising alongside the inflammatory ones is worth stating carefully, because it is easy to over-read. Both readouts move together in the same cross-sectional contrast, and that is all these data show. Whether the low-oxygen and inflammatory readouts stand in a causal relation, in either direction, is not something a synovial-fluid-versus-blood comparison can decide, and the joint imposition of temperature and low oxygen by an inflamed joint is not separable here either.

### Reading order

`gsea_all.csv` is the one table to start from — every test, both corrections, the family size, one row per set. `gsea_pooled_summary_by_db.csv` gives the collection-level counts the overview figure draws. The per-collection `gsea_<population>_<database>.csv` files carry the leading-edge gene lists that `gsea_all.csv` reduces to a count, and are the only place the duplicate copies of the six re-pinned Hallmark sets appear; `geneset_alias_map.csv` says which copy the pooled family kept. The `runsum_interactive_*` tables are substrate for an interactive running-sum comparison and are not meant to be read as tables.

## figures/_overview/pooled_overview_by_population.png

The synovial-fluid-versus-paired-blood contrast moves a great many
curated programs at once — 1,490 of 11,236 tests reach FDR < 0.05 in
Treg after pooling — and while 21 pooled-significant sets carry a
larger absolute NES than the mouse-derived WT_heat_up arm (NES +2.59,
pooled FDR 4e-12), 19 of those 21 are on the blood side, so among the
1,423 sets moving toward synovial fluid the arm ranks 3.

**How to read:** Columns are the three sorted populations side by side. Below the
dashed line each row is a collection, ordered by how many of its sets
reach significance; above it each row is one mouse-derived up arm. A
small point is a set at FDR < 0.05 after pooling across every test
asked of that ranking. Brown concentrates on the synovial-fluid side,
blue on paired blood. Horizontal position is the exact NES, clamped
to plus or minus 3.5. Yellow diamonds are the arms, filled when
significant. Grey text is that population's own count, or an arm's
NES and FDR. Read for calibration: far right in a dense row is strong
but ordinary. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; gsea_fdr=0.05; nes_cap=3.5; padj_pooled_method=BH` | `03_results/14_unbiased_enrichment/tables/{gsea_all,gsea_pooled_summary_by_db}.csv` |

## figures/_overview/treg_top_sets.png

In the JIA Treg contrast the largest shifts are downward rather than
upward: WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS reaches NES -3.46 toward
paired blood, against REACTOME_INTERLEUKIN_10_SIGNALING at +2.64
toward synovial fluid, so the niche difference this compartment reads
as an inflammatory gain is accompanied by an at least equally large
loss of translation and ribosomal programs.

**How to read:** One row per gene set, capped at the top ten in each direction by
absolute NES among sets at FDR < 0.05 after pooling; the subtitle
states how many are not shown, so the cap is not completeness.
Identifiers are shown with underscores as spaces and wrapped, never
truncated, each with its collection in brackets. Right of zero the
set's genes concentrate on the synovial-fluid side of this ranking,
left on paired blood. Point size is how many of the set's genes reach
the ranked list, so a large NES on a small point rests on few genes.
The grey number is the pooled FDR; a black ring marks a mouse-derived
arm. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment_viz.R` | `main` | `figures.top_n=20 (split 10 per direction); gsea_fdr=0.05; gsea_min_size=5; gsea_max_size=500` | `03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/_overview/progeny_activity_panel.png

Read without any gene-set list, the JIA Treg
synovial-fluid-versus-blood contrast carries its largest PROGENy
footprint in JAK-STAT (score +9.40, FDR 9e-20) while the Hypoxia
footprint scores +5.40 at FDR 2e-07, so the inflammatory and
low-oxygen readouts rise together in the same niche contrast and
neither can be read as the other's cause; 8 of the fourteen pathways
are significant on both tests in Treg.

**How to read:** One row per PROGENy pathway, ordered by its Treg score, with a grey
line joining the three populations so the spread within a row is the
between-population difference. Horizontal position is the model
activity score on the donor-pseudobulk moderated-t contrast
statistics: right of zero the pathway's footprint genes move with
synovial fluid, left with paired blood. Colour is the population. A
solid point reaches FDR < 0.05 on that test; a black ring marks one
that also reaches it in the independent donor-paired test, so ringed
and solid is corroborated twice and faded and unringed by neither. A
footprint is inferred from target-gene expression, not a measurement
of pathway activity. Correlative; no causal reading.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment_viz.R` | `main` | `progeny.organism=Human; progeny.top=500; progeny.minsize=5; gsea_fdr=0.05` | `03_results/14_unbiased_enrichment/tables/{progeny_activity,progeny_sf_vs_pb}.csv` |

## tables/geneset_manifest.csv

Eleven collections holding 13,880 sets were scored, and the KEGG collection that MSigDB 2026.1.Hs actually served is `CP:KEGG_LEGACY` (186 sets) rather than the `CP:KEGG` the configuration names.

**How to read:** One row per collection. `n_sets_in_source` is what the source returned and `n_sets_after_nominal_size_filter` what survived the 5-to-500 nominal bounds; that difference is a curation fact about the collection, and a set can still be declined later for its EFFECTIVE size after meeting a ranked list. `nominal_size_filter_applied` is FALSE for the three file-backed collections, so a small frozen set is reported as untestable with its size rather than removed here. `n_sets_aliased_out_of_pooling` counts sets another collection also holds and therefore pools instead (6 for `project_frozen`, 0 elsewhere) and `n_sets_offered_for_pooling` is the remainder; both are nominal, so the effective count is in `gsea_pooled_summary_by_db.csv`. `source` records the resolved collection, any substitution the loader made, and the msigdbr and MSigDB releases. Provenance table; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `load_msigdb, load_collectri, read_set_file, filter_by_size` | `gsea_min_size=5; gsea_max_size=500; project.species=Homo sapiens` | `msigdbr 26.1.0`, `../mouse_anchor/00_data/references/networks/CollecTRI_regulons_human.csv`, `00_data/references/{msigdb_hallmark,temp_hsr_lens}/*.txt`, `../mouse_anchor/03_results/human_projection/signatures/*/*_up.txt`, `../sting_positive_control/03_results/06_reference_axis/signatures/*_up.txt` |

## tables/ranked_list_keycheck.csv

All three ranked lists are keyed by HGNC symbol and none by Ensembl id (0% Ensembl-like keys), so the enrichment results below are not the silent near-zero-overlap null that an Ensembl-keyed list produces.

**How to read:** One row per sorted population. `frac_keys_ensembl_like` is the share of ranking keys matching an `ENSG…` pattern; the script stops hard above 0.5 with the diagnosis attached, because fgsea returns empty or NA rows for a failed symbol join rather than an error and the result reads as a biological null. `n_ranked` is the number of unique symbols, `first_key` and `last_key` the extremes of the descending signed moderated-`t` ranking, and `stat_min` / `stat_max` its range. Diagnostic; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `read_ranked, keycheck_row` | `paths.results=03_results/` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` |

## tables/geneset_overlap.csv

Every collection meets its ranked list at thousands of genes — the GO_BP union at 9,352 of the Treg ranking's 13,999 symbols — which is the positive evidence that the symbol join joined.

**How to read:** One row per population and collection. `n_overlap` is how many genes of the collection's whole union appear in that ranking; `frac_of_set_genes` and `frac_of_ranked` express it as a share of each side. Read this table alongside `ranked_list_keycheck.csv`: the key check says the ranking looks right and this table says the reference sets actually reach it. A GO_BP overlap under 2,000 stops the run. The mouse-derived and SAVI rows are small because those collections are small (145 and 72 genes in Treg), which is expected and not a failure. Diagnostic; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `dplyr::bind_rows over intersect (overlap section)` | `gsea_min_size=5; gsea_max_size=500` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/objects/14_genesets.rds` |

## tables/geneset_alias_map.csv

Six gene sets are held by two collections at once — the MSigDB Hallmark originals and this compartment's frozen re-pins of them — and this table records that each is scored twice, which copy the pooled multiplicity correction kept, and that the two copies agree to at most 5.0 × 10⁻³ in NES.

**How to read:** One row per population and duplicated set, eighteen rows. `kept_copy` is the collection whose copy entered the pooled family and `dropped_copy` the one excluded; the dropped copy is still scored and still in its own `gsea_<population>_<database>.csv`, so nothing is deleted, only counted once. `gene_content_identical` is verified set by set, not assumed — two collections holding one identifier over different genes stops the run, since the identifier would no longer name one hypothesis. `abs_nes_difference` is what to read this table for: the disagreement between copies is fgsea's stochastic normalisation and nothing more, and a value above the 0.01 tolerance also stops the run, because then dropping a copy would discard evidence rather than a duplicate. Bookkeeping plus a drift check on the frozen lists; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `is_pooling_alias (alias-resolution section) + the alias record in the pooling section` | `reproduction_check.nes_tolerance=0.01; padj_pooled_method=BH` | `03_results/objects/14_genesets.rds`, `03_results/objects/14_gsea/*.rds` |

## tables/wt_heat_up_reproduction.csv

Scored as an ordinary member of this sweep, `WT_heat_up` lands on the published NES in every population to within 1.1 × 10⁻⁴, 6.7 × 10⁻⁴ and 2.1 × 10⁻⁴ on identical effective set sizes, so the whole sweep is on the same footing as the targeted result it exists to calibrate.

**How to read:** One row per population. `nes_published` is read from the already-published table for that population and `nes_this_stage` from this sweep, off the same ranked list; `abs_nes_difference` must fall under `nes_tolerance` and the two effective set sizes must be equal for `reproduced` to be TRUE. A residual difference of order 10⁻⁴ is fgsea's own stochastic normalisation, not a method difference. The p-values differ slightly for the same reason and are shown so that agreement is not claimed on a statistic that was not compared. The run stops before touching the large collections if any row is FALSE, and was not tuned to agree. Verification table; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `run_cell (reproduction-gate section)` | `reproduction_check.signature=WT_heat_up; reproduction_check.nes_tolerance=0.01; gsea_seed=123; gsea_nperm=100000` | `03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv`, `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv` |

## tables/gsea_all.csv

Of 11,236 tests asked of the Treg ranking 1,490 reach FDR < 0.05 after pooling (1,423 toward synovial fluid, 67 toward blood), and the mouse-derived `WT_heat_up` arm ranks 31st among them by pooled q with 21 pooled-significant sets carrying a larger absolute NES.

**How to read:** The primary table — one row per population and gene set, 33,937 rows. **Each identifier appears exactly once per population**, so a row count is a usable rank denominator; the six Hallmark sets `project_frozen` re-pins are scored in both collections but pooled only under `Hallmark`, per `geneset_alias_map.csv`. `nes` is positive when the set's genes concentrate on the synovial-fluid side. Two corrections answer different questions: `padj` is Benjamini-Hochberg WITHIN one collection, what a single-collection run reports; `padj_pooled` corrects across all `n_tests_pooled` tests asked of that ranking. A row significant only under `padj` depends on not counting the rest of the sweep. Pooling is not uniformly stricter, so compare the two per collection rather than assuming a direction. The two per-collection denominators must not be swapped: `n_sets_scored_in_db` is what `padj` was corrected over, `n_tests_in_db` what the collection contributes to the pooled family. Correlative; calibration is the intended use.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `tidy_cell + stats::p.adjust (pooling section)` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000; gsea_fdr=0.05; padj_pooled_method=BH` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/objects/14_genesets.rds` |

## tables/gsea_pooled_summary_by_db.csv

Pooling across the whole family barely changes the Treg total (1,492 significant per-collection to 1,490 pooled) but redistributes it, costing the small deliberately chosen collections (Hallmark 31 to 28, KEGG 41 to 37, Reactome 264 to 249, the CollecTRI regulons 187 to 152) and gaining in GO_BP (663 to 694).

**How to read:** One row per population and collection. `sig_per_database` counts sets significant under the within-collection correction, `sig_pooled` under the family-wide one, and each is divided by its own denominator: `n_sets_scored_in_db` is what `padj` was corrected over, `n_tests_in_db` what the collection contributes to the pooled family after duplicates are resolved. For `project_frozen` those are 7 and 1, six of its sets being pooled under `Hallmark` (`n_sets_aliased_out_of_pooling`). `sig_lost_to_pooling` counts only sets lost by facing the wider family, and goes NEGATIVE where pooling is looser — a real property of Benjamini-Hochberg, since a set whose rank rises faster than the family grows comes out tighter. The row that matters for calibration is the small-collection one: a collection of fifty no longer benefits from having been tested alone. `min_pvalue` and `n_at_min_pvalue` expose a p-value floor with many ties. Enrichment statistics; correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `dplyr::summarise (pooling section)` | `gsea_fdr=0.05; padj_pooled_method=BH` | `03_results/14_unbiased_enrichment/tables/gsea_all.csv` (same in-memory table) |

## tables/gsea_&lt;population&gt;_&lt;database&gt;.csv

The per-collection results carrying the leading-edge gene lists: 33 files, one for each of the three sorted populations crossed with the eleven collections, from `gsea_treg_Hallmark.csv` (50 sets) to `gsea_treg_GO_BP.csv` (5,766 sets scored of 7,276 offered).

**How to read:** A family of files; `<population>` ranges over `treg`, `tcon`, `cd8` and `<database>` over the eleven collections in `geneset_manifest.csv`. Columns match `gsea_all.csv` except that these carry the full `core_enrichment` — slash-delimited leading-edge gene names — and the within-collection `padj` only, corrected over that file's own rows, which is what a standalone single-collection run reports. `is_pooled_alias` is TRUE for a set scored here but pooled under another collection: the six Hallmark sets in `gsea_<population>_project_frozen.csv`. Those are the only rows in this stage where one set's statistics appear twice, and their agreement is the drift check described above. Rows sort by p-value then descending absolute NES. Fewer sets appear than the manifest offers, because `clusterProfiler` declines a set whose EFFECTIVE size against that ranked list falls outside 5 to 500 — a coverage fact about the data. Use these when a leading edge is needed. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `run_cell, tidy_cell` | `gsea_min_size=5; gsea_max_size=500; gsea_seed=123; gsea_nperm=100000` | `03_results/03_pseudobulk/tables/ranked_{treg,tcon,cd8}.tsv`, `03_results/objects/14_genesets.rds` |

## tables/runsum_interactive_index.csv

Twenty-seven running-sum substrate tables were emitted — the three mouse-derived up arms in every population plus six curated comparators per population, split evenly between directions so an up-going curated set exists to compare against the up-going mouse arms.

**How to read:** One row per emitted substrate table, naming its file, collection, gene set, effective set size, NES and both corrections. Read `always_emitted` as the selection audit: TRUE means a mouse arm, present unconditionally; FALSE means a curated set that was PICKED on pooled p-value, so its presence carries a selection step and its NES is not the best in its collection. The even split between directions is deliberate — ranked on p-value alone every curated slot went to a paired-blood-side translation or ribosome set, leaving no synovial-fluid-side comparator. Index; no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `runsum_table (substrate section)` | `runsum_top_curated=6 (split 3 per direction); gsea_fdr=0.05` | `03_results/14_unbiased_enrichment/tables/gsea_all.csv` (same in-memory table), `03_results/objects/14_gsea/*.rds` |

## tables/runsum_interactive_&lt;population&gt;_&lt;set&gt;.csv

The gene-by-gene substrate behind each running-sum curve — for Treg, the three mouse-derived up arms at NES +2.59, +2.56 and +1.47 beside curated comparators spanning `REACTOME_NEUTROPHIL_DEGRANULATION` at +2.40 and `WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS` at −3.46.

**How to read:** A family of 27 files, one per row of `runsum_interactive_index.csv`; substrate for an interactive comparison, not a table to read. One row per gene of that population's ranking, in rank order: `stat` is the signed moderated `t`, `running_es` the weighted running enrichment score recomputed with the DOSE formula off the fitted object's own gene list and exponent, so the emitted curve IS the plotted one; `hit` marks a set member and `leading_edge` a member of the object's core enrichment. The nine columns and their order match the published substrate tables exactly, and the script asserts it. Filename-unsafe characters are substituted in the FILE name only. The curve shows WHERE a set concentrates, not how much. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `runsum_table` | `runsum_top_curated=6 (split 3 per direction)` | `03_results/objects/14_gsea/*.rds`, `03_results/objects/14_genesets.rds` |

## tables/progeny_activity.csv

Read with no gene-set list at all, every sorted population shows the same footprint pattern in the synovial-fluid-versus-blood contrast — JAK-STAT largest (+9.40 Treg, +8.89 Tcon, +10.85 CD8), then EGFR, TGFβ and Hypoxia, with WNT the only negative.

**⚠ The `nes` column in this file does NOT hold a normalized enrichment score.** It holds a `decoupleR::run_mlm` activity statistic and is **not comparable to the `nes` values in the sibling `gsea_*.csv` files**. Its magnitudes run much larger: the fgsea sweep spans −3.76 to +2.75 over 33,955 rows, these 42 rows span −3.87 to +10.85, and 15 of them exceed the largest absolute NES anywhere in the sweep. JAK-STAT at +9.4 in Treg is an ordinary MLM score and an impossible NES.

The column is named `nes` **deliberately**, for schema compatibility with `master_gsea_table` and the fgsea outputs it sits beside — not a bug. Never rank or threshold this file against a `gsea_*.csv` file.

**How to read:** One row per PROGENy pathway and population, from decoupleR MLM on the donor-pseudobulk moderated-`t` contrast statistics. `nes` (see the warning above) is positive when the pathway's weighted target genes move with synovial fluid. `padj` is Benjamini-Hochberg across the fourteen pathways within a population, so the correction is mild by construction. `set_size` and `core_enrichment` are contrast-invariant footprint membership — MLM has no leading edge, so they describe the model, not the result. The three populations do not share a gene universe, so a gene absent from one enters as `t = 0`; padding counts are logged. A footprint is inferred from target-gene expression, not measured. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `decoupleR::run_mlm (PROGENy contrast section)` | `progeny.organism=Human; progeny.top=500; progeny.minsize=5` | `03_results/03_pseudobulk/tables/de_SFvsPB_{treg,tcon,cd8}.csv`, `progeny::getModel` |

## tables/progeny_donor_activity.csv

Per-donor footprint activities for all 39 donor-by-tissue-by-population pseudobulk samples, which is what lets the direction in `progeny_activity.csv` be tested paired rather than read off a single contrast statistic.

**⚠ `activity` is a `decoupleR::run_mlm` statistic, not an enrichment score.** It carries the same units as the misleadingly named `nes` column of `progeny_activity.csv` and shares no scale with any `gsea_*.csv` column.

**How to read:** One row per pseudobulk sample and pathway, 546 rows. `activity` is the MLM score for that sample against the population's own centred log-CPM matrix; `pvalue` is decoupleR's per-sample value, NOT the SF-versus-PB test, which lives in `progeny_sf_vs_pb.csv`. Three choices behind the numbers: counts are joined to symbols through the committed Ensembl-to-HGNC map with duplicated symbols SUMMED, since these are counts from distinct ids of one gene; the retained genes are that population's own differential-expression universe, so both PROGENy arms see the same genes; and log-CPM is row-centred WITHIN a population. Activity is comparable within a population, never across. Correlative; annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `decoupleR::run_mlm (PROGENy donor-level section)` | `progeny.organism=Human; progeny.top=500; progeny.minsize=5; design.tissue_levels` | `03_results/03_pseudobulk/tables/{pseudobulk_counts,pseudobulk_coldata,gene_symbols}.csv`, `progeny::getModel` |

## tables/progeny_sf_vs_pb.csv

Tested paired across the six JIA donors who contribute both arms, ten of fourteen footprints separate synovial fluid from blood in Treg — EGFR most strongly (mean difference +10.9, FDR 2.7 × 10⁻⁴), with Hypoxia (+7.4, 1.2 × 10⁻³) and JAK-STAT (+9.0, 0.028) among them.

**⚠ Nothing in this file is an enrichment score.** `mean_difference` is a difference of `decoupleR::run_mlm` activities, sharing no scale with the `nes` column of any `gsea_*.csv` file — nor with the identically named column in `progeny_activity.csv`, which is also an MLM statistic.

**How to read:** One row per population and pathway. The test is a PAIRED `t` test on per-donor activities, matched donor for donor, because the design is paired; `n_paired_donors` is six in every population here. `mean_difference` is synovial fluid minus paired blood in that population's centred activity units with a 95% interval, and `padj` is Benjamini-Hochberg across the fourteen pathways. Read this AGAINST `progeny_activity.csv`, not instead of it: agreement between the contrast-statistic score and the six-donor paired test is what the figure's black rings encode, and a pathway significant on only one is weaker than either suggests. At six pairs a null is not evidence of absence. Correlative; annotation tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment.R` | `stats::t.test (PROGENy donor-level section)` | `progeny.minsize=5; design.tissue_levels.synovial_fluid=synovial_fluid; design.tissue_levels.peripheral_blood=peripheral_blood` | `03_results/14_unbiased_enrichment/tables/progeny_donor_activity.csv` (same in-memory table) |

## tables/_overview/&lt;figure stem&gt;.csv

The three source tables sitting beside the three figures, each holding exactly the rows its panel draws: `pooled_overview_by_population.csv`, `treg_top_sets.csv` and `progeny_activity_panel.csv`.

**How to read:** A family of three files, one per figure, written in the same call that writes the figure so a panel cannot exist without the numbers behind it. Each is a projection of a table above — the overview table is every pooled-significant row plus all mouse-arm rows with their collection-level counts joined on; the top-sets table is the twenty plotted Treg rows; the PROGENy table joins the contrast-statistic and donor-paired results per pathway with the two significance flags the glyphs encode. Read one when checking a specific mark on its figure; read the stage tables above for anything wider, since these are deliberately narrowed to what is drawn. Correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment_viz.R` | `save_overview` | `figures.top_n=20; gsea_fdr=0.05; nes_cap=3.5` | `03_results/14_unbiased_enrichment/tables/{gsea_all,gsea_pooled_summary_by_db,progeny_activity,progeny_sf_vs_pb}.csv` |

## figures/_overview/arm_nes_by_cell_state.png

The mouse 39 °C-derived up arm rises on the synovial-fluid side of the
paired contrast in all three sorted cell states, at NES 2.5916 in Treg
(119 of 199 arm genes reaching the ranked list), 2.6802 in Tcon (130)
and 2.0712 in CD8 (113), every one below pooled FDR 1e-4. The Treg
score sits between the Tcon and CD8 scores, so the separation reads as
pan-T. KO_heat_up tracks it row for row (NES 2.5645 in Treg on 132 of
218 genes). The 7-gene interaction arm reaches NES 1.4030 to 1.5425 on
6 testable genes and clears pooled FDR 0.05 in none of the three cell
states, so at that size it carries no direction here.

**How to read:** One dot per mouse-derived up arm and cell state, at the confirmatory
tier: donor-level pseudobulk within frozen sort labels, limma-voom
moderated t, then pre-ranked fgsea. Rows are the three arms; inside a
row the three cell states are offset vertically and coloured, each
with its own annotation line. The x position is the normalised
enrichment score for synovial fluid over paired blood. A filled dot
clears the config FDR threshold of 0.05 and an open dot sits above it.
The annotation column gives how many of the arm's genes reached that
population's ranked list, against how many the frozen arm holds, then
the adjusted p. Read that count with the score: an arm of 6 testable
genes and one of 130 are not equally resolved. The adjusted p is
Benjamini-Hochberg pooled across every set that population's sweep
tested; the same-stem source table also carries the per-collection
value over the three arms alone, and the two agree except for
Interaction_up in CD8, 0.035 per-collection against 0.172 pooled. A
score of this kind has no interval, so none is drawn. An arm rising
here means its gene content moves with the synovial-fluid side of this
ranking; naming follows how the arm was derived, from mouse iTreg 37
versus 39 °C contrasts, and the reading stays correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_unbiased_enrichment_viz.py` | `build_figure` | `thresholds.gsea_fdr = 0.05; gsea_min_size = 5; gsea_max_size = 500` | `03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/_overview/progeny_paired_forest.png

Tested one donor at a time, 10 of the fourteen PROGENy footprints
separate the JIA synovial Treg pool from paired blood at FDR < 0.05.
EGFR is the largest at +10.9 (FDR 2.7e-04), followed by NFkB at +8.3
(FDR 9.9e-04) and Hypoxia at +7.4 (FDR 1.2e-03) across all six paired
donors, while TNFa sits lower on the synovial-fluid side at -1.96
(FDR 0.028), so the footprints that separate split across both
directions.

**How to read:** One row per PROGENy pathway, rows ordered by the Treg paired
difference so the ordering belongs to one population. Inside a row
the three sorted populations are offset vertically, Treg above, Tcon
centre, CD8 below, so their markers and intervals stay separately
readable; colour repeats the same key. Horizontal position is the
mean within-donor difference in activity, synovial fluid minus that
donor's own blood, and the bar through it is the 95% interval of the
paired t-test. A filled marker reaches FDR < 0.05 in that population
and an open marker sits above that threshold, so the distinction
survives a greyscale print. Right of zero the footprint is higher in
synovial fluid, left of it higher in blood. Two limits on the
reading. The activity score is computed on expression centred within
a population, so a difference is comparable between pathways of the
same population and carries no meaning compared across populations; a
Treg point sitting further right than a CD8 point on the same row is
not a between-population effect size. And six pairs is a small n, so
an open marker here leaves the question open rather than settling it.
The six donors are the same count in every population but not the
same six people: Treg pairs JIA_patient_3 and skips 5, Tcon and CD8
do the reverse. Correlative; a footprint is inferred from target-gene
expression.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_progeny_paired_forest_viz.R` | `main` | `progeny.organism=Human; progeny.top=500; progeny.minsize=5; gsea_fdr=0.05; population_offset=0.24 row units` | `03_results/14_unbiased_enrichment/tables/progeny_sf_vs_pb.csv` |

## figures/by_contrast/ (per-database GSEA battery)

The full browse surface for the unbiased sweep: 108 panels across 3
sorted populations and 11 gene-set collections, one directory per
(population, collection) cell. Read against the whole family, the
Treg contrast carries 1,490 pooled-significant sets out of 11,236
tests, so the mouse-derived WT_heat_up arm (NES +2.5916, pooled FDR
4e-12, 119 of 199 genes in the ranked list) and HALLMARK_HYPOXIA (NES
+2.2563, pooled FDR 6e-08, 139 genes) both sit inside a broad
co-enrichment rather than standing alone. Three collections are too
small for the full battery and carry fewer panels: mouse_projection
(3 set(s), 2 panel type(s)); project_frozen (1 set(s), 1 panel
type(s)); sting_axes (2 set(s), 1 panel type(s)). The omission is a
redundancy judgement recorded per collection above, and no statistic
is withheld: every set of every collection is in gsea_all.csv.

**How to read:** LAYOUT.
figures/by_contrast/&lt;population&gt;/&lt;COLLECTION&gt;/{dotplot,facet,barplot,running_sum}.{pdf,png},
with the rows behind each panel in the mirrored path under
tables/by_contrast/. Population directories are Treg, Tcon and CD8,
and the contrast inside every one of them is the same donor-paired
synovial fluid versus peripheral blood comparison, published in each
CSV as SF_vs_PB_&lt;population&gt;. WHERE TO START. The three
Hallmark dotplots, one per population: fifty named programs on a
top-20 axis with both hypoxia and interferon among them. GLYPHS,
shared by every cell of the battery. dotplot: x = GeneRatio
(leading-edge genes divided by set size), point size = -log10(pooled
adjusted p), fill = NES with orange #B35806 positive and blue #2166AC
negative and the fill squished at plus or minus 3.5, black outline =
pooled FDR < 0.05. The dotplot SELECTS by adjusted p and ORDERS its
y-axis by GeneRatio descending, so vertical position there is a
gene-ratio ranking. facet: the same dotplot split into an NES > 0
block and an NES < 0 block. barplot: NES bars from zero, ordered by
NES. running_sum: three stacked panels, the running enrichment score
with its leading-edge peak on top, gene-hit ticks at each member's
rank in the middle, and the ranked moderated t at the bottom, with
the score clamped to [-1, 1] so curves stay comparable between
collections. SIGN. NES > 0 means the set's genes concentrate on the
synovial-fluid side of the ranking and NES < 0 on the paired-blood
side. ADJUSTED P. Every panel uses the Benjamini-Hochberg correction
across the whole family of tests asked of one population's ranked
list, which is stricter than a single-collection correction; the
per-collection value travels in each same-stem CSV under
padj_in_database. RANKING. The four panel types rank by different
metrics, adjusted p for dotplot and facet and |NES| for barplot and
running_sum, so read an absence against the rule named in that
panel's own subtitle before reading it as a null. This is a browse
surface, wide on purpose and privileging no set, and the claim spine
stays the donor-pseudobulk effect sizes. Correlative: enrichment
describes where a set's genes sit in a ranking. Claim tier: L3
(enrichment statistics), and no row of this battery reaches an
effect-size accumulator.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `emit_cell` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/*.rds + 03_results/14_unbiased_enrichment/tables/{gsea_all,geneset_manifest}.csv` |

## figures/_overview/named_sets_in_sweep.png

Scored with no favourite against all 11,236 sets tested in the JIA
Treg contrast, HALLMARK_HYPOXIA reaches NES +2.26 at pooled FDR
5.8e-08, rank 92 of 11,236 by pooled FDR, while the best-placed of
the six cGAS-STING sets reaches pooled FDR 0.195 at rank 2,934.
sting_specific_up does reach pooled FDR 0.018 in Tcon on 15 genes, so
the cGAS-STING reading turns on which population is read. Set size
tracks the outcome closely: 46.8% of Treg sets of 130 to 150 genes
reach pooled significance against 6.9% of sets of 10 to 22 genes, the
band five of the six cGAS-STING sets fall in, and 1,490 of 11,236
Treg tests are significant at all.

**How to read:** Columns are the three sorted populations on one shared row axis. Each
of the eleven upper rows is one named set, coloured by comparison
thread and ordered inside a thread by its Treg NES. Below the dashed
separator the bottom row is every set tested in that population, one
grey point each, which is the distribution a marker is read against.
Horizontal position is NES clamped to plus or minus 3.5; right of
zero the set's genes concentrate on the synovial-fluid side of the
ranking. A filled marker reaches pooled FDR < 0.05, an open marker
sits above it. Marker area is genes reaching the ranked list, so a
large NES on a small marker rests on few genes. Grey text gives
pooled FDR and rank within that population's whole sweep, ranked on
pooled FDR alone; by NES the same sets order differently, and the two
orderings answer different questions. A cell reading 'not tested' had
fewer than the minimum five of its genes in that ranked list, so it
carries no result. Four bounds on the reading. This contrast moves
many programs at once: 1,490 of 11,236 tests reach pooled FDR < 0.05
in Treg, 2,165 of 11,459 in Tcon and 1,027 of 11,242 in CD8. Set size
drives that rate: in Treg 73 of 156 sets of 130 to 150 genes are
pooled-significant (46.8%) against 231 of 3,327 sets of 10 to 22
genes (6.9%); five of the six cGAS-STING sets sit in that smaller
band and the sixth at five genes, while HALLMARK_HYPOXIA carries 139
testable genes. KO_heat_up is drawn beside WT_heat_up because that
comparator reaches pooled FDR 1.1e-12 against 3.7e-12, ranks 29
against 31, and the two lists share 182 genes. And the cGAS-STING
family's own signs disagree: two regulation-of terms carry opposite
sign and the positive-regulation term runs negative. The selection,
its reason per row, and the two excluded substring matches are
committed in tables/sweep_named_sets.csv. Correlative throughout.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_sweep_named_sets_viz.R` | `main` | `gsea_min_size=5; gsea_max_size=500; gsea_fdr=0.05; nes_cap=3.5; padj_pooled_method=BH; set_selection=tables/sweep_named_sets.csv` | `03_results/14_unbiased_enrichment/tables/{gsea_all,sweep_named_sets,sweep_named_sets_stats,sweep_setsize_baseline}.csv` |


## figures/by_contrast/&lt;population&gt;/GO_BP/*.png

GO_BP GSEA of the donor-pseudobulk synovial-fluid-versus-paired-blood
contrast in 3 sorted populations, drawn with the RNAseq-toolkit
plotters on the cached gseaResult from 03_results/objects/14_gsea/
with the adjusted p re-keyed to the sweep-wide pooled correction
published in gsea_all.csv. Gene Ontology biological process terms
from MSigDB C5, the largest collection in the battery. The gap
between a per-database adjusted p and the pooled adjusted p is
smallest here, because this collection supplies about half of the
pooled family. All four panels are drawn for this collection (5902
sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/GO_BP/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__GO_BP.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/GO_MF/*.png

GO_MF GSEA of the donor-pseudobulk synovial-fluid-versus-paired-blood
contrast in 3 sorted populations, drawn with the RNAseq-toolkit
plotters on the cached gseaResult from 03_results/objects/14_gsea/
with the adjusted p re-keyed to the sweep-wide pooled correction
published in gsea_all.csv. Gene Ontology molecular function terms
from MSigDB C5. All four panels are drawn for this collection (1407
sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/GO_MF/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__GO_MF.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/Reactome/*.png

Reactome GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. Reactome
canonical pathways from MSigDB C2. Set names here are the longest in
the battery and are wrapped rather than shortened, so the label
column is wide. All four panels are drawn for this collection (1600
sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/Reactome/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__Reactome.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/GO_CC/*.png

GO_CC GSEA of the donor-pseudobulk synovial-fluid-versus-paired-blood
contrast in 3 sorted populations, drawn with the RNAseq-toolkit
plotters on the cached gseaResult from 03_results/objects/14_gsea/
with the adjusted p re-keyed to the sweep-wide pooled correction
published in gsea_all.csv. Gene Ontology cellular component terms
from MSigDB C5. All four panels are drawn for this collection (861
sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/GO_CC/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__GO_CC.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/WikiPathways/*.png

WikiPathways GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. WikiPathways
canonical pathways from MSigDB C2. All four panels are drawn for this
collection (848 sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/WikiPathways/ listing the rows
it drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__WikiPathways.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/TF_Targets/*.png

TF_Targets GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. CollecTRI
regulons, one unsigned gene set per transcription factor with
activating and repressing targets pooled, so a set is that factor's
transcriptional neighbourhood. A set enriching says the factor's
targets move with one side of the contrast; it is a statement about
target-gene expression and carries no measurement of the factor's
activity. All four panels are drawn for this collection (600 sets in
the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/TF_Targets/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__TF_Targets.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/KEGG/*.png

KEGG GSEA of the donor-pseudobulk synovial-fluid-versus-paired-blood
contrast in 3 sorted populations, drawn with the RNAseq-toolkit
plotters on the cached gseaResult from 03_results/objects/14_gsea/
with the adjusted p re-keyed to the sweep-wide pooled correction
published in gsea_all.csv. KEGG canonical pathways from MSigDB C2.
The configuration asks for CP:KEGG and msigdbr 26 resolves that to
the CP:KEGG_LEGACY subcollection; the resolved name is recorded in
geneset_manifest.csv. All four panels are drawn for this collection
(185 sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/KEGG/ listing the rows it drew,
in draw order, with the rule that picked them and the per-collection
adjusted p under padj_in_database. A set enriching says its gene
content moves with one side of this contrast. Correlative. Claim
tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__KEGG.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/Hallmark/*.png

Hallmark GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. Fifty broad
MSigDB Hallmark programs, the collection a bench reader can name.
HALLMARK_HYPOXIA and both interferon-response sets are members, so
this is the one panel of the battery where the hypoxia and interferon
readings sit on a single axis with everything else. All four panels
are drawn for this collection (50 sets in the pooled family).

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/Hallmark/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. Worked example of
the two rankings diverging: in the Treg Hallmark cell
HALLMARK_HYPOXIA (NES +2.2563, pooled FDR 6e-08, 139 genes in the
ranked list) is 12 by |NES| and 10 by pooled adjusted p, so it sits
inside the dotplot's top 20 and outside the running sum's top 5. An
absence from one panel is a statement about that panel's ranking
metric. A set enriching says its gene content moves with one side of
this contrast. Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__Hallmark.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/mouse_projection/*.png

mouse_projection GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. The three
mouse-derived up arms projected onto human symbols: WT_heat_up (199
genes, 119 in the Treg ranked list), KO_heat_up and Interaction_up.
They are ordinary members of the sweep with no privilege, and
WT_heat_up doubles as the reproduction check against the published
targeted result on the same ranked list. Panels drawn for this
collection: dotplot, running_sum. The collection offers 3 set(s) to
the pooled family, and over a collection that small a dotplot, a
direction split and a significant-only barplot are the same few
points drawn several ways, so the redundant panels are left out. The
running sum is kept at every size because it is the only panel that
shows where in the ranking a set's genes sit.

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/mouse_projection/ listing the
rows it drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__mouse_projection.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/sting_axes/*.png

sting_axes GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. The two
frozen axes from the SAVI positive-control compartment,
sting_specific_up and ifn_only_up, which separate STING-attributable
content from generic type-I interferon content. Panels drawn for this
collection: running_sum. The collection offers 2 set(s) to the pooled
family, and over a collection that small a dotplot, a direction split
and a significant-only barplot are the same few points drawn several
ways, so the redundant panels are left out. The running sum is kept
at every size because it is the only panel that shows where in the
ranking a set's genes sit.

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/sting_axes/ listing the rows it
drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__sting_axes.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/by_contrast/&lt;population&gt;/project_frozen/*.png

project_frozen GSEA of the donor-pseudobulk
synovial-fluid-versus-paired-blood contrast in 3 sorted populations,
drawn with the RNAseq-toolkit plotters on the cached gseaResult from
03_results/objects/14_gsea/ with the adjusted p re-keyed to the
sweep-wide pooled correction published in gsea_all.csv. The frozen
curated lists this compartment owns. Six of the seven are re-pins of
MSigDB Hallmark sets with identical gene content, already drawn in
the Hallmark panel with the same statistics, and geneset_manifest.csv
records them as n_sets_aliased_out_of_pooling = 6. They are excluded
here, so HALLMARK_HYPOXIA appears once in this battery and HSR_core
is the only set drawn under this name. HSR_core is the curated
heat-shock-response lens, held independent of the mouse anchor and
general to proteotoxic stress. Panels drawn for this collection:
running_sum. The collection offers 1 set(s) to the pooled family, and
over a collection that small a dotplot, a direction split and a
significant-only barplot are the same few points drawn several ways,
so the redundant panels are left out. The running sum is kept at
every size because it is the only panel that shows where in the
ranking a set's genes sit.

**How to read:** SELECTION RULES, which govern every absence: dotplot top 20 by pooled
adjusted p; facet top 10 per direction by pooled adjusted p; barplot
sets at pooled FDR < 0.05 only, then top 20 of those by |NES|;
running_sum top 5 by |NES|. Glyphs, the sign convention and the
pooled correction are described once in the `figures/by_contrast/
(per-database GSEA battery)` section of this README. Each panel
writes its own same-stem CSV under
tables/by_contrast/&lt;population&gt;/project_frozen/ listing the
rows it drew, in draw order, with the rule that picked them and the
per-collection adjusted p under padj_in_database. A set enriching
says its gene content moves with one side of this contrast.
Correlative. Claim tier: L3 (enrichment statistics).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14c_gsea_battery_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; thresholds.gsea_min_size=5; thresholds.gsea_max_size=500; figures.top_n=20; facet_top_n_per_direction=10; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/objects/14_gsea/{treg,tcon,cd8}__project_frozen.rds + 03_results/14_unbiased_enrichment/tables/gsea_all.csv` |

## figures/_overview/program_nes_by_cell_state.png

All four oxygen-named sets rise on the synovial-fluid side of the
paired contrast in Treg, and they spread across the panel:
HALLMARK_HYPOXIA reaches NES 2.2563 at pooled FDR 6e-08 on 139 genes
while GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS reaches 1.2172 at 0.334
on 120, so a reading taken from one of the four would be a reading the
other three do not reproduce. In the second group the two interferon
sets carry the strongest rows in the figure, ifn_only_up reaching NES
2.6057 at 5e-10 in CD8. sting_specific_up clears pooled FDR 0.05 in
Tcon at NES 1.8745 on 15 genes and sits above the threshold in Treg
(1.5187, 0.195) and CD8 (1.6700, 0.102), and the three pathway-
database STING terms stay between NES -0.84 and 1.37 on 10 to 24
genes.

**How to read:** The companion to arm_nes_by_cell_state under
03_results/14_unbiased_enrichment/, drawn to the same geometry so the
two can be laid side by side. One dot per gene set and cell state, at
the confirmatory tier: donor-level pseudobulk within frozen sort
labels, limma-voom moderated t, then pre-ranked fgsea. Rows are gene
sets in two labelled groups, ordered inside a group by descending Treg
NES, and inside a row the three cell states are offset vertically and
coloured, each with its own annotation line. The x position is the
normalised enrichment score for synovial fluid over paired blood, with
a vertical rule at zero. A filled dot clears the config FDR threshold
of 0.05 and an open dot sits above it. A cell reading not tested had
fewer genes in that population's ranked list than gsea_min_size, so
the sweep never scored it and the cell carries an absence rather than
a null. Read every score against the gene count beside it. In this
sweep, size alone moves the odds a long way: in Treg a set of 130 to
150 genes clears pooled FDR 0.05 in 46.8% of the 156 such sets tested,
while a set of 10 to 22 genes clears it in 6.9% of 3,327. The four
oxygen-named sets are four differently built sets for one named
biology, and their scores differ, so the group is read as a family and
a single member's score is read as that member's. The cGAS-STING group
carries a positive result in Tcon, so the group is read set by set.
The complete six-member cGAS-STING family of this sweep, including the
two regulation-of terms whose signs disagree, is drawn in
named_sets_in_sweep under the same stage. Temperature and hypoxia are
both imposed by the inflamed joint and stay entangled in cross-
sectional human data, so these rows describe what the niche contrast
contains and the question of what drives which stays open. A score of
this kind has no interval, so none is drawn, and the reading stays
correlative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_program_nes_by_cell_state_viz.py` | `build_figure` | `thresholds.gsea_fdr = 0.05; gsea_min_size = 5; gsea_max_size = 500; row order = descending Treg NES within a group` | `03_results/14_unbiased_enrichment/tables/_overview/named_sets_in_sweep.csv, 03_results/14_unbiased_enrichment/tables/gsea_all.csv, 03_results/14_unbiased_enrichment/tables/sweep_setsize_baseline.csv` |

