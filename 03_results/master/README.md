# master/ — the accumulator tables

Cross-stage tables that gather one row per result from the individual stages. A stage computes.
These files collect. Read them for every effect in one place.

## effect_sizes_treg_arthritis.csv · master_effect_sizes.csv

**12 rows. This compartment's confirmatory results, each with a confidence interval.**

One row per (signature × cell state × contrast). The estimate is a donor-level effect: each
donor contributes one value per cell state, so the sample size is the number of donors.
`n_donors` and `n_cells` both appear, and `n_donors` is the one that sets the power.

`evidence_tier` says how far a row may be taken. `signoff_state` records whether a human has
accepted it. A row carries a claim when both permit it.

`master_effect_sizes.csv` is the same 12 rows with `stage` and `database` added, so a row can
be traced back to the stage that produced it.

| Column | Meaning |
|---|---|
| `signature` | the gene set scored, named for how it was derived |
| `cell_state` | the frozen label the score was computed within |
| `contrast` | synovial fluid against paired blood, within donor |
| `estimate`, `se`, `ci_low`, `ci_high` | the effect and its interval |
| `padj` | multiple-testing adjusted p-value |
| `evidence_tier` | how far the row may be read |
| `signoff_state` | whether a human has accepted it |

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures.py` | `write_effect_sizes` | `scoring.effect_metric` | `03_results/05_scoring/tables/` |

## master_gsea_pseudobulk.csv

**48 rows. Gene-set enrichment across the donor-level pseudobulk contrasts.**

One row per (gene set × contrast × direction). `nes` is the normalised enrichment score, positive
when the set sits toward the synovial-fluid end of the ranking. `set_size` counts the genes of
that set present in the ranked list, which is smaller than the curated set: a gene absent from
the ranking was either undetected or removed by the expression filter.

`core_enrichment` lists the genes carrying the score.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/05_score_signatures.py` | `accumulate_gsea` | `thresholds.gsea_min_size`, `thresholds.gsea_max_size` | `03_results/05_scoring/tables/gsea_pseudobulk_*.csv` |

## Reading a set size beside a statistic

A statistic and its set size come from the same run. Symbol-alias resolution has since raised
some set sizes, and those enrichments are re-read when the affected stages regenerate. Where a
caption gives both, it names the run each belongs to.
