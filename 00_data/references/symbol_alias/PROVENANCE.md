# symbol_alias — reference symbols resolved into this compartment's vocabulary

## What this is

`symbol_alias_map.csv` is a lookup from the HGNC symbol a reference gene set **ships**
to the symbol this compartment's count matrix **carries**, with one row per candidate
pair and an explicit resolution on each. `symbol_alias_provenance.csv` records what the
map was built from, so a stale map is visible rather than assumed.

## Why it exists

GSE160097 ships 40 per-GSM CellRanger outputs, 33 of them in the v2 layout, and its
feature union is the 32,738-gene **hg19** reference (`var['genome'] == 'hg19'` in
`03_results/objects/00_build.h5ad`). The matrix is therefore frozen to that build's HGNC
vintage:

| the gene | this matrix | current HGNC |
|---|---|---|
| cGAS | `MB21D1` | `CGAS` |
| STING | `TMEM173` | `STING1` |
| MARCHF5 | `MARCH5` | `MARCHF5` |
| MRE11 | `MRE11A` | `MRE11` |

Reference collections ship current symbols and every match in this compartment is an
exact string match, so a renamed gene leaves a set silently — and the loss reads as
biological absence. In the Treg synovial-fluid-versus-paired-blood ranked list
`TMEM173` sits at rank 265 and `MB21D1` at 458 of 13,999, the top 3.3%. The two genes
that name the cGAS-STING axis are its two strongest members here and they are invisible
to all six of its gene sets.

## How it was built

`Rscript 02_analysis/scripts/00_symbol_alias_map.R`, from the compartment root.

- **Arbiter:** org.Hs.eg.db (version in the provenance row).
- **Resolved into:** `03_results/03_pseudobulk/tables/gene_symbols.csv`, the 21,740-symbol
  post-QC vocabulary the pseudobulk seam emits. Deliberately not a ranked list: the map is
  a property of the matrix, so a pair whose target the expression filter later drops is
  reported by a consumer's ledger as expression-filtered rather than as vocabulary loss.
- **Reference universe:** the 20,149 unique symbols of everything this compartment
  consumes — seven MSigDB collections, the CollecTRI regulon targets, the toolkit's human
  MitoPathways build, and 33 frozen one-symbol-per-line lists (the projected mouse arms in
  **both** directions, the Hallmark re-pins, the curated HSR and TCR lenses, the SAVI axes,
  the eTreg sets).
- **The one-to-one condition:** a reference symbol qualifies when it is absent from the
  vocabulary, resolves to exactly one Entrez id, and exactly one alias of that same Entrez
  id is present in the vocabulary.
- **The ownership guard:** a candidate that is the official symbol of any *other* Entrez id
  is rejected. Many retired symbols were reassigned to a different gene — `ACOD1` carries
  `CAD`, which now names carbamoyl-phosphate synthetase; `IL17F` carries `IL17A`; `TTR`
  carries `TTN`; `PGF` carries `PIGF`; `THPO` carries `TPO`. Accepting one attaches one
  gene's expression to another gene's set membership. 120 candidates are refused this way
  and every refusal is a row, so the guard is auditable rather than implied.

## How to read `symbol_alias_map.csv`

One row per **candidate pair**. `resolution` is the decision:

| resolution | meaning | applied? |
|---|---|---|
| `accepted` | a clean 1:1 nomenclature update | yes |
| `flagged_for_review` | passes the guard, but is not a clean rename; withheld by human decision | **no** |
| `rejected_symbol_belongs_to_another_gene` | the candidate is another gene's official symbol | no |
| `rejected_multiple_aliases_in_vocabulary` | two aliases of the gene are in the vocabulary, so there is no unique target; both are named in `matrix_symbol`, `/`-joined | no |
| `rejected_reference_symbol_ambiguous_in_org_db` | the reference symbol carries more than one Entrez id | no |

Only `accepted` is ever applied. `n_aliases_in_vocabulary` is 1 for every accepted row by
construction. Reference symbols with no candidate at all are **not** rows — their counts
are in the provenance file, because carrying thousands of them would defeat the review
step this asset most needs: reading the map cold and recognising every accepted pair as a
nomenclature update.

`MIR4435-2HG -> MIR4435-1HG` is the one `flagged_for_review` pair. The two are distinct
lncRNA loci in most annotations; `-1HG` was reclassified and org.Hs.eg.db carries it as an
alias, but a lncRNA locus merge is not the same class of event as a protein-coding rename,
and this pair lands in `sting_specific_up`, which is a claim surface. It is excluded from
every tested set, which is why that set reaches 13 rather than 14 members in the Treg
ranked list. The exclusion is declared in `analysis_config.yaml` under
`symbol_alias.flagged_for_review`, and `00_symbol_alias_map.R` stops if a flagged pair is
no longer produced as a candidate — a stale exclusion must not sit silently inert.

## Consumers

Read through the helpers, never by re-deriving:

- **R** — `02_analysis/helpers/symbol_alias.R` (`build_alias_map`, `resolve_sets`,
  `symbol_ledger`, `accepted_pairs`).
- **Python** — `02_analysis/helpers/geneset_utils.py` (`load_alias_map`,
  `resolve_symbols`). Python cannot reach org.Hs.eg.db, which is why the map is a
  committed CSV both languages read rather than a live lookup.

Resolution runs **one way only**: the reference symbol is newer and the matrix symbol is
older, so a reference symbol is resolved *down* into the vocabulary the data carries.

## Rebuild triggers

Regenerate, and read the diff, when any of these moves:

1. `03_results/03_pseudobulk/tables/gene_symbols.csv` (its SHA-256 is in the provenance row).
2. The org.Hs.eg.db release.
3. A collection added to `unbiased_enrichment:` or a new frozen list under the directories
   named in `symbol_alias.extra_reference_dirs`.

Reporting rule, from the umbrella guardrail: matching is reported as a ledger with a
bucket per cause — matched, matched-via-alias, expression-filtered, below-detection,
absent-from-reference — and never as a recovery fraction or a pass/fail floor.
