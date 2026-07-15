# heat_leadingedge_taxonomy — provenance

`leadingedge_gene_taxonomy.csv` classifies the union of the SF-vs-PB `WT_heat_up` fgsea
leading-edge genes (66 genes across Treg/Tcon/CD8, stage `09_heat_hypoxia`) into biological
programs: `heat_shock_proteostasis`, `hypoxia_HIF`, `immediate_early_stress`,
`effector_activation`, `other`.

Produced 2026-07-14 by the `agy` (Gemini 3.1 Pro) CLI via the `delegate-cli` skill, prompted
to classify each gene by its known regulation (HSF1 / HIF1α / AP-1-immediate-early /
TCR-activation / effector-cytokine) with a one-line rationale (carried in the `evidence`
column). Genes in MSigDB HALLMARK_HYPOXIA were defaulted to `hypoxia_HIF`. Frozen here as a
committed reference so the leading-edge composition tally in `09_heat_hypoxia.py` is
reproducible. Regenerable by re-running the classification (agy) on the leading-edge union.

This is an annotation-tier interpretive aid, not a confirmatory statistic. The full agy output
(table + tally + synthesis) is retained in the compartment decision log.
