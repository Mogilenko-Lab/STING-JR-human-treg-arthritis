# 00_build — artifact captions

_**Abbreviations:** SF = synovial fluid (inflamed joint); PB = peripheral blood. The SF-vs-PB contrast is paired within each of the 7 JIA donors. Treg = CD4⁺CD127ˡᵒCD25⁺ regulatory; Tcon = CD4⁺CD25⁻ conventional; CD8 = CD8⁺CD45RO⁺ memory._

## figures/_overview/cells_per_gsm.png

All 7 donors contribute paired SF+PB Tregs; Tcon and CD8 lack a PB
sample for p3 (by design). SF-Treg p5 is the thinnest stratum.

**How to read:** Grouped bars = cells recovered per donor; orange = synovial fluid
(SF), blue = peripheral blood (PB); one facet per sorted population. A
missing PB bar (p3 in Tcon/CD8) is an intentionally-absent sample, not
a QC drop. Descriptive counts only — no claim tier.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/00_build_anndata_viz.py` | `main` | `design.populations = [CD4_Treg, CD4_Tcon, CD8]` | `03_results/00_build/tables/cells_per_gsm.csv` |
