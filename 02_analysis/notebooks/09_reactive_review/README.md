# 09_reactive_review — JIA synovial-Treg reactive review (marimo)

A single reactive [marimo](https://marimo.io) app that re-reads the JIA synovial-fluid
Treg result end to end and adds the per-cell hypoxia / UPR readouts. It mirrors the
recorded synovial-Treg review narrative and its tables, and adds a live dual-embedding
explorer. **The app only visualizes and narrates frozen tables — it never computes.**

## How to run

Live (reactive, in-browser):

```bash
python3 -m marimo run 02_analysis/notebooks/09_reactive_review/reactive_review.py \
  --headless --host 0.0.0.0 --port 2731 --no-token
# then open http://<host>:2731/  (lasso needs the live kernel)
```

Edit mode (to change cells):

```bash
python3 -m marimo edit 02_analysis/notebooks/09_reactive_review/reactive_review.py
```

Export a static snapshot (no live reactivity — the lasso panel is inert):

```bash
python3 -m marimo export html --no-include-code \
  02_analysis/notebooks/09_reactive_review/reactive_review.py \
  -o 02_analysis/notebooks/09_reactive_review/reactive_review.html
```

The `__marimo__/` cache and the exported `reactive_review.html` are gitignored
(regenerable from the `.py`).

## Frozen tables it reads

| Table | Role |
|---|---|
| `03_results/interactive/08_harvest_readout.parquet` | per-cell x/y, cell state, tissue, donor, %mt, hypoxia/UPR AUCell, WT_heat/eTreg/HSP scores — the embedding + readout substrate |
| `03_results/05_scoring/tables/gsea_pseudobulk_{treg,tcon,cd8}.csv` | primary donor-pseudobulk WT_heat NES |
| `03_results/master/effect_sizes_treg_arthritis.csv` (`effect_metric == percell_auc_smd`) | secondary per-cell AUCell SF-vs-PB SMD |
| `03_results/08_harvest_readout/tables/harvest_readout_summary.csv` | SF-vs-PB hypoxia / UPR readout summary by cell state |
| `03_results/07_embedding/tables/{or_union_membership,hook_per_lineage_summary}.csv` | current OR-gated drafted-subset context |
| `02_analysis/config/analysis_config.yaml::decisions.go_no_go` | decision echoed verbatim |

## Panels (top to bottom)

1. **Birds-eye framing** — the mouse→JIA question and the headline: the mouse `WT_heat`
   up-program enriches in synovial-fluid T cells versus paired blood, pan-T and not
   Treg-preferential. Correlative.
2. **Primary evidence** — `WT_heat` NES heatmap + table (Treg / Tcon / CD8), donor
   pseudobulk. NES ≈ 2.51 / 2.57 / 2.05, all FDR ≪ 1e-6. The primary tier.
3. **Secondary corroboration** — per-cell `WT_heat_up` AUCell SF-vs-PB donor SMD table.
   Corroborative tier, on its own scale.
4. **Reactive dual embedding** — two synced square UMAPs of one filtered+sampled cell
   set (caption above, lasso cluster-characteristics summary directly below). Left =
   cell-state reference; right = AUCell-only readout dropdown (WT_heat / eTreg / HSP /
   HALLMARK_HYPOXIA / HALLMARK_UPR / %mt). Lasso → composition + mean readout scores
   selected-vs-rest.
5. **Hypoxia + UPR readouts** — SF-vs-PB SMD by cell state + the hypoxia high-pocket
   table. Surfaces the synovial-Treg-leading hypoxia pocket (29.8% SF vs 1.7% PB Tregs
   in the top-decile pocket) and the flat UPR. Correlative readout, no HIF-causality claim.
6. **Harvest context (revisable)** — the current OR-gated drafted subset (≈ 33% of the
   compartment) from the embedding stage, framed as open to revision in light of the
   hypoxia readout.
7. **Decision** — the recorded continue call (GO, pan-T) plus an explicit note that the
   harvest strategy stays revisitable and is decided separately. Echoes
   `decisions.go_no_go` verbatim.

## Evidence-tier note

The **primary** tier is the donor-level pseudobulk `WT_heat` NES (panel 2) — the only
tier that carries the claim. The **secondary annotation** tier is everything per-cell:
the AUCell SMD (panel 3), the embedding readouts (panel 4), and the hypoxia / UPR
surfaces (panel 5). The two tiers are kept on separate scales and **never pooled**. The
hypoxia surface is a readout of the tissue state, read correlatively — not evidence that
HIF or fever drives it. The UMAP and the drafted subset are maps, never evidence; the
harvest is a current, revisable draft and the mouse `WT_heat` score annotates it but is
never a selection predicate.
