# Breakpoint notebooks

Two kinds, one per human-review inflection point.

## Static review notebooks (R · rendered to HTML)

`NN_<topic>_review/NN_<topic>_review.qmd` — the recorded steering gate. Reads the stage's
published CSV tables + embeds the figures, opens with the narrative questions, and terminates
in a recorded `analysis_config.yaml::decisions.*` claim + branch. Render with:

```bash
Rscript 02_analysis/notebooks/render.R 02_analysis/notebooks/05_gonogo_review/05_gonogo_review.qmd html
```

- `01_qc_review` → `decisions.qc` · `02_annotation_review` → `decisions.treg_gate`
- `05_gonogo_review` → `decisions.go_no_go`

## Interactive explorers (Python · jscatter · live kernel)

`NN_<topic>_explore/NN_<topic>_explore.qmd` — GPU-accelerated interactive UMAP
([jupyter-scatter](https://github.com/flekschas/jupyter-scatter): pan / zoom / brush / hover /
linked selection) for *operating* at each inflection point. **These need a live kernel** — open
in VS Code (Quarto + Jupyter extensions) or JupyterLab with the **Python (STING base)** kernel
and run the cells. A static HTML render only captures a snapshot.

They load compact tables from `03_results/interactive/*.parquet` (coords + obs + a few
gene/score columns — the live kernel never opens the multi-GB checkpoints). Build/refresh them
whenever the checkpoints change:

```bash
python 02_analysis/scripts/export_explorers.py
```

| Explorer | Loads | Brush to answer |
|---|---|---|
| `01_qc_explore` | `01_qc_explore.parquet` | sort separation, contamination, where the high-mito (stressed) cells sit |
| `02_annotation_explore` | `02_annotation_explore.parquet` | are `sort_consistent=False` cells still FOXP3⁺ (activated) or true mis-sorts |
| `05_gonogo_explore` | `05_gonogo_explore.parquet` | is the per-cell WT_heat score Treg-preferential or pan-T; SF vs PB shift |

**Setup (one-time, already done in this environment):** `jupyter-scatter` + the Jupyter kernel
stack are installed in the base env and a `Python (STING base)` kernelspec is registered. If the
container is rebuilt, re-run:

```bash
pip install jupyter-scatter ipykernel jupyterlab
python -m ipykernel install --user --name sting-base --display-name "Python (STING base)"
```
