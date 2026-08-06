# interactive/ — per-cell substrates for the explorer views

Six Parquet files, one per stage that produced an explorable per-cell table. Each holds one row
per cell with the coordinates, labels and scores that stage computed, flattened so a viewer can
load it without the AnnData object.

These are inputs to interactive views. The figures and tables under each `NN_<stage>/` directory
are the deliverables, and these files let someone re-plot the same cells themselves.

| File | Cells | What it carries |
|---|---|---|
| `01_qc_explore.parquet` | all barcodes before filtering | per-cell QC metrics and the keep/drop decision |
| `02_annotation_explore.parquet` | post-QC cells | frozen cell-state labels with the marker scores behind them |
| `05_gonogo_explore.parquet` | post-QC cells | per-cell signature scores |
| `08_harvest_readout.parquet` | post-QC cells | population-of-interest hook membership, carried as annotation |
| `16_narrative_embedding.parquet` | 99,915 | full-object embedding coordinates with AUCell scores |
| `17_treg_reembedding.parquet` | 27,175 | Treg-only embedding coordinates with the same score columns |

Written by `07_embedding.py`, `08_harvest_readout.py` and `16_narrative_scoring.py`.

## Two things to know before using them

An embedding is a map. Proximity on it reflects the neighbourhood graph that produced it. Biology
is tested on donor-level pseudobulk within the frozen labels.

The population-of-interest hooks in `08_harvest_readout.parquet` are a resource for generating
hypotheses. Any result worth reporting is re-derived in the donor-level analysis first.

## Reading them

```python
import pandas as pd
df = pd.read_parquet("03_results/interactive/17_treg_reembedding.parquet")
```

Total 47 MB. These are excluded when the results are packaged for sharing, since the per-stage
tables carry the same numbers in a form that opens in a spreadsheet.
