> **Submodule of the STING-JR super-repo.** For cross-project context and shared conventions, read the umbrella `../AGENTS.md` first.

# AGENTS.md — human_treg_arthritis

**Created:** 2026-06-29 · **Type:** analysis · **Species:** Homo sapiens (GRCh38)

This file is the single source of truth for AI agent behavior in this project.
`CLAUDE.md` imports it via `@AGENTS.md`. Edit here, not there.

---

## Project context

The canonical scientific framing — question, datasets, hypotheses, goals — lives in
`docs/_internal/scientific-context.md`. Read it before acting. Session handoffs are written
to `docs/_internal/sessions/YYYY-MM-DD_<slug>.md`.

## Active role

_No role activated yet. Run `sciagent activate <role>` to populate this section._

Run `sciagent status` to see the current stack, effective skills, agents, and commands.

## Critical rules

1. **Config, not hardcoding.** Paths and parameters come from
   `02_analysis/config/analysis_config.yaml`. Never inline file paths or thresholds.
2. **Normalize, then visualize.** Compute scores/embeddings and checkpoint them under
   `03_results/objects/` before any visualization step.
3. **Read-only input data.** `00_data/` is read-only. Outputs go to `03_results/`.
4. **Cache expensive operations.** Anything over a minute goes through a checkpoint object.
5. **Pseudobulk DE runs in R, across an explicit seam.** Aggregation is Python
   (`NNa_pseudobulk_export.py`, counts only, no statistics); the model is R
   (`NNb_pseudobulk_de.R`, edgeR/limma-voom). The export **must** emit `gene_map.csv`
   (`ensembl_id,gene_symbol`) and the R side must join on it before ranking — count matrices are
   keyed by Ensembl id while every reference gene set matches on HGNC symbol, so skipping the map
   yields ranked lists that intersect the references at ~zero and fail **silently**, looking like a
   biological null. Collapse to one row per symbol on max `|t|`, rank on the sign-preserving
   moderated `t`. Full contract and rationale: umbrella `../AGENTS.md`.

## Directory structure

```
human_treg_arthritis/
├── 00_data/{raw,processed,references}   # READ-ONLY input data
├── 01_modules/                          # toolkits (submodules)
├── 02_analysis/{config,helpers,scripts,notebooks}
├── 03_results/
│   ├── objects/                         # checkpoint state (.h5ad, .rds) — NOT a deliverable
│   ├── 01_qc/{tables,figures}/          # QC-phase deliverables
│   ├── 02_eda/{tables,figures}/         # EDA-phase deliverables
│   ├── master/                          # cross-stage accumulator tables
│   ├── interactive/                     # standalone HTML dashboards
│   └── _scratch/                        # throwaway plots/tables
├── docs/                                # public-facing: plan, guides
└── logs/
```

Each `NN_<slug>/` stage holds the tables **and** figures for that phase. Stage ids are declared
in `analysis_config.yaml:stages` (two-digit prefix, lowercase snake_case); append a stage there
before writing to it. Checkpoints live in `objects/`, cross-stage tables in `master/`.

## Artifact captions

Every file in `03_results/<stage>/` needs a caption entry in that stage's `README.md`, written
as the artifact is produced: `## <filename>`, a one-sentence scientific **finding** (what it
shows, not what it is), then a `Script | Function | Config | Input` table. Function is the exact
name read from the script (never guessed); Config is the active `key = value`. Never cite
`docs/_internal/` in these READMEs.

## Figures

Viz scripts source `02_analysis/helpers/figure_style.{R,py}` (the project shim, which loads the
toolkit-symlinked `figure-style` lib) and save figures via `save_overview()` or
`save_figure()` — never inline `ggsave`/`plt.savefig` or `theme()`. See the
`figure-style` skill for the full contract (one legible tier, dual-format PDF + PNG from one plot object, font floors, captioning).

## Craft conventions

Standing craft conventions for this repo are toolkit-managed in the `<!-- SCIAGENT:CRAFT -->`
block rendered by `sciagent activate`. That block is authoritative — do not hand-edit it.

---

## Documentation namespace

`docs/_internal/` is agent-facing (not the published deliverable). Route every output here:

| Output kind | Directory |
|-------------|-----------|
| session handoff | `docs/_internal/sessions/` |
| research note | `docs/_internal/research/` |
| decision log | `docs/_internal/reasoning/` |
| scientific framing | `docs/_internal/scientific-context.md` |
| public phased plan | `docs/_internal/plans/{date-slug}/` |

Naming conventions: see `docs/_internal/README.md`.

**One-way reference rule:** public `docs/` must NOT reference `_internal/` paths.

**Never-grow-this-file rule:** if `AGENTS.md` exceeds ~150 lines, split content into
referenced files under `docs/_internal/`.

---

## Notes

**Raw mounts shadow their `.gitkeep` — do not commit the deletion.** Each
`00_data/<accession>/raw/` is a read-only bind mount of the staged bytes on
`/data2`, so the tracked `00_data/**/raw/.gitkeep` placeholder disappears from the
working tree whenever the container is up and git reports it as deleted. It is not
deleted: the placeholder must stay in HEAD or a fresh clone has no `raw/` directory
to mount into, and the mount is read-only so it cannot be restored in place.
Silence the phantom deletion once per clone, alongside `bash hooks/install.sh`:

```bash
git update-index --skip-worktree 00_data/*/raw/.gitkeep
```

This is local index state, never committed. The other compartments carry the same
mount and the same fix.

<!-- BEGIN SCIAGENT:ROLES v1 hash=cef2142357f46125eef5599c5524670bd467f933 -->
# Active roles

Stack (in order, last-wins on name collisions):
1. **base** — `roles/base.yaml` — Default bioinformatics analysis role with full agent suite
2. **pathway-signature** — `roles/pathway-signature.yaml` — Pathway/TF/signature functional interpretation — GSEA + decoupleR + CoReSh + pathway-explorer  *(overlay)*

## Skills (effective)
- `anndata` — (pathway-signature) [shadows base]
- `scanpy` — (pathway-signature) [shadows base]
- `single-cell-rna-qc` — (base)
- `anndatar-seurat-scanpy-conversion` — (pathway-signature) [shadows base]
- `bulk-rnaseq-gsea` — (pathway-signature) [shadows base]
- `bulk-rnaseq-activity-inference` — (pathway-signature) [shadows base]
- `bulk-rnaseq-pathway-explorer` — (pathway-signature) [shadows base]
- `gatom-metabolomic-predictions` — (pathway-signature) [shadows base]
- `coresh-signature-search` — (pathway-signature) [shadows base]
- `starsolo-spliced-unspliced` — (base)
- `rna-velocity-trajectory` — (base)
- `genenmf-metaprogram-discovery` — (base)
- `pycistopic-atac-topic-modeling` — (base)
- `crescendo-scatac-cre-analysis` — (base)
- `chromvar-motif-accessibility` — (base)
- `tf-footprint-differential-analysis` — (base)
- `scenic-grn-inference` — (base)
- `scvi-framework` — (base)
- `scvi-basic` — (base)
- `scvi-scanvi` — (base)
- `scvi-multivi` — (base)
- `scvi-peakvi` — (base)
- `scvi-mrvi` — (base)
- `scvi-contrastivevi` — (base)
- `scvi-linearscvi` — (base)
- `scvi-lda` — (base)
- `scvi-hub-models` — (base)
- `scvi-scarches-reference-mapping` — (base)
- `treearches-hierarchy-learning` — (base)
- `scglue-unpaired-multiomics-integration` — (base)
- `figure-style` — (base)
- `scrna-pipeline-conventions` — (base)
- `cellranger-multi-to-anndata` — (base)
- `scrna-cxg-host` — (base)
- `consensus-nmf-multirun` — (base)
- `skill-creator` — (base)
- `reasoning-trace` — (base)

## Skills (inherited via requires:)
- `tobias-footprint-bindetect` — (via `tf-footprint-differential-analysis`)
- `hint-atac-differential-footprint` — (via `tf-footprint-differential-analysis`)
- `signac-footprint-visualization` — (via `tf-footprint-differential-analysis`)

## Sub-agents (effective, Claude-only)
- `docs-librarian` — (pathway-signature) [shadows base]
- `bio-interpreter` — (pathway-signature) [shadows base]
- `insight-explorer` — (pathway-signature) [shadows base]
- `captions` — (pathway-signature) [shadows base]
- `figure-audit` — (base)
- `doc-curator` — (pathway-signature) [shadows base]
- `code-reviewer` — (pathway-signature) [shadows base]
- `handoff` — (pathway-signature) [shadows base]

## Slash commands (effective, Claude-only)
- `/commit` — (pathway-signature) [shadows base]
- `/decompose` — (pathway-signature) [shadows base]
- `/implement` — (pathway-signature) [shadows base]
<!-- END SCIAGENT:ROLES -->

<!-- BEGIN SCIAGENT:CRAFT v1 hash=2399c4db2eb105a2de6c58e227fc7983dcf81fcc -->
# Craft standards

Toolkit-managed standing conventions for this repo — do not hand-edit.

- **Figures** — legible both shrunk in a journal column and projected to the back of a room: bigger, fewer, bolder (one legible tier, base >= 14pt; emit a vector PDF + raster PNG from one plot object). Style only via the project theme entry point (no inline `theme()`/`ggsave(width=)`/raw hex); cap to top-N; never truncate axis labels; disambiguate glyphs; prefer the residualized channel.
- **Results** — every artifact under `03_results/<stage>/{figures,tables}/` with `by_contrast/<c>/` + `_overview/`; a figure's source table is its same-stem neighbor. Compute never plots; viz never computes.
- **README** — a task is unfinished until the sibling `README.md` captions every `03_results/` file you create/edit/delete, including *how to read* it (glyphs, sign convention, Δρ, claim tier).
- **Planning** — plans in `docs/_internal/plans/{date-slug}/` as `00_INDEX.md` + `NN_<slug>.md`; one phase == one script == one implementer (~35% of context); review every 3 phases that code runs and produces its artifacts.
- **Reproducibility** — no ephemeral scripts: every `03_results/` artifact reproducible from a committed `02_analysis/scripts/NN_*`; log non-trivial decisions to `docs/_internal/reasoning/` before proceeding (`_scratch/` is the only sanctioned throwaway zone).
- **Memory & traceability** — durable project memory lives only in tracked locations (`AGENTS.md`, `docs/_internal/`, committed scripts). Ephemeral/auto memory of vendored coding harnesses (e.g. Claude Code auto-memory, untracked `~/.claude` state written outside the repo) is prohibited as project state — `autoMemoryEnabled` is off; record decisions in `docs/_internal/reasoning/` and session handoffs in `docs/_internal/sessions/` so context stays reproducible, traceable, and observable.
<!-- END SCIAGENT:CRAFT -->
