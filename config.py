"""
config.py — the single source of paths + parameters for the Python pipeline.
=============================================================================
Per `scrna-pipeline-conventions`: nothing downstream hardcodes a path or a
threshold. Every numbered script does:

    from config import PATHS, PARAMS, DESIGN, CONFIG

- `PATHS`  — frozen directory contract; accessors mkdir-on-access so a stage
             dir exists the first time it is written to.
- `PARAMS` — analysis thresholds, read from `analysis_config.yaml::thresholds`.
- `DESIGN` — experimental design (tissue levels, populations, contrast).
- `CONFIG` — the full parsed yaml (thread to figure-style as `config=`).

All relative paths resolve against the COMPARTMENT ROOT (this file's dir), so
scripts work regardless of the shell's CWD as long as they `from config import`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import yaml

# --- compartment root = this file's directory --------------------------------
ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "02_analysis" / "config" / "analysis_config.yaml"

with open(CONFIG_PATH) as _fh:
    CONFIG: Dict[str, Any] = yaml.safe_load(_fh) or {}

_paths = CONFIG.get("paths", {}) or {}


def _abs(rel: str) -> Path:
    """Resolve a config-relative path against the compartment root."""
    p = Path(rel)
    return p if p.is_absolute() else (ROOT / p)


@dataclass(frozen=True)
class _Paths:
    """Directory contract. Dir accessors mkdir-on-access; file accessors do not."""

    root: Path = ROOT
    results: Path = _abs(_paths.get("results", "03_results/"))
    objects: Path = _abs(_paths.get("objects", "03_results/objects/"))
    master: Path = _abs(_paths.get("master", "03_results/master/"))
    interactive: Path = _abs(_paths.get("interactive", "03_results/interactive/"))
    scratch: Path = _abs(_paths.get("scratch", "03_results/_scratch/"))

    # --- inputs (read-only) ---
    raw: Path = _abs(_paths.get("raw_gse160097", "00_data/GSE160097_JIA-SF-Treg/raw/"))
    samples_csv: Path = _abs(_paths.get("samples_gse160097", "00_data/GSE160097_JIA-SF-Treg/samples.csv"))
    signature_contract: Path = _abs(_paths.get("signature_contract", "../mouse_anchor/03_results/human_projection/"))

    _figures_subdir: str = _paths.get("stage_figures_subdir", "figures")
    _tables_subdir: str = _paths.get("stage_tables_subdir", "tables")

    # -- stage dirs (mkdir-on-access) --
    def stage_dir(self, stage: str) -> Path:
        d = self.results / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def figures(self, stage: str) -> Path:
        d = self.stage_dir(stage) / self._figures_subdir
        d.mkdir(parents=True, exist_ok=True)
        return d

    def tables(self, stage: str) -> Path:
        d = self.stage_dir(stage) / self._tables_subdir
        d.mkdir(parents=True, exist_ok=True)
        return d

    # -- root resources (mkdir-on-access) --
    def objects_dir(self) -> Path:
        self.objects.mkdir(parents=True, exist_ok=True)
        return self.objects

    def object(self, name: str) -> Path:
        """Checkpoint path `objects/<name>.h5ad` (dir ensured; file not created)."""
        self.objects_dir()
        stem = name if name.endswith((".h5ad", ".rds", ".parquet", ".csv")) else f"{name}.h5ad"
        return self.objects / stem

    def master_dir(self) -> Path:
        self.master.mkdir(parents=True, exist_ok=True)
        return self.master

    def master_file(self, name: str) -> Path:
        self.master_dir()
        return self.master / name

    def interactive_dir(self) -> Path:
        self.interactive.mkdir(parents=True, exist_ok=True)
        return self.interactive

    def scratch_dir(self) -> Path:
        self.scratch.mkdir(parents=True, exist_ok=True)
        return self.scratch


PATHS = _Paths()


@dataclass(frozen=True)
class _Params:
    """Analysis thresholds, sourced from `analysis_config.yaml::thresholds`."""

    _t: Dict[str, Any]

    def __getattr__(self, key: str) -> Any:  # dataclass __getattr__ fallback
        t = object.__getattribute__(self, "_t")
        if key in t:
            return t[key]
        raise AttributeError(f"threshold '{key}' not in analysis_config.yaml::thresholds")

    def get(self, key: str, default: Any = None) -> Any:
        return object.__getattribute__(self, "_t").get(key, default)


PARAMS = _Params(_t=CONFIG.get("thresholds", {}) or {})

# --- design (frozen, read straight from config) ------------------------------
DESIGN: Dict[str, Any] = CONFIG.get("design", {}) or {}
PROJECT: Dict[str, Any] = CONFIG.get("project", {}) or {}

# Convenience constants used across the pipeline.
SPECIES: str = PROJECT.get("species", "Homo sapiens")
SPECIES_DB: str = PROJECT.get("species_db", "HS")
POPULATIONS = list(DESIGN.get("populations", ["CD4_Treg", "CD4_Tcon", "CD8"]))
TISSUE_KEY = DESIGN.get("tissue_key", "tissue")
DONOR_KEY = DESIGN.get("donor_key", "donor")
POPULATION_KEY = DESIGN.get("population_key", "population")
TISSUE_NUM = DESIGN.get("tissue_levels", {}).get("synovial_fluid", "synovial_fluid")
TISSUE_DEN = DESIGN.get("tissue_levels", {}).get("peripheral_blood", "peripheral_blood")

# Short label map used in figures / frozen coarse labels.
COARSE_LABEL = {"CD4_Treg": "Treg", "CD4_Tcon": "Tcon", "CD8": "CD8"}


# --- one population palette for every figure ---------------------------------
class _PopulationColors(dict):
    """`{population: hex}` for the three sorted populations.

    Iteration, `len()` and `.keys()` see ONLY the three canonical figure labels
    (`Treg`, `Tcon`, `CD8`), so a script that builds a legend by walking the palette gets
    three entries rather than one per spelling. The long FACS gate names (`CD4_Treg`,
    `CD4_Tcon`) resolve on lookup instead, so a stage that carries them in an `obs` column
    or a table needs no translation at the call site.
    """

    def __init__(self, canonical: Dict[str, str], aliases: Dict[str, str]) -> None:
        super().__init__(canonical)
        self._aliases = dict(aliases)

    def __missing__(self, key: str) -> str:
        target = self._aliases.get(key)
        if target is None:
            raise KeyError(
                f"'{key}' is not a sorted population. Known: {sorted(self)} "
                f"(also accepted: {sorted(self._aliases)}).")
        return self[target]

    def __contains__(self, key: object) -> bool:
        return dict.__contains__(self, key) or key in self._aliases

    def get(self, key, default=None):  # dict.get bypasses __missing__
        try:
            return self[key]
        except KeyError:
            return default


def _resolve_named_hues(block: str) -> Dict[str, str]:
    """`{level: hex}` for one `colors.<block>` mapping of level -> `okabe_ito` key.

    Shared by the two categorical axes so both resolve identically and neither can be
    given a hex literal in the config.
    """
    colors = CONFIG.get("colors", {}) or {}
    okabe = colors.get("okabe_ito", {}) or {}
    named = colors.get(block, {}) or {}
    if not named:
        raise KeyError(
            f"analysis_config.yaml::colors.{block} is absent — every figure that colours "
            f"by {block} reads it, so a missing block is a config error rather than "
            "something to fall back from.")
    out: Dict[str, str] = {}
    for level, hue in named.items():
        if hue not in okabe:
            raise KeyError(
                f"colors.{block}.{level} names '{hue}', which is not a key of "
                f"colors.okabe_ito ({sorted(okabe)}). A level names a palette entry, "
                "never a hex literal.")
        out[level] = okabe[hue]
    return out


def population_colors() -> _PopulationColors:
    """The one population palette, resolved from `colors.populations`."""
    canonical = _resolve_named_hues("populations")
    aliases = {long_name: short_name for long_name, short_name in COARSE_LABEL.items()
               if short_name in canonical and long_name != short_name}
    return _PopulationColors(canonical, aliases)


def tissue_colors() -> Dict[str, str]:
    """The one tissue palette, resolved from `colors.tissue`.

    Keys are checked against `design.tissue_levels`, so renaming a tissue level cannot
    silently leave its colour behind under the old name.
    """
    resolved = _resolve_named_hues("tissue")
    levels = set((DESIGN.get("tissue_levels", {}) or {}).values())
    if levels and set(resolved) != levels:
        raise KeyError(
            f"colors.tissue keys {sorted(resolved)} do not match design.tissue_levels "
            f"{sorted(levels)}. The two axes of this compartment are declared together; "
            "renaming a tissue level means renaming its colour.")
    return resolved


# The two categorical axes. They must stay disjoint: several figures draw both, and one
# puts a tissue panel and a cell-state panel in the same row, where a shared hue would
# make the two axes indistinguishable. Asserted at import so the constraint cannot be
# broken by a config edit alone.
POPULATION_COLORS = population_colors()
TISSUE_COLORS: Dict[str, str] = tissue_colors()

_shared_hues = set(POPULATION_COLORS.values()) & set(TISSUE_COLORS.values())
if _shared_hues:
    raise ValueError(
        f"colors.populations and colors.tissue share {sorted(_shared_hues)}. The two "
        "categorical axes of this compartment must stay disjoint, because figures draw "
        "them side by side. Give one of the two axes a different okabe_ito entry.")
