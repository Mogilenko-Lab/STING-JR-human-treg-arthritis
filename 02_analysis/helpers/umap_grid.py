"""
02_analysis/helpers/umap_grid.py — one geometry for every UMAP strip in this compartment.
=========================================================================================
Both embedding viz scripts draw their panel rows through here, so strips from either stage
stack with panels at identical size.

Lengths are inches, divided by the canvas at draw time: a panel holds one fixed aspect and
the canvas width follows from the column count. Callers own which columns to draw, their
titles, their colour limits and the prose.

A score row draws through `scatter_cont(rescale=True, colorbar=False)` and hangs one
`row_colorbar` at the end, so the strip carries a single gradient. The grid reserves the same
colourbar slack either way, so a one-bar strip still stacks against a six-bar one.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# --- panel geometry, inches ---------------------------------------------------
# PANEL_ASPECT is width over height, held fixed for every panel of every strip, and each
# map's bounding box is padded out to match it. Panels therefore stay one size across both
# stages while the data keeps equal aspect, so a UMAP reads undistorted wherever it appears.
PANEL_ASPECT = 0.90
PANEL_H_IN = 3.60
PANEL_W_IN = PANEL_H_IN * PANEL_ASPECT
CB_PAD_IN = 0.10         # panel to colourbar
CB_W_IN = 0.14           # colourbar width
GUTTER_IN = 0.42         # colourbar to next panel
MARGIN_IN = 0.30         # canvas left/right margin

TITLE_BAND_IN = 0.95     # wrapped panel title, with clearance for a group label above it
LEGEND_BAND_IN = 0.80    # categorical legends, clear of the next row's labels
TOP_BAND_IN = 1.00       # figure title and subtitle
FOOTER_BAND_IN = 0.86    # standing lines, up to three

CELL_IN = PANEL_W_IN + CB_PAD_IN + CB_W_IN + GUTTER_IN
ROW_UNIT_IN = TITLE_BAND_IN + PANEL_H_IN + LEGEND_BAND_IN

# --- drawing parameters shared by both stages --------------------------------
CMAP = "viridis"         # one sequential scale for every continuous panel
CLIP = (2, 98)           # robust percentile clip
BOX_PAD = 0.03           # bounding-box margin, as a fraction of the box
SAMPLE_SEED = 0

UNIT_AUCELL = "AUCell"
UNIT_EXPR = "log-norm expression"
UNIT_MODULE = "module score"

# A rescaled bar is [0, 1] across the panel's clip, not the channel's units, and says so: it
# cannot be read for level. The clip it stands for goes in the caption's config line.
UNIT_RESCALED = "rescaled to panel clip"


def canvas(n_col: int, n_row: int = 1) -> tuple[float, float]:
    """Canvas inches for an `n_row` x `n_col` panel grid."""
    return (2 * MARGIN_IN + n_col * CELL_IN,
            TOP_BAND_IN + n_row * ROW_UNIT_IN + FOOTER_BAND_IN)


def fs(cfg: dict, key: str) -> float:
    """One font size, read from the config `figures:` block."""
    return float((cfg.get("figures", {}) or {})[key])


TITLE_CHARS = 22         # identifier characters a panel column carries on one line


def _wrap_ident(name: str, max_chars: int) -> list:
    """Break a set identifier on underscores to fit a column."""
    lines, cur = [], ""
    for part in name.split("_"):
        cand = part if not cur else f"{cur}_{part}"
        if len(cand) <= max_chars or not cur:
            cur = cand
        else:
            lines.append(cur)
            cur = part
    lines.append(cur)
    return lines


def set_titles(names: list, sizes: dict, max_chars: int = TITLE_CHARS) -> list:
    """Titles of uniform height: identifier wrapped to the column, then the size scored.

    Short identifiers take leading blank lines, so every panel title in a strip occupies the
    same number of lines and the row keeps one top edge.
    """
    wrapped = [_wrap_ident(n, max_chars) for n in names]
    depth = max(len(w) for w in wrapped)
    return ["\n" * (depth - len(w)) + "\n".join(w) + f"\n{sizes[n]} genes scored"
            for n, w in zip(names, wrapped)]


# =============================================================================
# Frame and axes
# =============================================================================
def sample_frame(df: pd.DataFrame, n: int | None, seed: int = SAMPLE_SEED) -> pd.DataFrame:
    """The one frame a figure's panels share, row order shuffled so groups paint evenly."""
    d = df if (n is None or n >= len(df)) else df.sample(n=n, random_state=seed)
    return d.sample(frac=1.0, random_state=seed).reset_index(drop=True)


def data_box(d: pd.DataFrame, pad: float = BOX_PAD):
    """One bounding box over the frame, padded out to PANEL_ASPECT.

    The shorter axis takes the padding, so the cloud fills its panel, the data keeps equal
    aspect, and every panel of every strip shares one box.
    """
    x = d["x"].to_numpy(dtype=float)
    y = d["y"].to_numpy(dtype=float)
    cx, cy = (x.min() + x.max()) / 2.0, (y.min() + y.max()) / 2.0
    xr, yr = x.max() - x.min(), y.max() - y.min()
    if xr / yr < PANEL_ASPECT:
        xr = yr * PANEL_ASPECT
    else:
        yr = xr / PANEL_ASPECT
    hx, hy = xr / 2.0 * (1.0 + pad), yr / 2.0 * (1.0 + pad)
    return (cx - hx, cx + hx), (cy - hy, cy + hy)


def grid(fig, n_col: int, n_row: int = 1) -> list:
    """Panel boxes at PANEL_ASPECT, row-major from the top. One axes list per row."""
    w, h = fig.get_size_inches()
    side_x, side_y = PANEL_W_IN / w, PANEL_H_IN / h
    rows = []
    for r in range(n_row):
        top_in = TOP_BAND_IN + r * ROW_UNIT_IN + TITLE_BAND_IN
        bottom = 1.0 - (top_in + PANEL_H_IN) / h
        rows.append([fig.add_axes([(MARGIN_IN + i * CELL_IN) / w, bottom, side_x, side_y])
                     for i in range(n_col)])
    return rows


def _frame_panel(ax, xlim, ylim, title: str, cfg: dict):
    """Tick-free UMAP panel on the shared bounding box."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=fs(cfg, "base_size"), pad=8, linespacing=1.25)


def _colorbar(fig, ax, sc, unit: str, cfg: dict):
    """Colourbar in the slack beside its panel, labelled by its unit."""
    w, _ = fig.get_size_inches()
    pos = ax.get_position()
    cax = fig.add_axes([pos.x1 + CB_PAD_IN / w, pos.y0 + 0.08 * pos.height,
                        CB_W_IN / w, 0.84 * pos.height])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label(unit, fontsize=fs(cfg, "axis_title_size"))
    cb.ax.tick_params(labelsize=fs(cfg, "axis_text_size"))
    cb.outline.set_visible(False)


def clip_of(frame: pd.DataFrame, cols: list) -> tuple[float, float]:
    """One 2nd-to-98th percentile clip over pooled columns, putting sets on one scale."""
    pooled = np.concatenate([frame[c].to_numpy(dtype=float) for c in cols])
    lo, hi = np.nanpercentile(pooled, list(CLIP))
    return float(lo), float(hi)


def point_size(n_cells: int) -> float:
    """Point size matched to the drawn density, keeping panels evenly filled."""
    return 1.6 if n_cells > 40_000 else 2.4


# =============================================================================
# Panels
# =============================================================================
def scatter_cont(ax, d: pd.DataFrame, col: str, title: str, xlim, ylim, cfg: dict,
                 vlim=None, unit: str = UNIT_AUCELL, rescale: bool = False,
                 colorbar: bool = True):
    """Continuous panel: robust clip, high values drawn last.

    `vlim` supplies limits from elsewhere, which is what makes a subset panel readable
    against its full-object twin.

    `rescale` maps the panel onto [0, 1] across those same limits, so panels of different
    range share one bar. The clip is unchanged, so the picture is too; what is lost is level,
    since equal brightness now means equal position within unequal clips.

    `colorbar=False` suppresses the per-panel bar and returns the mappable for `row_colorbar`.
    """
    _frame_panel(ax, xlim, ylim, title, cfg)
    v = d[col].to_numpy(dtype=float)
    lo, hi = vlim if vlim is not None else np.nanpercentile(v, list(CLIP))
    if rescale:
        if not hi > lo:
            raise ValueError(
                f"{col}: colour limits {lo} to {hi} leave no span, so the panel cannot be "
                "rescaled onto a shared bar.")
        v = (v - lo) / (hi - lo)
        lo, hi = 0.0, 1.0
    order = np.argsort(v)
    sc = ax.scatter(d["x"].to_numpy()[order], d["y"].to_numpy()[order], c=v[order],
                    s=point_size(len(d)), cmap=CMAP, vmin=lo, vmax=hi,
                    linewidths=0, rasterized=True)
    if colorbar:
        _colorbar(ax.figure, ax, sc, unit, cfg)
    return sc


def row_colorbar(fig, ax, sc, unit: str, cfg: dict):
    """One bar for a row, in the slot a per-panel bar would take, so panel positions are
    untouched and a one-bar strip still stacks against a six-bar one."""
    _colorbar(fig, ax, sc, unit, cfg)


def scatter_cat(ax, d: pd.DataFrame, col: str, palette: dict, title: str, xlim, ylim,
                cfg: dict, order: list, labels: dict | None = None, ncol: int = 2):
    """Categorical panel: one shuffled scatter call, legend in the band below."""
    _frame_panel(ax, xlim, ylim, title, cfg)
    values = d[col].astype(str).to_numpy()
    colours = np.array([palette[v] for v in values])
    ax.scatter(d["x"].to_numpy(), d["y"].to_numpy(), c=colours, s=point_size(len(d)),
               linewidths=0, rasterized=True)
    present = [c for c in order if c in set(values)]
    show = labels or {}
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[c],
                      markeredgecolor="none", markersize=8, label=show.get(c, c))
               for c in present]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.015),
              ncol=ncol, frameon=False, fontsize=fs(cfg, "legend_text_size"),
              handletextpad=0.25, columnspacing=1.0, borderpad=0.0)


# =============================================================================
# Dressing
# =============================================================================
def pin_canvas(fig):
    """Pin the saved extent to the declared canvas.

    `save_figure` writes with `bbox_inches="tight"`, which crops to ink, so two strips whose
    clouds reach different edges come out different sizes and their panel grids land at
    different offsets. An edgeless full-canvas rectangle draws nothing and holds the extent,
    so every strip saves at its declared size and stacks exactly.
    """
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure,
                                 fill=False, edgecolor="none", linewidth=0))


def dress(fig, suptitle: str, subtitle: str, footer_lines: list, cfg: dict,
          colour: str = "#000000"):
    """Title block above the grid, standing lines at the foot, canvas extent pinned."""
    pin_canvas(fig)
    _, h = fig.get_size_inches()
    fig.suptitle(suptitle, fontsize=fs(cfg, "title_size"), y=1.0 - 0.16 / h)
    fig.text(0.5, 1.0 - 0.46 / h, subtitle, ha="center", va="top",
             fontsize=fs(cfg, "subtitle_size"))
    y = (FOOTER_BAND_IN - 0.16) / h
    for line in footer_lines:
        fig.text(0.5, y, line, ha="center", va="center",
                 fontsize=fs(cfg, "axis_title_size"), color=colour)
        y -= 0.235 / h


def group_bands(fig, axes: list, bands: list, cfg: dict, colour: str = "#000000"):
    """Label spans of columns and rule them off from each other.

    `bands` is [(label, first_col, last_col, rule_after)]. The rules keep provenance visible
    in a row of six panels.
    """
    _, h = fig.get_size_inches()
    for label, i0, i1, rule_after in bands:
        p0, p1 = axes[i0].get_position(), axes[i1].get_position()
        fig.text((p0.x0 + p1.x1) / 2.0, p0.y1 + (TITLE_BAND_IN - 0.14) / h, label,
                 ha="center", va="bottom", fontsize=fs(cfg, "axis_title_size"),
                 color=colour, style="italic")
        if rule_after and i1 + 1 < len(axes):
            x = (p1.x1 + axes[i1 + 1].get_position().x0) / 2.0
            fig.add_artist(Line2D([x, x], [p0.y0, p0.y1 + (TITLE_BAND_IN - 0.20) / h],
                                  color=colour, linewidth=0.8, alpha=0.5))
