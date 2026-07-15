import marimo

# JIA synovial-Treg reactive review.
# Reads frozen tables under 03_results (compute is done upstream by the scoring,
# embedding, and harvest-readout scripts). This app only visualizes and narrates.
# The mouse WT_heat overlap in JIA synovial-fluid T cells reads as consistent with
# the mouse 39 C stress axis, correlative only. The hypoxia surface is a readout,
# never a claim that HIF or fever drives the human state. The per-cell scores are a
# secondary annotation tier, held apart from the donor-level pseudobulk NES.

app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    import yaml
    from plotly.subplots import make_subplots

    return Path, go, make_subplots, np, pd, px, yaml


@app.cell
def _(Path):
    # Resolve the project root by walking up until 03_results/interactive appears,
    # so the app loads its tables whether launched from the repo root, the app dir,
    # or under the umbrella container.
    def _find_root():
        bases = []
        try:
            bases.append(Path(__file__).resolve())
        except NameError:
            pass
        bases.append(Path.cwd().resolve())
        for b in bases:
            p = b
            for _ in range(8):
                if (p / "03_results" / "interactive").exists():
                    return p
                p = p.parent
        raise FileNotFoundError("could not locate 03_results/interactive")

    PROJECT_ROOT = _find_root()
    INTERACTIVE = PROJECT_ROOT / "03_results" / "interactive"
    RESULTS = PROJECT_ROOT / "03_results"
    CONFIG = PROJECT_ROOT / "02_analysis" / "config" / "analysis_config.yaml"
    return CONFIG, INTERACTIVE, RESULTS


@app.cell
def _(CONFIG, INTERACTIVE, RESULTS, pd, yaml):
    # The per-cell embedding + readout substrate: one row per cell, x / y UMAP,
    # frozen cell-state label, tissue, donor, mitochondrial fraction, and the AUCell
    # hypoxia / UPR surfaces alongside the WT_heat, eTreg, and HSP per-cell scores.
    cells = pd.read_parquet(INTERACTIVE / "08_harvest_readout.parquet")

    # Primary tier: donor-pseudobulk fgsea NES for the mouse WT_heat set on the
    # SF-vs-PB signed ranked list, one table per sorted population.
    _gsea_tag = {"Treg": "treg", "Tcon": "tcon", "CD8": "cd8"}
    _gsea_parts = []
    for _pop, _tag in _gsea_tag.items():
        _f = RESULTS / "05_scoring" / "tables" / f"gsea_pseudobulk_{_tag}.csv"
        _df = pd.read_csv(_f)
        _df["cell_state"] = _pop
        _gsea_parts.append(_df)
    gsea = pd.concat(_gsea_parts, ignore_index=True)

    # Secondary tier: per-cell AUCell SF-vs-PB donor-level standardized mean
    # difference, carried on its own scale and never pooled with the NES above.
    eff = pd.read_csv(RESULTS / "master" / "effect_sizes_treg_arthritis.csv")
    smd = eff[eff["effect_metric"] == "percell_auc_smd"].copy()

    # The SF-vs-PB hypoxia / UPR readout summary by cell state.
    readout = pd.read_csv(
        RESULTS / "08_harvest_readout" / "tables" / "harvest_readout_summary.csv"
    )

    heat_hypoxia_dir = RESULTS / "09_heat_hypoxia" / "tables"
    heat_hypoxia_purge = pd.read_csv(
        heat_hypoxia_dir / "gene_purge_nes_comparison.csv"
    )
    heat_hypoxia_coloc = pd.read_csv(
        heat_hypoxia_dir / "heat_hypoxia_colocalization.csv"
    )
    heat_hypoxia_leadingedge = pd.read_csv(
        heat_hypoxia_dir / "leadingedge_composition.csv"
    )

    # The current OR-gated drafted subset context, open to revision.
    or_union = pd.read_csv(
        RESULTS / "07_embedding" / "tables" / "or_union_membership.csv"
    )
    hook_lineage = pd.read_csv(
        RESULTS / "07_embedding" / "tables" / "hook_per_lineage_summary.csv"
    )
    return (
        cells,
        eff,
        gsea,
        heat_hypoxia_coloc,
        heat_hypoxia_leadingedge,
        heat_hypoxia_purge,
        hook_lineage,
        or_union,
        readout,
        smd,
    )


@app.cell
def _(CONFIG, yaml):
    with open(CONFIG) as _fh:
        _cfg = yaml.safe_load(_fh)
    go_decision = _cfg.get("decisions", {}).get("go_no_go", {})
    return (go_decision,)


@app.cell
def _(mo):
    mo.md(
        r"""
        # The mouse 39 °C Treg stress axis in JIA synovial-fluid T cells

        I built this one reactive review to re-read the JIA synovial-Treg result
        end to end: does the mouse 39 °C fever program (`WT_heat`, derived from the
        internal mouse iTreg experiment) show interpretable enrichment in sorted
        synovial-fluid (SF) T cells relative to paired peripheral blood (PB) in
        GSE160097, and does that enrichment concentrate in Tregs or read across the
        whole T compartment?

        **The headline read.** The mouse `WT_heat` up-program is **enriched** in JIA
        synovial-fluid T cells versus paired blood, and the enrichment is **broad
        across the sorted populations — pan-T, not Treg-preferential**. The running
        sum packs the core set at the top of every ranked list: NES ≈ 2.5 in Treg,
        2.6 in Tcon, and 2.1 in CD8, all clearing FDR far below 1e-6. This is a
        positive, reproducible signal read as **consistent with** the mouse stress
        axis, a correlative observation and not evidence that fever or HIF drives the
        human state.

        **What this review adds.** Beyond the primary enrichment I surface two
        per-cell readouts harvested from the same frozen labels — the MSigDB Hallmark
        **hypoxia** and **unfolded-protein-response (UPR)** AUCell surfaces — as
        exploratory annotation. The hypoxia surface carries a synovial-Treg-leading
        high pocket. The UPR surface reads flat. Both are readouts of the tissue
        state, correlative, and carry no HIF-causality claim.

        **How to read the evidence tiers.** The donor-pseudobulk NES table is the
        primary tier. The per-cell AUCell scores in the embedding explorer and the
        hypoxia / UPR readouts are a secondary, corroborative annotation tier. I keep
        the two tiers on separate scales and never pool them into one ranking. The
        UMAP below is a map of the cells for visualization, read as annotation and
        never as the evidence for a claim.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Primary evidence — `WT_heat` NES on donor-pseudobulk SF-vs-PB ranked lists

        Pre-ranked fgsea of the mouse `WT_heat` up / down sets against the
        SF-vs-PB donor-pseudobulk signed-Wald ranked list, one column per sorted
        population. Positive NES (red) means the set packs toward the SF-high end.
        The `WT_heat_up` row is the number that settles the claim: it is strongly
        positive and significant in Treg, Tcon, and CD8 alike, so the signal is
        pan-T rather than Treg-restricted. Stars mark BH-FDR (`***` < 0.001,
        `**` < 0.01, `*` < 0.05).
        """
    )
    return


@app.cell
def _(go, gsea, mo):
    def _stars(p):
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    _rows = ["WT_heat_up", "WT_heat_down"]
    _cols = ["Treg", "Tcon", "CD8"]
    _z, _t = [], []
    for _r in _rows:
        _zr, _tr = [], []
        for _c in _cols:
            _m = gsea[(gsea["pathway_id"] == _r) & (gsea["cell_state"] == _c)]
            if len(_m):
                _nes = float(_m["nes"].iloc[0])
                _padj = float(_m["padj"].iloc[0])
                _zr.append(_nes)
                _tr.append(f"{_nes:+.2f}<br>{_stars(_padj)}")
            else:
                _zr.append(None)
                _tr.append("")
        _z.append(_zr)
        _t.append(_tr)

    _fig = go.Figure(
        go.Heatmap(
            z=_z,
            x=_cols,
            y=_rows,
            text=_t,
            colorscale="RdBu",
            reversescale=True,
            zmid=0,
            zmin=-2.7,
            zmax=2.7,
            texttemplate="%{text}",
            textfont=dict(size=16),
            xgap=3,
            ygap=3,
            colorbar=dict(title="NES"),
        )
    )
    _fig.update_yaxes(autorange="reversed", tickfont=dict(size=15))
    _fig.update_xaxes(tickfont=dict(size=15))
    _fig.update_layout(
        height=340,
        template="plotly_white",
        title=dict(text="WT_heat NES — SF vs PB, by sorted population", x=0.5),
        margin=dict(l=10, r=10, t=60, b=10),
    )

    # A same-numbers markdown table sits beside the heatmap so the reader can lift the
    # exact NES / FDR / set size for any cell.
    _tbl_rows = []
    for _, _g in gsea[gsea["pathway_id"].isin(_rows)].iterrows():
        _tbl_rows.append(
            f"| {_g['cell_state']} | {_g['pathway_id']} | {float(_g['nes']):+.3f} "
            f"| {float(_g['pvalue']):.2g} | {float(_g['padj']):.2g} | "
            f"{int(_g['set_size'])} |"
        )
    _tbl = (
        "| population | signature | NES | p | FDR | set size |\n"
        "|---|---|---|---|---|---|\n" + "\n".join(_tbl_rows)
    )
    _cap = mo.md(
        "**Reading the panel.** `WT_heat_up` runs strongly SF-high in every sorted "
        "population — Treg NES +2.51, Tcon +2.57, CD8 +2.05, each at FDR far below "
        "1e-6 — so the enrichment is pan-T, not Treg-preferential. `WT_heat_down` "
        "sits near zero and non-significant, consistent with the up-limb carrying the "
        "signal. This donor-pseudobulk NES is the primary evidence the review rests "
        "on.\n\n" + _tbl
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Heat-vs-hypoxia check inside JIA SF

        Synovial fluid is a hypoxic tissue niche, so the SF-vs-PB `WT_heat` enrichment
        could reflect the hypoxic microenvironment rather than a heat program. I ask
        that here three ways, all within JIA. Does the enrichment survive removing the
        genes `WT_heat_up` shares with `HALLMARK_HYPOXIA` (primary donor-pseudobulk
        tier)? Do per-cell heat and hypoxia scores co-localize in SF cells, and which
        biological programs do the enriching leading-edge genes actually represent
        (both secondary annotation)? These do not fully deconfound tissue context —
        heat and hypoxia separate cleanly only where an elevated temperature occurs
        without a hypoxic niche — but they do separate a hypoxia-overlap signal from
        the rest.
        """
    )
    return


@app.cell
def _(go, heat_hypoxia_coloc, heat_hypoxia_leadingedge, heat_hypoxia_purge, mo):
    _order = ["Treg", "Tcon", "CD8"]
    _p = heat_hypoxia_purge.set_index("population").loc[_order].reset_index()
    _c = heat_hypoxia_coloc[
        (heat_hypoxia_coloc["level"] == "cell")
        & (heat_hypoxia_coloc["method"] == "spearman")
    ].set_index("population").loc[_order].reset_index()
    _l = heat_hypoxia_leadingedge.set_index("population").loc[_order].reset_index()

    _fig = go.Figure()
    _fig.add_trace(
        go.Bar(
            x=_p["population"],
            y=_p["NES_full"],
            name="full WT_heat_up",
            marker_color="#B35806",
        )
    )
    _fig.add_trace(
        go.Bar(
            x=_p["population"],
            y=_p["NES_purged"],
            name="hypoxia-overlap purged",
            marker_color="#0072B2",
        )
    )
    _fig.update_layout(
        barmode="group",
        height=360,
        template="plotly_white",
        title=dict(text="WT_heat_up NES before and after hypoxia-gene purge", x=0.5),
        yaxis_title="NES",
        xaxis_title="",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0.5, xanchor="center"),
        margin=dict(l=20, r=20, t=80, b=30),
    )

    _rows = []
    for _, _r in _p.iterrows():
        _pop = _r["population"]
        _cc = _c[_c["population"] == _pop].iloc[0]
        _ll = _l[_l["population"] == _pop].iloc[0]
        _act = float(_ll["frac_effector_activation"]) + float(_ll["frac_immediate_early_stress"])
        _rows.append(
            f"| {_pop} | {float(_r['NES_full']):+.2f} | "
            f"{float(_r['NES_purged']):+.2f} | {float(_r['padj_purged']):.2g} | "
            f"{float(_cc['r']):+.2f} | {int(_cc['n']):,} | "
            f"{100 * _act:.0f}% | "
            f"{100 * float(_ll['frac_hypoxia_HIF']):.0f}% | "
            f"{100 * float(_ll['frac_heat_shock_proteostasis']):.0f}% |"
        )
    _tbl = (
        "| population | full NES | purged NES | purged FDR | SF cell Spearman r | SF cells | activation + immediate-early | hypoxia | heat-shock |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
        + "\n".join(_rows)
    )
    _cap = mo.md(
        "The `WT_heat_up` enrichment is hypoxia-independent but not thermal-specific. "
        "Removing the 18 genes it shares with `HALLMARK_HYPOXIA` barely dents the NES — "
        "Treg +2.38, Tcon +2.39, CD8 +1.90, all at FDR < 5e-5 — so the SF-vs-PB signal "
        "is not carried by the hypoxia-overlap genes. Within SF the per-cell heat and "
        "hypoxia scores correlate only weakly (Spearman 0.08 to 0.20), so heat-high and "
        "hypoxia-high are largely different cells. But the leading edge that drives the "
        "enrichment is mostly generic T-cell activation and immediate-early genes "
        "(about two-thirds together), with a hypoxia minority and only a trace of "
        "classic heat-shock or proteostasis genes (HSPA1A, HSPH1, CLU). I read this "
        "honestly: the mouse `WT_heat` overlap with inflamed-joint T cells is real and "
        "is not a hypoxia artefact, but at gene resolution it reads as an activation "
        "and immediate-early program more than a thermal-specific one, so heat is a "
        "loose label here for what is largely in-vivo T-cell activation.\n\n" + _tbl
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Secondary (corroborative) — per-cell `WT_heat_up` AUCell SF-vs-PB SMD

        The per-cell lens is AUCell — area under each cell's gene-recovery curve for
        the up-set, rank-based and bounded in [0, 1], robust to library size and
        composition. I summarize it as a donor-level standardized mean difference
        (SMD) of SF versus PB per population, with a 95% confidence interval. It
        stays on its own footing and is never pooled with the pseudobulk NES above.
        """
    )
    return


@app.cell
def _(mo, smd):
    _order = {"Treg": 0, "Tcon": 1, "CD8": 2}
    _s = smd.copy()
    _s["_k"] = _s["cell_state"].map(_order).fillna(9)
    _s = _s.sort_values("_k")
    _rows = []
    for _, _r in _s.iterrows():
        _rows.append(
            f"| {_r['cell_state']} | {_r['signature']} | "
            f"{float(_r['estimate']):+.2f} | "
            f"[{float(_r['ci_low']):.2f}, {float(_r['ci_high']):.2f}] | "
            f"{float(_r['pvalue']):.2g} | {int(_r['n_donors'])} | "
            f"{int(_r['n_cells']):,} |"
        )
    _tbl = (
        "| population | signature | SMD | 95% CI | p | donors | cells |\n"
        "|---|---|---|---|---|---|---|\n" + "\n".join(_rows)
    )
    mo.md(
        "Per-cell `WT_heat_up` AUCell SF-vs-PB donor-level SMD, the corroborative "
        "tier. Every population shifts SF-high in the same direction as the primary "
        "NES, and the effect is largest in Treg (SMD +5.53) yet clearly present in "
        "Tcon and CD8, echoing the pan-T reading. This annotates the primary result "
        "and does not stand as independent evidence.\n\n" + _tbl
    )
    return


@app.cell
def _(mo):
    mo.md(r"""## Reactive embedding explorer""")
    return


@app.cell
def _(cells, mo):
    # I expose one curated, AUCell-only signature surface on the right panel rather
    # than every scored column. The hypoxia / UPR UCell twins are dropped so the two
    # engines never sit side by side, and the surface groups the WT_heat, eTreg, and
    # HSP per-cell scores with the two Hallmark AUCell readouts and the mitochondrial
    # fraction as an annotation channel.
    _curated = [
        ("mouse anchor — WT_heat (updown)", "WT_heat_updown"),
        ("effector — eTreg", "score_eTreg"),
        ("stress — HSP", "score_HSP"),
        ("readout — HALLMARK_HYPOXIA (AUCell)", "HALLMARK_HYPOXIA_AUCell"),
        (
            "readout — HALLMARK_UPR (AUCell)",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
        ),
        ("annotation — pct_counts_mt", "pct_counts_mt"),
    ]
    color_options = {}
    for _lab, _col in _curated:
        if _col in cells.columns:
            color_options[_lab] = _col

    # The left reference panel is always the frozen cell-state label.
    _meta_cat = {"cell state": "coarse_label", "tissue": "tissue", "donor": "donor"}
    left_color_by = mo.ui.dropdown(
        options=_meta_cat,
        value="cell state",
        label="Left panel — reference (categorical)",
    )
    # The right exploration panel defaults to the hypoxia readout so the new surface
    # is visible out of the box, cell state on the left against hypoxia on the right.
    color_by = mo.ui.dropdown(
        options=color_options,
        value="readout — HALLMARK_HYPOXIA (AUCell)",
        label="Right panel — color by",
    )

    _states = sorted(cells["coarse_label"].unique().tolist())
    ct_filter = mo.ui.multiselect(
        options=_states, value=_states, label="Cell states"
    )
    _tissues = sorted(cells["tissue"].unique().tolist())
    tissue_filter = mo.ui.multiselect(
        options=_tissues, value=_tissues, label="Tissue"
    )

    _n = len(cells)
    cap_slider = mo.ui.slider(
        start=5000,
        stop=_n,
        step=5000,
        value=min(40000, _n),
        label="Max points drawn (WebGL)",
        show_value=True,
    )

    controls = mo.vstack(
        [
            mo.hstack([left_color_by, color_by], justify="start", gap=1.5),
            mo.hstack([ct_filter, tissue_filter], justify="start", gap=1.5),
            cap_slider,
        ]
    )
    controls
    return cap_slider, color_by, ct_filter, left_color_by, tissue_filter


@app.cell
def _(
    cap_slider,
    cells,
    color_by,
    ct_filter,
    left_color_by,
    mo,
    px,
    tissue_filter,
):
    # ONE filtered, sampled cell set feeds BOTH panels. The cell-state and tissue
    # multiselects and the point cap all act on this single frame, so excluding a
    # population drops it from both maps at once and the two panels always depict the
    # identical cells at identical UMAP coordinates. I sample once here — never twice
    # — so the point sets can't drift apart.
    _d = cells[
        cells["coarse_label"].isin(ct_filter.value)
        & cells["tissue"].isin(tissue_filter.value)
    ]
    _total = len(_d)
    _drawn = _d
    if _total > cap_slider.value:
        _drawn = _d.sample(n=cap_slider.value, random_state=0)
    _drawn = _drawn.reset_index(drop=True)

    # Shared square framing. I read one bounding box off the FULL object's x / y so
    # both panels carry the identical ranges and stay directly comparable, holding
    # 1 unit on x equal to 1 unit on y so each map renders as a true square.
    _x = cells["x"]
    _y = cells["y"]
    _padx = 0.03 * (float(_x.max()) - float(_x.min()))
    _pady = 0.03 * (float(_y.max()) - float(_y.min()))
    _shared_x = [float(_x.min()) - _padx, float(_x.max()) + _padx]
    _shared_y = [float(_y.min()) - _pady, float(_y.max()) + _pady]

    def _cat_scatter(col):
        _plot = _drawn.assign(_color=_drawn[col].astype(str))
        _cat_order = sorted(cells[col].astype(str).unique().tolist())
        return px.scatter(
            _plot,
            x="x",
            y="y",
            color="_color",
            render_mode="webgl",
            labels={"_color": col},
            category_orders={"_color": _cat_order},
        )

    def _cont_scatter(col):
        return px.scatter(
            _drawn,
            x="x",
            y="y",
            color=col,
            render_mode="webgl",
            color_continuous_scale="viridis",
        )

    # Left = the cell-state reference map, always categorical.
    _left_col = left_color_by.value
    _fig_l = _cat_scatter(_left_col)
    _fig_l.update_layout(title=dict(text=f"reference ({_left_col})", x=0.5))

    # Right = the exploration map, driven by the color-by control. A continuous
    # signature by default, still categorical if a metadata field is chosen.
    _right_col = color_by.value
    _right_is_cat = _right_col in ("coarse_label", "tissue", "donor")
    _fig_r = _cat_scatter(_right_col) if _right_is_cat else _cont_scatter(_right_col)
    _fig_r.update_layout(title=dict(text=f"exploration ({_right_col})", x=0.5))

    for _fig in (_fig_l, _fig_r):
        _fig.update_traces(marker=dict(size=3, opacity=0.55))
        _fig.update_layout(
            width=600,
            height=660,
            template="plotly_white",
            # Legend and colorbar BOTH at the bottom, so neither eats horizontal
            # room. The left and right panels then keep an identical plotting width
            # and, with the square aspect below, render as two equal squares.
            legend=dict(
                title="",
                itemsizing="constant",
                orientation="h",
                yanchor="top",
                y=-0.16,
                x=0.5,
                xanchor="center",
            ),
            coloraxis_colorbar=dict(
                orientation="h",
                yanchor="top",
                y=-0.16,
                x=0.5,
                xanchor="center",
                thickness=12,
                len=0.85,
                title="",
            ),
            margin=dict(l=10, r=10, t=48, b=110),
            dragmode="lasso",
            xaxis=dict(title="UMAP 1", range=_shared_x, constrain="domain"),
            yaxis=dict(
                title="UMAP 2",
                range=_shared_y,
                scaleanchor="x",
                scaleratio=1,
                constrain="domain",
            ),
            # One shared uirevision, never recomputed from the selection, holds zoom /
            # pan and legend trace-visibility across figure updates. The same value on
            # both panels keeps their views in step, coherent because they share
            # coordinates.
            uirevision="jia-embed",
        )

    # Map each drawn cell to a stable position, then record which drawn cells sit in
    # each plotly trace. A lasso returns (curveNumber, pointIndex) pairs, and these
    # maps turn those pairs back into rows of the shown frame.
    _coord_pos = {}
    for _pos, (_a, _b) in enumerate(
        zip(_drawn["x"].tolist(), _drawn["y"].tolist())
    ):
        _coord_pos[(_a, _b)] = _pos

    def _curve_map(fig):
        _m = []
        for _tr in fig.data:
            _m.append(
                [
                    _coord_pos.get((float(_px), float(_py)))
                    for _px, _py in zip(list(_tr.x), list(_tr.y))
                ]
            )
        return _m

    curve_map_l = _curve_map(_fig_l)
    curve_map_r = _curve_map(_fig_r)
    drawn = _drawn

    embed_l = mo.ui.plotly(_fig_l)
    embed_r = mo.ui.plotly(_fig_r)

    _caption = mo.md(
        f"**{len(_drawn):,}** of **{_total:,}** filtered cells drawn "
        f"(cap {cap_slider.value:,}, WebGL `Scattergl`). Two square UMAPs of the same "
        f"cells on one shared scale, so one unit on UMAP 1 matches one unit on UMAP 2 "
        f"and both panels share the ranges read off the full object. The left is the "
        f"frozen cell-state reference (categorical palette), the right is the selected "
        f"readout (viridis continuous, AUCell-only surface). Both panels draw the "
        f"identical sampled cells at identical coordinates under one cell-state / "
        f"tissue filter and one point cap, so a region that lights up on the right "
        f"reads off the cell-state colors at the same position on the left. Zoom, pan, "
        f"and legend show/hide persist and stay coordinated across both panels as I "
        f"switch readouts. Lasso a region on either map with the lasso tool in the plot "
        f"toolbar and the panel directly below describes those cells. This UMAP is a "
        f"map for visualization, read as annotation and never as evidence. Whenever the "
        f"drawn count sits below the filtered total, the view is a random sample at the "
        f"visible cap above."
    )
    # The controls sit directly above these maps (previous cell), so the maps render
    # first here — nothing separates a control from the figure it drives. The
    # how-to-read narrative folds into a collapsed panel just below the maps, and the
    # lasso cluster-characteristics summary still renders in the cell directly below.
    mo.vstack(
        [
            mo.hstack([embed_l, embed_r], widths="equal", gap=1),
            mo.accordion({"How to read these two maps": _caption}),
        ]
    )
    return curve_map_l, curve_map_r, drawn, embed_l, embed_r


@app.cell
def _(curve_map_l, curve_map_r, drawn, embed_l, embed_r, go, mo):
    # Lasso either embedding above -> this cell describes the lassoed cells. It reads
    # each plotly selection back to rows of the shown frame through the trace maps,
    # then summarizes composition and per-cell readout means. Everything here is a
    # descriptive, exploratory read over the frozen per-cell columns, a secondary
    # annotation tier held apart from the donor-level pseudobulk enrichment.
    def _sel_positions(value, cmap):
        _found = set()
        for _pt in value or []:
            _c = _pt.get("curveNumber")
            _i = _pt.get("pointIndex")
            if _i is None:
                _i = _pt.get("pointNumber")
            if _c is None or _i is None:
                continue
            if 0 <= _c < len(cmap) and 0 <= _i < len(cmap[_c]):
                _p = cmap[_c][_i]
                if _p is not None:
                    _found.add(int(_p))
        return _found

    _pos = _sel_positions(embed_l.value, curve_map_l) | _sel_positions(
        embed_r.value, curve_map_r
    )

    # The per-cell readouts, shown as means over the selection. The mitochondrial
    # fraction stays on the color-by control only, kept off this bar because it is a
    # percentage rather than a bounded score and would flatten the other bars.
    _score_cols = {
        "WT_heat (updown)": "WT_heat_updown",
        "eTreg": "score_eTreg",
        "HSP": "score_HSP",
        "HALLMARK_HYPOXIA": "HALLMARK_HYPOXIA_AUCell",
        "HALLMARK_UPR": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_AUCell",
    }
    _score_cols = {k: v for k, v in _score_cols.items() if v in drawn.columns}

    if not _pos:
        _summary = mo.md(
            "**Cluster characteristics — exploratory.** Lasso a region on either map "
            "above to read that region's makeup here. I pick the lasso tool in the "
            "plot toolbar, draw around a cluster, and this panel reports how many "
            "cells fall inside, their cell-state, tissue, and donor composition, and "
            "each per-cell readout's mean in the selection against the other shown "
            "cells. The summary is a descriptive read over the frozen per-cell scores, "
            "a secondary annotation tier held apart from the donor-level pseudobulk "
            "enrichment, so it annotates the map rather than standing as evidence."
        )
    else:
        _sel = drawn.iloc[sorted(_pos)]
        _rest = drawn.drop(index=_sel.index)
        _n_sel = len(_sel)
        _n_shown = len(drawn)
        _pct = 100.0 * _n_sel / _n_shown if _n_shown else 0.0

        def _comp(col):
            _s = _sel[col].value_counts(normalize=True)
            _o = drawn[col].value_counts(normalize=True)
            _rows = [
                f"| {_k} | {100 * _s.get(_k, 0.0):.1f}% | {100 * _o.get(_k, 0.0):.1f}% |"
                for _k in _o.index
            ]
            return "\n".join(_rows)

        _ct_tbl = (
            "| cell state | selected | shown |\n|---|---|---|\n"
            + _comp("coarse_label")
        )
        _ts_tbl = (
            "| tissue | selected | shown |\n|---|---|---|\n" + _comp("tissue")
        )
        _dn_tbl = "| donor | selected | shown |\n|---|---|---|\n" + _comp("donor")

        _labels = list(_score_cols.keys())
        _sel_means = [float(_sel[c].mean()) for c in _score_cols.values()]
        _rest_means = [
            float(_rest[c].mean()) if len(_rest) else 0.0
            for c in _score_cols.values()
        ]
        _fig = go.Figure()
        _fig.add_trace(
            go.Bar(
                x=_labels,
                y=_sel_means,
                name="lassoed cells",
                marker_color="#c0392b",
            )
        )
        _fig.add_trace(
            go.Bar(
                x=_labels,
                y=_rest_means,
                name="other shown cells",
                marker_color="#7f8c8d",
            )
        )
        _fig.update_layout(
            barmode="group",
            height=360,
            template="plotly_white",
            yaxis_title="mean per-cell score",
            legend=dict(orientation="h", y=1.15, x=0),
            margin=dict(l=10, r=10, t=30, b=10),
            font=dict(size=13),
        )

        _summary = mo.vstack(
            [
                mo.md(
                    f"**Cluster characteristics — exploratory.** **{_n_sel:,}** cells "
                    f"lassoed, **{_pct:.1f}%** of the **{_n_shown:,}** shown. The bars "
                    f"give each per-cell readout's mean over the lassoed cells against "
                    f"the other shown cells, so I can contrast a hypoxia-high pocket "
                    f"with an eTreg-high or WT_heat-high one directly. The tables give "
                    f"the cell-state, tissue, and donor composition of the selection "
                    f"against the full shown set. This is a descriptive read over "
                    f"frozen per-cell scores, a secondary annotation tier held apart "
                    f"from the donor-level pseudobulk enrichment."
                ),
                mo.ui.plotly(_fig),
                mo.hstack(
                    [mo.md(_ct_tbl), mo.md(_ts_tbl), mo.md(_dn_tbl)],
                    justify="start",
                    gap=1.5,
                ),
            ]
        )
    _summary
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Hypoxia and UPR readouts — SF vs PB by cell state

        With labels frozen I harvested two Hallmark AUCell readouts on the same cells
        and summarized each as a SF-vs-PB standardized mean difference per sorted
        population. The bars below give that SMD. A positive value means the SF cells
        carry the higher readout. This is a secondary annotation-tier read, correlative
        throughout, and the hypoxia surface is a readout of the tissue state — not a
        claim that HIF or fever drives it.
        """
    )
    return


@app.cell
def _(go, mo, readout):
    _states = ["Treg", "Tcon", "CD8"]
    _reads = {
        "HALLMARK_HYPOXIA": ("hypoxia", "#c0392b"),
        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": ("UPR", "#2c7fb8"),
    }
    # One SF-vs-PB SMD per (cell state, readout). The summary carries the paired value
    # on both tissue rows, so I read it off the SF row per population.
    _sub = readout[readout["tissue"] == "synovial_fluid"]

    _fig = go.Figure()
    for _key, (_lab, _color) in _reads.items():
        _y = []
        for _st in _states:
            _m = _sub[(_sub["coarse_label"] == _st) & (_sub["readout"] == _key)]
            _y.append(float(_m["sf_minus_pb_smd"].iloc[0]) if len(_m) else 0.0)
        _fig.add_trace(
            go.Bar(
                x=_states,
                y=_y,
                name=_lab,
                marker_color=_color,
                text=[f"{v:+.2f}" for v in _y],
                textposition="outside",
                textfont=dict(size=14),
            )
        )
    _fig.update_layout(
        barmode="group",
        height=420,
        template="plotly_white",
        title=dict(text="SF-vs-PB readout SMD by cell state", x=0.5),
        yaxis_title="SF − PB standardized mean difference",
        xaxis_title="sorted population",
        legend=dict(orientation="h", y=1.12, x=0),
        margin=dict(l=10, r=10, t=60, b=10),
        font=dict(size=14),
    )

    # The hypoxia high-pocket: the fraction of cells above the global 90th-percentile
    # of the hypoxia AUCell surface, SF versus PB, per population.
    _hyp = readout[readout["readout"] == "HALLMARK_HYPOXIA"]
    _prows = []
    for _st in _states:
        _pb = _hyp[
            (_hyp["coarse_label"] == _st) & (_hyp["tissue"] == "peripheral_blood")
        ]
        _sf = _hyp[
            (_hyp["coarse_label"] == _st) & (_hyp["tissue"] == "synovial_fluid")
        ]
        if len(_pb) and len(_sf):
            _prows.append(
                f"| {_st} | {100 * float(_pb['frac_above_p90'].iloc[0]):.1f}% | "
                f"{100 * float(_sf['frac_above_p90'].iloc[0]):.1f}% |"
            )
    _ptbl = (
        "| population | PB in hypoxia-high pocket | SF in hypoxia-high pocket |\n"
        "|---|---|---|\n" + "\n".join(_prows)
    )

    _cap = mo.md(
        "**The plain reading.** Hypoxia runs SF-high in every population and **leads "
        "in Treg** (SF−PB SMD +1.67, versus +1.37 in Tcon and +0.83 in CD8). The "
        "high-pocket table sharpens it: **29.8%** of synovial-fluid Tregs sit in the "
        "top-decile hypoxia pocket against **1.7%** of blood Tregs, the widest "
        "SF-over-PB pocket enrichment of the three populations. The UPR readout is "
        "**flat**: near zero in Treg (+0.16) and mildly negative in Tcon (−0.22) and "
        "CD8 (−0.56), so it carries no coherent SF-vs-PB shift. Read correlatively, a "
        "synovial-Treg hypoxia pocket is the one readout worth carrying forward, and "
        "the hypoxia surface remains a readout of the tissue state, not a HIF-causality "
        "claim. Anything worth claiming from it is re-derived in the donor-level "
        "pseudobulk first.\n\n" + _ptbl
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## Harvest context — the current, revisable OR-gated drafted subset

        Alongside the primary claim I draft a bounded populations-of-interest subset
        off the same frozen labels, for later cross-dataset and cross-species
        questions. Membership is an OR-union of anchor-orthogonal hooks — a lineage
        hook, an effector hook, and a viable mitochondrial-high hook — under a
        viability gate. The mouse `WT_heat` score is carried as annotation only and
        is **never a selection predicate**. Everything here is a **current, revisable
        harvest**, open to revision in light of the hypoxia readout above, and the
        latent maps it feeds are maps, never evidence.
        """
    )
    return


@app.cell
def _(go, hook_lineage, mo, or_union):
    # OR-union membership breakdown across the whole compartment.
    _cat_order = [
        "baseline (not in union)",
        "lineage only",
        "effector only",
        "lineage + effector",
        "mt-hi viable pocket",
    ]
    _colors = {
        "baseline (not in union)": "#bdc3c7",
        "lineage only": "#2c7fb8",
        "effector only": "#27ae60",
        "lineage + effector": "#8e44ad",
        "mt-hi viable pocket": "#c0392b",
    }
    _um = or_union.set_index("membership_category")
    _labels = [c for c in _cat_order if c in _um.index]
    _vals = [float(_um.loc[c, "frac_all_cells"]) * 100 for c in _labels]
    _union_pct = float(or_union["union_frac_all_cells"].iloc[0]) * 100

    _fig = go.Figure(
        go.Bar(
            x=_labels,
            y=_vals,
            marker_color=[_colors[c] for c in _labels],
            text=[f"{v:.1f}%" for v in _vals],
            textposition="outside",
            textfont=dict(size=14),
        )
    )
    _fig.update_layout(
        height=420,
        template="plotly_white",
        title=dict(
            text=f"OR-union membership (drafted subset ≈ {_union_pct:.0f}% of compartment)",
            x=0.5,
        ),
        yaxis_title="% of all cells",
        xaxis_title="membership category",
        margin=dict(l=10, r=10, t=60, b=10),
        font=dict(size=13),
    )

    # Per-lineage hook fractions.
    _hrows = []
    for _, _r in hook_lineage.iterrows():
        _hrows.append(
            f"| {_r['lineage']} | {int(_r['n_cells']):,} | "
            f"{100 * float(_r['frac_hook_lineage']):.1f}% | "
            f"{100 * float(_r['frac_hook_effector']):.1f}% | "
            f"{100 * float(_r['frac_hook_mthi_viable']):.1f}% | "
            f"{100 * float(_r['frac_hook_or_union']):.1f}% |"
        )
    _htbl = (
        "| lineage | cells | lineage hook | effector hook | mt-hi viable hook | "
        "OR-union |\n|---|---|---|---|---|---|\n" + "\n".join(_hrows)
    )

    _cap = mo.md(
        f"**Reading the harvest.** The current OR-gated draft takes about "
        f"**{_union_pct:.0f}%** of the compartment. The lineage hook carries the bulk "
        f"of it (all Tregs by construction), the effector hook adds an eTreg-high "
        f"minority across all three populations, and the viable mitochondrial-high "
        f"hook is a small bounded pocket. Each hi arm keeps its matched-lo baseline so "
        f"downstream contrasts stay factorial. This draft is a **current, revisable** "
        f"design preview, not a committed cohort and not statistical evidence — the "
        f"hypoxia readout above is exactly the kind of signal that may reshape which "
        f"hooks I keep.\n\n" + _htbl
    )
    mo.vstack([mo.ui.plotly(_fig), _cap])
    return


@app.cell
def _(go_decision, mo):
    import json as _json

    _blob = _json.dumps(go_decision, indent=2, default=str)
    mo.md(
        f"""
        ## Decision

        **Where the analysis goes next: continue.** The mouse 39 °C `WT_heat`
        up-program is enriched in JIA synovial-fluid T cells relative to paired blood,
        and the enrichment is **broad across the sorted populations (pan-T), not
        Treg-preferential** — NES ≈ 2.51 in Treg, 2.57 in Tcon, 2.05 in CD8, all
        clearing FDR far below 1e-6. This is a positive, reproducible signal read as
        **consistent with** the mouse stress axis, correlative and not causal, and I
        state plainly that it is not Treg-specific.

        **Harvest strategy remains revisitable.** The OR-gated drafted subset above is
        a current design preview, not a committed cohort. The synovial-Treg hypoxia
        high-pocket surfaced in this review is the kind of orthogonal readout that may
        reshape which hooks I keep. That harvest call is taken separately later. This
        review makes no new claim beyond the pan-T enrichment, and the hypoxia and UPR
        readouts stay a secondary annotation tier, never pooled with the pseudobulk NES.

        Machine-readable echo of `decisions.go_no_go` (carried verbatim from config,
        unchanged here):

        ```json
        {_blob}
        ```
        """
    )
    return


if __name__ == "__main__":
    app.run()
