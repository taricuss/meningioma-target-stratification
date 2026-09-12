from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns

from .config import PLOT_STYLE, PATHS
from .stats import stars_for_p

logger = logging.getLogger(__name__)

sns.set_theme(
    context=PLOT_STYLE.get("context", "talk"),
    style=PLOT_STYLE.get("style", "whitegrid"),
    palette=PLOT_STYLE.get("palette", "viridis"),
    font_scale=float(PLOT_STYLE.get("font_scale", 1.1)),
)

DPI = int(PLOT_STYLE.get("figure_dpi", 300))
FORMATS = PLOT_STYLE.get("figure_format", ["png", "pdf"])


def _save_fig(fig: plt.Figure, stem: str, out_dir: Optional[Path] = None) -> List[Path]:
    if out_dir is None:
        out_dir = PATHS.results_figures
    out_dir.mkdir(parents=True, exist_ok=True)
    saved = []
    for fmt in FORMATS:
        p = out_dir / f"{stem}.{fmt}"
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)
    plt.close(fig)
    return saved


def fig1_study_design(out_dir: Optional[Path] = None) -> List[Path]:
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 6)
    ax.axis("off")
    boxes = [
        (0.5, 4.0, 3.0, 1.5, "Anchor papers\nCCR 2023 + StM 2026\n5 compounds / 4 mechanisms", "#66c2a5"),
        (4.0, 4.0, 3.0, 1.5, "Two independent classification systems\nNassiri 2021 + Bi lab 2023\nDiscovery cohorts", "#8da0cb"),
        (7.5, 4.0, 3.0, 1.5, "Aims 0-3\nCohort lock · Concordance ·\nTarget-gene scoring + replication", "#fc8d62"),
        (11.0, 4.0, 3.0, 1.5, "H5 Optional survival\nexploratory", "#e78ac3"),
        (0.5, 1.0, 13.5, 2.0, "Deliverable: testable patient-selection hypotheses\nreturned to Warta/Heidelberg for validation", "#a6d854"),
    ]
    for x, y, w, h, txt, c in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", facecolor=c, edgecolor="k", lw=1.2, alpha=0.85)
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=10, weight="bold", color="black")
    arrows = [((3.5, 4.75), (4.0, 4.75)), ((7.0, 4.75), (7.5, 4.75)), ((10.5, 4.75), (11.0, 4.75)), ((7.25, 4.0), (7.25, 3.0))]
    for (x1, y1), (x2, y2) in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=2))
    ax.set_title("Fig. 1 — Study design: bridging Heidelberg drug screens to subgroup-stratified patient-selection hypotheses", pad=15, fontsize=14, weight="bold")
    return _save_fig(fig, "fig1_study_design", out_dir)


def fig2_subgroup_target_boxplots(
    program_scores_z: pd.DataFrame,
    metadata: pd.DataFrame,
    subgroup_col: str = "nassiri_group",
    subgroup_title: str = "Nassiri 2021 molecular groups",
    out_dir: Optional[Path] = None,
    subgroup_order: Optional[List[str]] = None,
) -> List[Path]:
    programs = program_scores_z.columns.tolist()
    n_prog = len(programs)
    ncols = 2
    nrows = int(np.ceil(n_prog / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows), sharex=False)
    axes = np.atleast_1d(axes).flatten()
    for idx, prog in enumerate(programs):
        ax = axes[idx]
        df = pd.DataFrame({"score": program_scores_z[prog], "group": metadata[subgroup_col]}).dropna()
        if subgroup_order:
            df["group"] = pd.Categorical(df["group"], categories=subgroup_order, ordered=True)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            sns.violinplot(
                data=df, x="group", y="score", hue="group", ax=ax,
                inner="quartile", linewidth=1.2, cut=0, palette="Set2", legend=False,
            )
            sns.boxplot(
                data=df, x="group", y="score", hue="group", ax=ax,
                width=0.15, showcaps=False, legend=False,
                boxprops=dict(facecolor="white", alpha=0.5),
                whiskerprops=dict(alpha=0.0),
                showfliers=False,
            )
        ax.set_xlabel("")
        ax.set_ylabel("ssGSEA z-score")
        ax.set_title(prog.replace("_", "\n"), fontsize=11, weight="bold")
        ax.tick_params(axis="x", rotation=20)
    for j in range(n_prog, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f"Fig. 2 — Target-gene program scores by {subgroup_title}", fontsize=15, weight="bold", y=1.01)
    fig.tight_layout()
    return _save_fig(fig, f"fig2_{subgroup_col}_target_boxplots", out_dir)


def fig3_nf2_forest(
    h3_result: pd.DataFrame,
    title: str = "H3 · NF2 Mutant/Loss vs Intact — per-program effect",
    out_dir: Optional[Path] = None,
) -> List[Path]:
    df = h3_result.copy()
    if df.empty:
        return []
    fig, ax = plt.subplots(figsize=(11, 0.5 + 0.6 * len(df)))
    y = np.arange(len(df))[::-1]
    x = df["delta_median_a_minus_b"].values
    p = df["p_bh_fdr_family"].values if "p_bh_fdr_family" in df else df["p_raw"].values
    cliffs = df["cliffs_delta"].values
    ax.errorbar(x, y, xerr=0.0, fmt="none", ecolor="gray", alpha=0.0)
    colors = ["#d7191c" if pv < 0.05 else "#2c7bb6" for pv in p]
    ax.scatter(x, y, s=140, c=colors, edgecolors="k", linewidths=0.8, zorder=3)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.program}\nΔ={r.delta_median_a_minus_b:+.2f} | Cliff's δ={cliffs_delta:.2f} [{stars_for_p(pv)}]"
                        for r, pv, cliffs_delta in zip(df.itertuples(), p, cliffs)])
    ax.set_xlabel("Δ median z-score (Intact − Mutant/Loss)")
    ax.set_title(title, fontsize=13, weight="bold")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    return _save_fig(fig, "fig3_nf2_forest", out_dir)


def fig4_replication_forest(
    replication_effects: pd.DataFrame,
    out_dir: Optional[Path] = None,
) -> List[Path]:
    if replication_effects.empty:
        return []
    df = replication_effects.copy()
    y = np.arange(len(df))[::-1]
    fig, (ax_disc, ax_rep) = plt.subplots(1, 2, figsize=(13, 0.5 + 0.5 * len(df)), sharey=True)
    for ax, col, t in [
        (ax_disc, "delta_discovery", "Discovery"),
        (ax_rep, "delta_replication", "Replication"),
    ]:
        x = df[col].values
        same = df["same_direction"].values if "same_direction" in df else [True] * len(df)
        colors = ["#1a9641" if s else "#ca0020" for s in same]
        ax.scatter(x, y, s=110, c=colors, edgecolors="k", linewidths=0.8)
        ax.axvline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_title(t, weight="bold")
        ax.set_xlabel("Δ median z-score")
        ax.grid(True, axis="x", alpha=0.3)
    ax_disc.set_yticks(y)
    ax_disc.set_yticklabels(df.apply(lambda r: f"{r['analysis'][:30]}: {r['program']}", axis=1).values)
    fig.suptitle("Fig. 4 — H4 · External replication: discovery vs independent cohort effect directions",
                 fontsize=14, weight="bold", y=1.01)
    fig.tight_layout()
    return _save_fig(fig, "fig4_replication_forest", out_dir)


def fig5_hdac8_resistance_axis(
    program_scores_z: pd.DataFrame,
    metadata: pd.DataFrame,
    subgroup_col: str = "nassiri_group",
    out_dir: Optional[Path] = None,
) -> List[Path]:
    if "HDAC8_TGFb_EMT_Resistance" not in program_scores_z.columns:
        return []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    df = pd.DataFrame({
        "resistance": program_scores_z["HDAC8_TGFb_EMT_Resistance"],
        "hdac": program_scores_z["HDAC_Panobinostat_Romidepsin"],
        "group": metadata[subgroup_col],
    }).dropna()
    ax = axes[0]
    sns.scatterplot(data=df, x="hdac", y="resistance", hue="group", style="group", palette="Set2", ax=ax, s=80, alpha=0.85)
    for g in df["group"].unique():
        sub = df[df["group"] == g]
        if len(sub) > 5:
            z = np.polyfit(sub["hdac"], sub["resistance"], 1)
            xlin = np.linspace(sub["hdac"].min(), sub["hdac"].max(), 50)
            ax.plot(xlin, np.polyval(z, xlin), linewidth=2, alpha=0.5)
    ax.set_xlabel("HDAC1/2 (panobinostat target) z-score")
    ax.set_ylabel("HDAC8/TGFβ/EMT resistance z-score")
    ax.set_title("StM 2026 resistance hook — target vs axis", weight="bold")

    ax = axes[1]
    order = sorted(df["group"].unique().tolist())
    data_to_plot = [df.loc[df["group"] == g, "resistance"].values for g in order]
    bp = ax.boxplot(data_to_plot, labels=order, patch_artist=True, widths=0.55, vert=True, showfliers=False)
    cmap = plt.get_cmap("Set2")
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(cmap(i % 8))
    for i, (g, vals) in enumerate(zip(order, data_to_plot)):
        from scipy.stats import kruskal
    ax.set_ylabel("HDAC8/TGFβ/EMT z-score")
    ax.set_title("Resistance-axis expression by subgroup", weight="bold")
    ax.tick_params(axis="x", rotation=20)
    fig.suptitle("Fig. 5 — Panobinostat resistance axis (StM 2026 HDAC8→TGFβ→EMT) across subgroups",
                 fontsize=13, weight="bold", y=1.02)
    fig.tight_layout()
    return _save_fig(fig, "fig5_hdac8_resistance_axis", out_dir)


def heatmap_concordance(cross: pd.DataFrame, out_dir: Optional[Path] = None, title: str = "H1 Concordance") -> List[Path]:
    cross_clean = cross.drop(index=["Total"], columns=["Total"], errors="ignore")
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * cross_clean.shape[1] + 3), max(5, 0.5 * cross_clean.shape[0] + 2)))
    pct = cross_clean.div(cross_clean.sum(axis=1), axis=0).round(3) * 100
    sns.heatmap(pct, annot=True, fmt=".1f", cmap="Blues", cbar_kws={"label": "Row %"}, linewidths=1.0, ax=ax)
    ax.set_xlabel("Bi-lab group")
    ax.set_ylabel("Nassiri group")
    ax.set_title(title, weight="bold")
    fig.tight_layout()
    return _save_fig(fig, "fig_s1_concordance_heatmap", out_dir)


def all_figures(ctx: dict, out_dir: Optional[Path] = None) -> Dict[str, List[Path]]:
    out = {}
    out["fig1_study_design"] = fig1_study_design(out_dir)
    if "aim2" in ctx and "nassiri" in ctx["aim2"]:
        r = ctx["aim2"]["nassiri"]
        meta = ctx["metadata_discovery"]
        out["fig2a_nassiri_boxplots"] = fig2_subgroup_target_boxplots(
            r["program_scores_z"], meta, subgroup_col="nassiri_group",
            subgroup_title="Nassiri 2021 groups", out_dir=out_dir,
            subgroup_order=["Immunogenic", "NF2-inactivated canonical", "Hypermetabolic"],
        )
        out["fig5_resistance"] = fig5_hdac8_resistance_axis(
            r["program_scores_z"], meta, subgroup_col="nassiri_group", out_dir=out_dir,
        )
        out["fig3_nf2"] = fig3_nf2_forest(r["h3_nf2"], out_dir=out_dir)
    if "aim2" in ctx and "bi" in ctx["aim2"]:
        r = ctx["aim2"]["bi"]
        meta = ctx["metadata_discovery_bi"]
        out["fig2b_bi_boxplots"] = fig2_subgroup_target_boxplots(
            r["program_scores_z"], meta, subgroup_col="bi_group",
            subgroup_title="Bi-lab 3-tier groups", out_dir=out_dir,
            subgroup_order=["Merlin-intact", "Immune-enriched", "Hypermitotic"],
        )
    if "concordance" in ctx:
        out["fig_s1_concordance"] = heatmap_concordance(ctx["concordance"].cross_tabulation, out_dir=out_dir)
    if "replication" in ctx:
        out["fig4_replication"] = fig4_replication_forest(ctx["replication"]["per_effect_forest"], out_dir=out_dir)
    return out
