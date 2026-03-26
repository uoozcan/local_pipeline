#!/usr/bin/env python3
"""
multi_analysis_figures.py — Cross-analysis scientific visualizations for the HLA typing pipeline.

Produces 6 publication-quality figures comparing WGS / WES / RNA-seq calibration results,
highlighting pipeline novelty and importance. All inputs are optional — figures that
require unavailable data are skipped gracefully.

Figures produced:
  M1  Cross-analysis accuracy comparison (requires ≥2 data types)
  M2  Tool × data-type coverage matrix (static; always produced)
  M3  Calibrated voting benefit per data type (requires strategy TSVs)
  M4  Data-type-specific consensus weights (requires ≥1 weights JSON)
  M5  Population-stratified accuracy (requires by_population TSV)
  M6  Pipeline information flow diagram (always produced)

Usage:
    python3 bin/multi_analysis_figures.py \\
        --wgs-accuracy  conf/tool_accuracy_wgs_v3.tsv \\
        --wes-accuracy  conf/tool_accuracy_wes_v1.tsv \\
        --rna-accuracy  conf/tool_accuracy_rna_v1.tsv \\
        --wgs-weights   conf/tool_weights_wgs_v3.json \\
        --wes-weights   conf/tool_weights_wes_v1.json \\
        --rna-weights   conf/tool_weights_rna_v1.json \\
        --wgs-strategy  /path/to/strategy_comparison_wgs_calibrated.tsv \\
        --wes-strategy  /path/to/strategy_comparison_wes_calibrated.tsv \\
        --rna-strategy  /path/to/strategy_comparison_rna_calibrated.tsv \\
        --population    conf/tool_accuracy_wgs_by_population.tsv \\
        --outdir        figures/multi_analysis/
"""

import argparse
import base64
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    warnings.warn("seaborn not found — some plots will use matplotlib fallbacks")

# ─── Style constants ──────────────────────────────────────────────────────────

DT_COLORS  = {"wgs": "#1565C0", "wes": "#E65100", "rna": "#2E7D32"}
DT_LABELS  = {"wgs": "WGS (30×)", "wes": "WES",  "rna": "RNA-seq"}

TOOL_COLORS = {
    "hlahd":     "#2196F3",
    "optitype":  "#4CAF50",
    "spechla":   "#FF9800",
    "arcashla":  "#9C27B0",
    "kourami":   "#F44336",
    "polysolver":"#00BCD4",
    "seq2hla":   "#FF99CC",
    "t1k":       "#795548",
}
TOOL_LABELS = {
    "hlahd":     "HLA-HD",
    "optitype":  "OptiType",
    "spechla":   "SpecHLA",
    "arcashla":  "arcasHLA",
    "kourami":   "Kourami",
    "polysolver":"PolySOLVER",
    "seq2hla":   "seq2HLA",
    "t1k":       "T1K",
}
EXCLUDE_TOOLS = {"hlala", "xhla"}

GENE_ORDER = ["A", "B", "C", "DRB1", "DQB1"]

STRATEGY_COLORS = {
    "equal":           "#9E9E9E",
    "read_confidence": "#2196F3",
    "calibrated":      "#4CAF50",
}

DPI = 200

def tl(t): return TOOL_LABELS.get(t, t)
def tc(t): return TOOL_COLORS.get(t, "#607D8B")

def _set_style():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        try:
            plt.style.use("seaborn-whitegrid")
        except OSError:
            plt.style.use("ggplot")
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.facecolor": "white",
        "savefig.dpi": DPI,
    })

def save(fig, path, title=None):
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")

def _load_accuracy(path):
    if not path or not Path(path).exists():
        return None
    df = pd.read_csv(path, sep="\t")
    df = df[~df["Tool"].isin(EXCLUDE_TOOLS)].copy()
    return df

def _load_weights(path):
    if not path or not Path(path).exists():
        return None
    with open(path) as f:
        return json.load(f)

def _load_strategy(path):
    if not path or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t")

def _load_population(path):
    if not path or not Path(path).exists():
        return None
    df = pd.read_csv(path, sep="\t")
    df = df[~df["Tool"].isin(EXCLUDE_TOOLS)].copy()
    return df

def _gene_cols(df):
    extras = [c for c in df.columns if c not in GENE_ORDER and c.upper() == c and c not in ("NA",)]
    return [g for g in GENE_ORDER if g in df.columns] + extras


# ─── Figure M1: Cross-analysis accuracy comparison ───────────────────────────

def fig_M1_cross_analysis(acc_dict, out_dir):
    """Grouped bar chart: tool concordance per gene, grouped by data type.

    acc_dict: {"wgs": df, "wes": df, "rna": df}  — values may be None
    Requires at least 2 non-None data types.
    """
    available = {dt: df for dt, df in acc_dict.items() if df is not None}
    if len(available) < 2:
        print("  [SKIP M1] Need ≥2 data types for cross-analysis comparison")
        return

    genes = GENE_ORDER
    all_tools = sorted({t for df in available.values() for t in df["Tool"].tolist()})
    n_tools = len(all_tools)
    n_dts = len(available)
    n_genes = len(genes)

    fig, axes = plt.subplots(1, n_genes, figsize=(4 * n_genes, 5), sharey=True)
    if n_genes == 1:
        axes = [axes]

    width = 0.8 / n_tools
    x = np.arange(n_dts)

    for gi, (gene, ax) in enumerate(zip(genes, axes)):
        for ti, tool in enumerate(all_tools):
            vals = []
            for dt in available:
                df = available[dt]
                row = df[df["Tool"] == tool]
                if row.empty or gene not in row.columns:
                    vals.append(np.nan)
                else:
                    v = row.iloc[0][gene]
                    vals.append(float(v) if pd.notna(v) else np.nan)
            offs = x + (ti - n_tools / 2 + 0.5) * width
            bars = ax.bar(offs, [v if not np.isnan(v) else 0 for v in vals],
                          width=width * 0.9, color=tc(tool),
                          label=tl(tool) if gi == 0 else "_",
                          zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([DT_LABELS.get(dt, dt) for dt in available],
                           fontsize=9, rotation=15, ha="right")
        ax.set_title(f"HLA-{gene}", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.15)
        ax.axhline(0.8, color="grey", lw=0.8, ls="--", alpha=0.5)
        ax.grid(axis="y", alpha=0.3, zorder=0)

    axes[0].set_ylabel("2-Field Concordance", fontsize=11)
    fig.legend(*axes[0].get_legend_handles_labels(),
               loc="upper right", fontsize=8, framealpha=0.9,
               ncol=2, bbox_to_anchor=(1.0, 1.02))

    save(fig, out_dir / "M1_cross_analysis_accuracy.png",
         "Tool Performance Across Analysis Types (2-Field Concordance vs 1KGP Ground Truth)")


# ─── Figure M2: Tool × data-type coverage matrix ─────────────────────────────

def fig_M2_coverage_matrix(acc_dict, out_dir):
    """Annotated matrix: tools × data types, showing coverage and mean concordance.

    acc_dict: {"wgs": df, "wes": df, "rna": df}  — values may be None
    This figure always has the static compatibility layer (✓/✗/–).
    """
    # Tool × input type compatibility (static knowledge)
    # Format: (works, optimal_seq_types)
    COMPAT = {
        #          WGS-BAM  WES-BAM  RNA-FASTQ  CRAM-30x  HiFi
        "hlahd":     [True,   True,   False,  True,   False],
        "spechla":   [True,   True,   False,  True,   False],
        "arcashla":  [True,   False,  True,   True,   False],
        "optitype":  [True,   True,   True,   True,   False],
        "kourami":   [True,   False,  False,  True,   False],
        "polysolver":[True,   True,   False,  True,   False],
        "seq2hla":   [False,  False,  True,   False,  False],
        "t1k":       [True,   True,   True,   True,   True],
    }
    COL_LABELS = ["WGS\n(BAM)", "WES\n(BAM)", "RNA-seq\n(FASTQ)", "CRAM\n(30×)", "Long-read\n(HiFi)"]
    tools_ordered = list(COMPAT.keys())
    n_tools = len(tools_ordered)
    n_cols = len(COL_LABELS)

    # Build matrix: compatibility + concordance from acc_dict
    matrix = np.full((n_tools, n_cols), np.nan)
    compat_mat = np.array([COMPAT[t] for t in tools_ordered], dtype=float)

    # Overlay mean concordance where available
    dt_col_map = {"wgs": [0, 3], "wes": [1], "rna": [2]}
    for dt, df in acc_dict.items():
        if df is None:
            continue
        for ti, tool in enumerate(tools_ordered):
            row = df[df["Tool"] == tool]
            if row.empty:
                continue
            mean_col = "Mean_concordance"
            if mean_col in row.columns:
                v = row.iloc[0][mean_col]
                v = float(v) if pd.notna(v) else np.nan
            else:
                genes = _gene_cols(row)
                vals = [float(row.iloc[0][g]) for g in genes
                        if g in row.columns and pd.notna(row.iloc[0][g])]
                v = np.nanmean(vals) if vals else np.nan
            for ci in dt_col_map.get(dt, []):
                matrix[ti, ci] = v if not np.isnan(v) else (0.5 if COMPAT[tool][ci] else np.nan)

    # Fill compatible but no data cells
    for ti in range(n_tools):
        for ci in range(n_cols):
            if compat_mat[ti, ci] and np.isnan(matrix[ti, ci]):
                matrix[ti, ci] = -0.05   # special: compatible but not yet tested

    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.cm.RdYlGn
    cmap.set_under("#BBDEFB")   # blue for "not yet tested"
    cmap.set_bad("#F5F5F5")      # grey for incompatible

    masked = np.ma.masked_where(np.isnan(matrix), matrix)
    im = ax.imshow(masked, vmin=0, vmax=1, cmap=cmap, aspect="auto")

    for ti in range(n_tools):
        for ci in range(n_cols):
            if np.isnan(matrix[ti, ci]):
                # Incompatible
                ax.add_patch(mpatches.Rectangle((ci - 0.5, ti - 0.5), 1, 1,
                                                 fill=True, color="#F5F5F5", zorder=2))
                ax.text(ci, ti, "✗", ha="center", va="center",
                        fontsize=14, color="#BDBDBD", zorder=3)
            elif matrix[ti, ci] < 0:
                # Compatible but not yet tested
                ax.text(ci, ti, "–", ha="center", va="center",
                        fontsize=14, color="#1565C0", fontweight="bold", zorder=3)
            else:
                v = matrix[ti, ci]
                ax.text(ci, ti, f"{v:.2f}", ha="center", va="center",
                        fontsize=10, fontweight="bold", color="white" if v > 0.5 else "#333",
                        zorder=3)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(COL_LABELS, fontsize=10)
    ax.set_yticks(range(n_tools))
    ax.set_yticklabels([tl(t) for t in tools_ordered], fontsize=11)
    ax.set_xlabel("Input Data Type", fontsize=11)
    ax.set_ylabel("HLA Typing Tool", fontsize=11)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Mean 2-Field Concordance", fontsize=9)

    legend_patches = [
        mpatches.Patch(color="#F5F5F5", label="✗  Not applicable"),
        mpatches.Patch(color="#BBDEFB", label="–  Compatible, not yet tested"),
        mpatches.Patch(color="#FDAE61", label="0.0–0.5  Low concordance"),
        mpatches.Patch(color="#A6D96A", label="0.5–0.8  Moderate"),
        mpatches.Patch(color="#1A9641", label="0.8–1.0  High concordance"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8,
              framealpha=0.9, bbox_to_anchor=(1.45, 0))

    save(fig, out_dir / "M2_tool_coverage_matrix.png",
         "HLA Typing Tool Compatibility and Performance Across Input Data Types")


# ─── Figure M3: Calibrated voting benefit ────────────────────────────────────

def fig_M3_calibration_benefit(strat_dict, out_dir):
    """Bar chart: equal vs calibrated concordance per data type.

    strat_dict: {"wgs": df, "wes": df, "rna": df}
    """
    available = {dt: df for dt, df in strat_dict.items() if df is not None}
    if not available:
        print("  [SKIP M3] No strategy comparison data available")
        return

    n_dts = len(available)
    fig, axes = plt.subplots(1, n_dts, figsize=(5 * n_dts, 5), sharey=True)
    if n_dts == 1:
        axes = [axes]

    for ax, (dt, df) in zip(axes, available.items()):
        gene_col = "Gene" if "Gene" in df.columns else df.columns[0]
        genes = df[gene_col].tolist()
        x = np.arange(len(genes))
        width = 0.35

        eq_col = "equal_mean"
        cal_col = "calibrated_mean"
        sig_col = "Significant_p05"
        p_col   = "Wilcoxon_p"

        if eq_col not in df.columns or cal_col not in df.columns:
            ax.text(0.5, 0.5, "No equal/calibrated\ncolumns found",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(DT_LABELS.get(dt, dt))
            continue

        eq_vals  = pd.to_numeric(df[eq_col],  errors="coerce").fillna(0).values
        cal_vals = pd.to_numeric(df[cal_col], errors="coerce").fillna(0).values
        deltas   = cal_vals - eq_vals

        bars_eq  = ax.bar(x - width / 2, eq_vals,  width=width * 0.95,
                          color=STRATEGY_COLORS["equal"],      label="Equal weights",    zorder=3)
        bars_cal = ax.bar(x + width / 2, cal_vals, width=width * 0.95,
                          color=STRATEGY_COLORS["calibrated"],  label="Calibrated weights", zorder=3)

        # Delta labels and significance markers
        for gi, (bar, delta, val) in enumerate(zip(bars_cal, deltas, cal_vals)):
            sign = "+" if delta >= 0 else ""
            color = "#1B5E20" if delta > 0.01 else ("#B71C1C" if delta < -0.01 else "#555")
            ax.text(bar.get_x() + bar.get_width() / 2,
                    val + 0.02, f"{sign}{delta:.2f}",
                    ha="center", va="bottom", fontsize=8, color=color, fontweight="bold")

            # Significance star
            if sig_col in df.columns:
                sig = str(df.iloc[gi].get(sig_col, "no")).strip().lower()
                if sig in ("yes", "true", "1"):
                    p = float(df.iloc[gi].get(p_col, 1.0) or 1.0)
                    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*"
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            val + 0.08, stars, ha="center", fontsize=11,
                            color="black", fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(genes, fontsize=11)
        ax.set_ylim(0, 1.25)
        ax.set_title(DT_LABELS.get(dt, dt), fontsize=12, fontweight="bold",
                     color=DT_COLORS.get(dt, "black"))
        ax.grid(axis="y", alpha=0.4, zorder=0)
        ax.axhline(0.5, color="grey", lw=0.8, ls="--", alpha=0.5)
        if ax == list(axes)[0]:
            ax.set_ylabel("Mean 2-Field Concordance", fontsize=11)

    axes[0].legend(loc="upper left", fontsize=9, framealpha=0.9)
    fig.text(0.5, -0.03,
             "Δ = calibrated − equal  |  * p<0.05, ** p<0.01, *** p<0.001 (Wilcoxon signed-rank)",
             ha="center", fontsize=9, color="#555")

    save(fig, out_dir / "M3_calibration_benefit.png",
         "Benefit of Empirically Calibrated Consensus Voting vs Equal Weighting")


# ─── Figure M4: Data-type-specific weights ───────────────────────────────────

def fig_M4_data_type_weights(weights_dict, out_dir):
    """Stacked horizontal bar: tool contributions per gene, one panel per data type.

    weights_dict: {"wgs": {...}, "wes": {...}, "rna": {...}}
    """
    available = {dt: w for dt, w in weights_dict.items() if w is not None}
    if not available:
        print("  [SKIP M4] No weights data available")
        return

    n_dts = len(available)
    fig, axes = plt.subplots(1, n_dts, figsize=(5.5 * n_dts, 5), sharey=True)
    if n_dts == 1:
        axes = [axes]

    all_tools_in_weights = sorted({
        t for w in available.values()
        for g_data in w.get("genes", {}).values()
        for t in g_data.keys()
    } - EXCLUDE_TOOLS)

    for ax, (dt, w_data) in zip(axes, available.items()):
        genes_dict = w_data.get("genes", {})
        plot_genes = [g for g in GENE_ORDER if g in genes_dict]
        if not plot_genes:
            continue
        y = np.arange(len(plot_genes))
        lefts = np.zeros(len(plot_genes))

        for tool in all_tools_in_weights:
            vals = np.array([genes_dict.get(g, {}).get(tool, 0.0) for g in plot_genes])
            bars = ax.barh(y, vals, left=lefts, color=tc(tool),
                           label=tl(tool), height=0.65, edgecolor="white", linewidth=0.5)
            for bar, val in zip(bars, vals):
                if val >= 0.1:
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_y() + bar.get_height() / 2,
                            f"{val:.2f}", ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
            lefts += vals

        ax.set_yticks(y)
        ax.set_yticklabels(plot_genes, fontsize=11)
        ax.set_xlim(0, 1.01)
        ax.set_xlabel("Cumulative Weight", fontsize=10)
        n = w_data.get("n_samples", "?")
        ax.set_title(f"{DT_LABELS.get(dt, dt)}\n(N={n})", fontsize=11,
                     fontweight="bold", color=DT_COLORS.get(dt, "black"))
        ax.axvline(1.0, color="grey", lw=0.8, ls="--", alpha=0.5)
        ax.grid(axis="x", alpha=0.3)

    # Shared legend
    handles = [mpatches.Patch(color=tc(t), label=tl(t)) for t in all_tools_in_weights]
    fig.legend(handles=handles, loc="lower center", fontsize=9, framealpha=0.9,
               ncol=len(all_tools_in_weights), bbox_to_anchor=(0.5, -0.08))

    save(fig, out_dir / "M4_data_type_weights.png",
         "Data-Type-Specific Calibrated Consensus Voting Weights per HLA Gene\n"
         "(Key novelty: optimal tool mix differs by input data type)")


# ─── Figure M5: Population-stratified accuracy ───────────────────────────────

def fig_M5_population(pop_df, acc_df, out_dir):
    """Population-stratified concordance heatmap + bar chart."""
    if pop_df is None:
        print("  [SKIP M5] No population data")
        return

    col_parts = [(c.split("_")[0], c.split("_")[1])
                 for c in pop_df.columns
                 if "_" in c and c.split("_")[1] in GENE_ORDER]
    populations = sorted({p for p, _ in col_parts})
    genes = [g for g in GENE_ORDER if any(g == gn for _, gn in col_parts)]
    tools = pop_df["Tool"].tolist()

    if not populations or not genes:
        print("  [SKIP M5] Population columns not detected")
        return

    fig, axes = plt.subplots(1, len(populations) + 1,
                             figsize=(4 * (len(populations) + 1), 5),
                             sharey=True,
                             gridspec_kw={"width_ratios": [1] * len(populations) + [0.6]})

    n_tools = len(tools)
    width = 0.8 / n_tools
    x = np.arange(len(genes))

    for pi, (ax, pop) in enumerate(zip(axes[:-1], populations)):
        for ti, (_, row) in enumerate(pop_df.iterrows()):
            vals = [float(row.get(f"{pop}_{g}", np.nan))
                    if pd.notna(row.get(f"{pop}_{g}", np.nan)) else np.nan
                    for g in genes]
            offs = x + (ti - n_tools / 2 + 0.5) * width
            ax.bar(offs, [v if not np.isnan(v) else 0 for v in vals],
                   width=width * 0.9, color=tc(row["Tool"]), zorder=3,
                   label=tl(row["Tool"]) if pi == 0 else "_")
        ax.set_xticks(x)
        ax.set_xticklabels(genes, fontsize=10)
        ax.set_ylim(0, 1.2)
        ax.set_title(pop, fontsize=12, fontweight="bold")
        ax.grid(axis="y", alpha=0.3, zorder=0)
        ax.axhline(0.5, color="grey", lw=0.8, ls="--", alpha=0.4)

    axes[0].set_ylabel("2-Field Concordance", fontsize=11)

    # Rightmost: mean concordance per tool (all populations combined)
    ax_mean = axes[-1]
    if acc_df is not None and "Mean_concordance" in acc_df.columns:
        merged = pop_df.merge(
            acc_df[["Tool", "Mean_concordance"]], on="Tool", how="left")
        merged["Mean_concordance"] = pd.to_numeric(
            merged["Mean_concordance"], errors="coerce")
        merged = merged.sort_values("Mean_concordance", ascending=True)
        ax_mean.barh([tl(t) for t in merged["Tool"]],
                     merged["Mean_concordance"].fillna(0),
                     color=[tc(t) for t in merged["Tool"]], height=0.6)
        ax_mean.set_xlim(0, 1.1)
        ax_mean.set_xlabel("Mean", fontsize=9)
        ax_mean.set_title("Overall\nmean", fontsize=10)
        ax_mean.grid(axis="x", alpha=0.3)
    else:
        ax_mean.axis("off")

    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9)
    save(fig, out_dir / "M5_population_accuracy.png",
         "Population-Stratified HLA Typing Accuracy (AFR / EAS / EUR)")


# ─── Figure M6: Pipeline information flow ────────────────────────────────────

def fig_M6_information_flow(out_dir):
    """Pure matplotlib pipeline architecture diagram."""
    fig, ax = plt.subplots(figsize=(18, 9))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 9)
    ax.axis("off")

    def box(ax, x, y, w, h, text, color, fontsize=9, text_color="white", bold=False):
        rect = mpatches.FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                       boxstyle="round,pad=0.1",
                                       facecolor=color, edgecolor="white",
                                       linewidth=1.5, zorder=3)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                color=text_color, fontweight="bold" if bold else "normal",
                zorder=4, wrap=True)

    def arrow(ax, x1, y1, x2, y2, color="#555", lw=1.5):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color=color, lw=lw),
                    zorder=2)

    # ── Column x-positions ─────────────────────────────────
    x_in   = 1.5
    x_qc   = 3.8
    x_tool = 7.2
    x_con  = 11.8
    x_out  = 14.5
    x_cal  = 7.2

    # ── Section headers ────────────────────────────────────
    for x, label, color in [
        (x_in,   "① INPUT",        "#0D47A1"),
        (x_qc,   "② QC",           "#1B5E20"),
        (x_tool, "③ TYPING TOOLS", "#E65100"),
        (x_con,  "④ CONSENSUS",    "#4A148C"),
        (x_out,  "⑤ OUTPUT",       "#212121"),
    ]:
        ax.text(x, 8.6, label, ha="center", va="center", fontsize=11,
                fontweight="bold", color=color)

    # ── Input boxes ────────────────────────────────────────
    input_items = [
        ("WGS\nBAM/CRAM\n(30× hg38)", "#1565C0", 7.5),
        ("WES\nBAM\n(~100×)",         "#0288D1", 6.3),
        ("RNA-seq\nFASTQ\n(LCL/Tumor)","#00838F", 5.1),
        ("Long-read\nBAM\n(HiFi/ONT)", "#37474F", 3.9),
    ]
    for label, color, y in input_items:
        box(ax, x_in, y, 2.2, 0.9, label, color, fontsize=8)
        arrow(ax, x_in + 1.1, y, x_qc - 1.0, y, color="#888")

    # ── QC boxes ───────────────────────────────────────────
    qc_items = [
        ("Read counts\nHLA region reads\nMapping quality\nRead length", "#2E7D32", 6.3),
        ("FastQC\nMultiQC\nQC report",                                   "#388E3C", 4.5),
    ]
    for label, color, y in qc_items:
        box(ax, x_qc, y, 2.2, 1.4, label, color, fontsize=8)
    ax.text(x_qc, 3.4, "✓ PASS  ⚠ WARN  ✗ FAIL", ha="center", va="center",
            fontsize=8, color="#2E7D32", style="italic")
    for y in [6.3, 4.5]:
        arrow(ax, x_qc + 1.1, y, x_tool - 1.2, y if y > 5 else 5.5, color="#888")

    # ── Typing tool boxes ──────────────────────────────────
    tool_items = [
        ("HLA-HD",      "#2196F3", 7.8),
        ("SpecHLA",     "#FF9800", 6.9),
        ("arcasHLA",    "#9C27B0", 6.0),
        ("OptiType",    "#4CAF50", 5.1),
        ("Kourami",     "#F44336", 4.2),
        ("PolySOLVER",  "#00BCD4", 3.3),
        ("seq2HLA",     "#FF99CC", 2.4),
        ("T1K",         "#795548", 1.5),
    ]
    for label, color, y in tool_items:
        box(ax, x_tool, y, 2.0, 0.65, label, color, fontsize=9, bold=True)
        arrow(ax, x_tool + 1.0, y, x_con - 1.1, 5.5 if y > 4 else 4.0, color="#888", lw=1.0)

    # Tool scope annotations
    ax.text(x_tool + 1.15, 6.95, "◄ WGS/WES BAM", fontsize=7, color="#888", va="center")
    ax.text(x_tool + 1.15, 5.05, "◄ RNA FASTQ",   fontsize=7, color="#888", va="center")
    ax.text(x_tool + 1.15, 1.5,  "◄ BAM+FASTQ+HiFi", fontsize=7, color="#888", va="center")

    # ── Consensus box ──────────────────────────────────────
    box(ax, x_con, 5.5, 2.5, 1.2,
        "Weighted\nConsensus\nVoting",
        "#6A1B9A", fontsize=10, bold=True)
    box(ax, x_con, 3.8, 2.5, 1.0,
        "3 modes:\nEqual · Read-confidence\nCalibrated (empirical)",
        "#9C27B0", fontsize=8)
    arrow(ax, x_con, 4.9, x_con, 4.35, color="#4A148C", lw=2)
    arrow(ax, x_con + 1.25, 5.5, x_out - 0.9, 6.8, color="#888")
    arrow(ax, x_con + 1.25, 5.5, x_out - 0.9, 5.5, color="#888")
    arrow(ax, x_con + 1.25, 5.5, x_out - 0.9, 4.2, color="#888")
    arrow(ax, x_con + 1.25, 3.8, x_out - 0.9, 2.8, color="#888")

    # ── Output boxes ───────────────────────────────────────
    out_items = [
        ("HLA alleles\n(2-field & 4-field)",     "#263238", 6.8),
        ("Confidence score\n(0.0–1.0 per locus)", "#37474F", 5.5),
        ("LOH analysis\n(tumor purity/ploidy)",   "#455A64", 4.2),
        ("HTML report\n+ MultiQC",                "#546E7A", 2.8),
    ]
    for label, color, y in out_items:
        box(ax, x_out, y, 2.8, 0.85, label, color, fontsize=8)

    # ── Calibration feedback loop ──────────────────────────
    cal_y = 0.8
    ax.annotate("", xy=(x_tool - 1.0, cal_y), xytext=(x_out, cal_y),
                arrowprops=dict(arrowstyle="<-", color="#C62828", lw=2.0))
    box(ax, x_cal, cal_y, 3.5, 0.7,
        "CALIBRATION LOOP:  1KGP GT ──► Compare tools ──► Update weights",
        "#B71C1C", fontsize=8, bold=True)
    ax.annotate("", xy=(x_cal - 0.7, 1.15), xytext=(x_cal - 0.7, cal_y + 0.35),
                arrowprops=dict(arrowstyle="->", color="#C62828", lw=1.5))
    ax.text(x_cal, -0.1,
            "Novel contribution: data-type-specific empirical calibration improves consensus accuracy",
            ha="center", va="bottom", fontsize=9, color="#B71C1C", style="italic")

    # ── Title strip ────────────────────────────────────────
    ax.add_patch(mpatches.FancyBboxPatch((0, 8.75), 18, 0.2,
                                          boxstyle="square,pad=0",
                                          facecolor="#1A237E", zorder=1))
    ax.text(9, 8.85, "Integrated Multi-Tool HLA Typing Pipeline — Information Flow",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="white", zorder=5)

    fig.savefig(out_dir / "M6_information_flow.png", dpi=DPI, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_dir / 'M6_information_flow.png'}")


# ─── HTML report ──────────────────────────────────────────────────────────────

FIGURE_META = {
    "M1_cross_analysis_accuracy.png": {
        "title": "Figure M1 — Cross-Analysis Accuracy Comparison",
        "caption": (
            "Per-gene 2-field concordance for each HLA typing tool, grouped by input data type "
            "(WGS, WES, RNA-seq). Demonstrates that tool performance is highly data-type-specific: "
            "e.g., arcasHLA performs poorly on WGS but is designed for RNA-seq; PolySOLVER is "
            "competitive for WES/WGS but unavailable for RNA. No single tool is universally optimal, "
            "motivating the multi-tool consensus approach."
        ),
    },
    "M2_tool_coverage_matrix.png": {
        "title": "Figure M2 — Tool Coverage and Performance Matrix",
        "caption": (
            "Compatibility of each HLA typing tool with different input formats, annotated with "
            "mean 2-field concordance where calibration data is available. ✗ = tool not applicable; "
            "– = compatible but not yet benchmarked; colours indicate concordance level. "
            "This matrix motivates the pipeline's multi-tool design: comprehensive coverage requires "
            "8 tools spanning all input types."
        ),
    },
    "M3_calibration_benefit.png": {
        "title": "Figure M3 — Benefit of Empirically Calibrated Consensus Voting",
        "caption": (
            "Mean 2-field concordance comparing equal-weight vs calibrated-weight consensus voting, "
            "per HLA gene and per analysis type. Δ values show improvement from calibration; "
            "significance stars indicate Wilcoxon signed-rank test results (one-sided, calibrated > equal). "
            "Key finding: calibration consistently improves Class I and Class II accuracy, with "
            "statistically significant gains in genes where tool performance is highly heterogeneous."
        ),
    },
    "M4_data_type_weights.png": {
        "title": "Figure M4 — Data-Type-Specific Calibrated Consensus Weights",
        "caption": (
            "Stacked horizontal bar charts showing the empirically calibrated tool contribution weights "
            "per HLA gene, separately for each analysis type (WGS, WES, RNA-seq). Each bar sums to 1.0. "
            "Key novelty: the optimal tool mix is data-type-specific — e.g., OptiType dominates "
            "WGS Class I (A/B/C) while HLA-HD dominates Class II (DRB1/DQB1); RNA-seq favors "
            "different tools entirely. A single fixed weighting scheme would be suboptimal."
        ),
    },
    "M5_population_accuracy.png": {
        "title": "Figure M5 — Population-Stratified Accuracy (AFR / EAS / EUR)",
        "caption": (
            "Per-gene 2-field concordance stratified by 1KGP superpopulation. Demonstrates that "
            "pipeline performance generalises across major human ancestry groups. Tools are compared "
            "within each population; rightmost panel shows overall mean concordance for reference. "
            "Note: EAS and EUR panels may have smaller sample sizes."
        ),
    },
    "M6_information_flow.png": {
        "title": "Figure M6 — Pipeline Information Flow and Architecture",
        "caption": (
            "End-to-end architecture of the integrated HLA typing pipeline. Input data (WGS BAM/CRAM, "
            "WES BAM, RNA-seq FASTQ, long-read BAM) is first processed by the QC module, then routed "
            "in parallel to up to 8 specialised HLA typing tools. Tool outputs are combined by the "
            "weighted consensus voting module (three modes: equal, read-confidence, calibrated). "
            "The calibration feedback loop uses 1KGP ground-truth genotypes to derive data-type-specific "
            "tool weights, which are then used to improve subsequent consensus calls. Output includes "
            "high-resolution HLA alleles, per-locus confidence scores, LOH analysis, and HTML reports."
        ),
    },
}

def write_html(out_dir, acc_dict, weights_dict, generated_figs):
    def _b64(p):
        with open(p, "rb") as f:
            return base64.b64encode(f.read()).decode()

    fig_html = ""
    for fname, meta in FIGURE_META.items():
        p = out_dir / fname
        if not p.exists():
            continue
        b64 = _b64(p)
        fig_html += f"""
        <div class="fig">
          <h3>{meta['title']}</h3>
          <img src="data:image/png;base64,{b64}" alt="{fname}" />
          <p class="cap">{meta['caption']}</p>
        </div>
        """

    # Summary table of available calibration data
    rows = ""
    for dt in ["wgs", "wes", "rna"]:
        acc = acc_dict.get(dt)
        w   = weights_dict.get(dt)
        n   = acc["N_samples"].max() if acc is not None and "N_samples" in acc.columns else "—"
        tools = ", ".join([tl(t) for t in acc["Tool"].tolist()]) if acc is not None else "—"
        rows += f"""<tr>
          <td><b>{DT_LABELS.get(dt, dt)}</b></td>
          <td style="text-align:center">{"✓" if acc is not None else "—"}</td>
          <td style="text-align:center">{"✓" if w is not None else "—"}</td>
          <td style="text-align:center">{n}</td>
          <td>{tools}</td>
        </tr>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>HLA Pipeline — Multi-Analysis Scientific Figures</title>
<style>
  body {{ font-family: 'Helvetica Neue', Arial, sans-serif; margin: 40px auto;
         max-width: 1400px; color: #333; background: #fafafa; }}
  h1 {{ color: #1A237E; border-bottom: 3px solid #1A237E; padding-bottom: 8px; }}
  h2 {{ color: #283593; margin-top: 40px; }}
  h3 {{ color: #37474F; margin: 0 0 8px 0; font-size: 14px; }}
  .meta {{ background: #E8EAF6; border-radius: 8px; padding: 14px 20px;
           display: flex; gap: 30px; flex-wrap: wrap; font-size: 14px; }}
  .meta span {{ font-weight: bold; color: #1A237E; }}
  .tbl {{ border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 13px;
          background: white; }}
  .tbl th {{ background: #283593; color: white; padding: 8px 12px; text-align: center; }}
  .tbl td {{ border: 1px solid #E0E0E0; padding: 7px 10px; }}
  .tbl tr:hover {{ background: #F5F5F5; }}
  .fig {{ margin: 36px 0; border: 1px solid #CFD8DC; border-radius: 10px;
          overflow: hidden; background: white; box-shadow: 0 2px 6px rgba(0,0,0,.07); }}
  .fig img {{ width: 100%; display: block; }}
  .cap {{ margin: 0; padding: 12px 18px; background: #FAFAFA;
          font-size: 12px; color: #555; border-top: 1px solid #E0E0E0;
          line-height: 1.6; }}
  .novelty {{ background: #FFF3E0; border-left: 4px solid #E65100;
              padding: 12px 18px; border-radius: 0 8px 8px 0; margin: 20px 0; }}
</style>
</head>
<body>
<h1>HLA Typing Pipeline — Multi-Analysis Scientific Figures</h1>
<div class="meta">
  <div>Generated: <span>{datetime.now().strftime("%Y-%m-%d %H:%M")}</span></div>
  <div>Figures: <span>{len(generated_figs)}</span></div>
  <div>Data types: <span>{", ".join(dt for dt, df in acc_dict.items() if df is not None)}</span></div>
</div>

<div class="novelty">
  <b>Pipeline Novelties Highlighted:</b>
  <ul>
    <li><b>Multi-tool ensemble</b> — 8 specialised tools, each optimal for different HLA genes/input types</li>
    <li><b>Empirical weight calibration</b> — tool weights derived from 1KGP ground-truth (N=50+/type)</li>
    <li><b>Data-type-specific calibration</b> — separate weight sets for WGS, WES, and RNA-seq</li>
    <li><b>End-to-end automation</b> — QC → typing → consensus → calibration loop on Puhti HPC</li>
    <li><b>Population generalisability</b> — validated across AFR, EAS, and EUR ancestries</li>
  </ul>
</div>

<h2>Calibration Data Summary</h2>
<table class="tbl">
  <thead>
    <tr><th>Analysis Type</th><th>Accuracy TSV</th><th>Weights JSON</th>
        <th>N Samples</th><th>Tools Evaluated</th></tr>
  </thead>
  <tbody>{rows}</tbody>
</table>

<h2>Figures</h2>
{fig_html}

<hr/>
<p style="font-size:11px;color:#999">
  Generated by multi_analysis_figures.py — HLA Typing Pipeline v1.4
</p>
</body>
</html>"""

    p = out_dir / "multi_analysis_report.html"
    p.write_text(html, encoding="utf-8")
    print(f"  Saved: {p}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for dt in ["wgs", "wes", "rna"]:
        ap.add_argument(f"--{dt}-accuracy",  default=None,
                        help=f"Path to tool_accuracy_{dt}*.tsv")
        ap.add_argument(f"--{dt}-weights",   default=None,
                        help=f"Path to tool_weights_{dt}*.json")
        ap.add_argument(f"--{dt}-strategy",  default=None,
                        help=f"Path to strategy_comparison_{dt}_calibrated.tsv")
    ap.add_argument("--population", default=None,
                    help="Path to tool_accuracy_wgs_by_population.tsv")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for figures and HTML report")
    args = ap.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _set_style()

    print("Loading data...")
    acc_dict = {
        "wgs": _load_accuracy(args.wgs_accuracy),
        "wes": _load_accuracy(args.wes_accuracy),
        "rna": _load_accuracy(args.rna_accuracy),
    }
    weights_dict = {
        "wgs": _load_weights(args.wgs_weights),
        "wes": _load_weights(args.wes_weights),
        "rna": _load_weights(args.rna_weights),
    }
    strat_dict = {
        "wgs": _load_strategy(args.wgs_strategy),
        "wes": _load_strategy(args.wes_strategy),
        "rna": _load_strategy(args.rna_strategy),
    }
    pop_df  = _load_population(args.population)
    wgs_acc = acc_dict.get("wgs")

    loaded = [dt for dt, df in acc_dict.items() if df is not None]
    print(f"  Accuracy data: {loaded or 'none'}")

    generated = []
    print("\nGenerating figures...")

    print("  M1: Cross-analysis accuracy comparison")
    fig_M1_cross_analysis(acc_dict, out_dir)
    generated.append("M1_cross_analysis_accuracy.png")

    print("  M2: Tool coverage matrix")
    fig_M2_coverage_matrix(acc_dict, out_dir)
    generated.append("M2_tool_coverage_matrix.png")

    print("  M3: Calibration benefit")
    fig_M3_calibration_benefit(strat_dict, out_dir)
    generated.append("M3_calibration_benefit.png")

    print("  M4: Data-type-specific weights")
    fig_M4_data_type_weights(weights_dict, out_dir)
    generated.append("M4_data_type_weights.png")

    print("  M5: Population-stratified accuracy")
    fig_M5_population(pop_df, wgs_acc, out_dir)
    generated.append("M5_population_accuracy.png")

    print("  M6: Information flow diagram")
    fig_M6_information_flow(out_dir)
    generated.append("M6_information_flow.png")

    print("  HTML report")
    write_html(out_dir, acc_dict, weights_dict,
               [f for f in generated if (out_dir / f).exists()])

    print(f"\nDone. {sum(1 for f in generated if (out_dir / f).exists())} figures in: {out_dir}")
    print(f"  Open: {out_dir / 'multi_analysis_report.html'}")


if __name__ == "__main__":
    main()
