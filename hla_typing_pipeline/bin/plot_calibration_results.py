#!/usr/bin/env python3
"""
plot_calibration_results.py — Visualize HLA tool calibration results.

Produces 5–6 publication-quality plots and a self-contained HTML report
from the output files of calibrate_tool_weights.py.

Usage:
    python3 plot_calibration_results.py \\
        --conf-dir  /path/to/hla_typing_pipeline/conf \\
        --output-dir /path/to/plots \\
        [--strategy-file /path/to/strategy_comparison_calibrated.tsv] \\
        [--title "1KGP 30x WGS – 133 samples"]
"""

import argparse
import base64
import json
import os
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    warnings.warn("seaborn not found — heatmap will use matplotlib imshow")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TOOL_COLORS = {
    "hlahd":     "#2196F3",   # blue
    "optitype":  "#4CAF50",   # green
    "spechla":   "#FF9800",   # orange
    "arcashla":  "#9C27B0",   # purple
    "kourami":   "#F44336",   # red
    "polysolver":"#00BCD4",   # cyan
    "seq2hla":   "#FF99CC",   # pink
    "t1k":       "#795548",   # brown
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
GENE_ORDER = ["A", "B", "C", "DRB1", "DQB1"]
EXCLUDE_TOOLS = {"hlala", "xhla"}   # hlala: N=1 all-zero; xhla: removed from pipeline

STRATEGY_COLORS = {
    "equal":            "#9E9E9E",
    "read_confidence":  "#2196F3",
    "calibrated":       "#4CAF50",
}

DPI = 180
FIGSIZE_HEAT  = (10, 5)
FIGSIZE_BAR   = (13, 6)
FIGSIZE_RANK  = (10, 5)
FIGSIZE_POP   = (16, 6)
FIGSIZE_STACK = (12, 6)
FIGSIZE_STRAT = (13, 6)
FIGSIZE_BOX   = (11, 5)

# Mapping from Nextflow process name → short tool key (for resource plots)
PROCESS_TO_TOOL = {
    "HLAHD":          "hlahd",
    "HLAHD_FASTQ":    "hlahd",
    "SPECHLA":        "spechla",
    "SPECHLA_FASTQ":  "spechla",
    "ARCASHLA":       "arcashla",
    "ARCASHLA_FASTQ": "arcashla",
    "OPTITYPE":       "optitype",
    "OPTITYPE_FASTQ": "optitype",
    "XHLA":           "xhla",
    "XHLA_FASTQ":     "xhla",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def tool_color(t):
    return TOOL_COLORS.get(t, "#607D8B")

def tool_label(t):
    return TOOL_LABELS.get(t, t)

def pct(v):
    return f"{v:.0%}" if pd.notna(v) else "NA"

def _set_style():
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        try:
            plt.style.use("seaborn-whitegrid")
        except OSError:
            plt.style.use("ggplot")

def save(fig, path, title=None):
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")

def _gene_cols(df):
    """Return gene columns present in df (ordered by GENE_ORDER then extras)."""
    extras = [c for c in df.columns if c not in GENE_ORDER and c.upper() == c
              and c not in ("NA",)]
    return [g for g in GENE_ORDER if g in df.columns] + extras

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_accuracy(conf_dir, data_type="wgs", explicit_path=None):
    if explicit_path:
        candidates = [Path(explicit_path)]
    else:
        dt = data_type or "wgs"
        if dt == "wgs":
            names = ["tool_accuracy_wgs_v3.tsv", "tool_accuracy_wgs_v2.tsv", "tool_accuracy_wgs.tsv"]
        elif dt == "wes":
            names = ["tool_accuracy_wes_v1.tsv", "tool_accuracy_wes.tsv"]
        elif dt == "rna":
            names = ["tool_accuracy_rna_v1.tsv", "tool_accuracy_rna.tsv"]
        else:
            names = [f"tool_accuracy_{dt}_v1.tsv", f"tool_accuracy_{dt}.tsv"]
        candidates = [conf_dir / n for n in names]
    for p in candidates:
        if Path(p).exists():
            df = pd.read_csv(p, sep="\t")
            df = df[~df["Tool"].isin(EXCLUDE_TOOLS)].copy()
            df["Tool_label"] = df["Tool"].map(tool_label)
            print(f"  Loaded accuracy: {Path(p).name}  ({len(df)} tools)")
            return df
    raise FileNotFoundError(f"No tool_accuracy_{data_type}*.tsv found in {conf_dir}")

def load_weights(conf_dir, data_type="wgs", explicit_path=None):
    if explicit_path:
        candidates = [Path(explicit_path)]
    else:
        dt = data_type or "wgs"
        if dt == "wgs":
            names = ["tool_weights_wgs_v3.json", "tool_weights_wgs_v2.json", "tool_weights_wgs.json"]
        elif dt == "wes":
            names = ["tool_weights_wes_v1.json", "tool_weights_wes.json"]
        elif dt == "rna":
            names = ["tool_weights_rna_v1.json", "tool_weights_rna.json"]
        else:
            names = [f"tool_weights_{dt}_v1.json", f"tool_weights_{dt}.json"]
        candidates = [conf_dir / n for n in names]
    for p in candidates:
        if Path(p).exists():
            with open(p) as f:
                d = json.load(f)
            print(f"  Loaded weights: {Path(p).name}")
            return d
    return None

def load_population(conf_dir, data_type="wgs", explicit_path=None):
    if explicit_path:
        p = Path(explicit_path)
    else:
        p = conf_dir / f"tool_accuracy_{data_type}_by_population.tsv"
    if p.exists():
        df = pd.read_csv(p, sep="\t")
        df = df[~df["Tool"].isin(EXCLUDE_TOOLS)].copy()
        print(f"  Loaded population: {p.name}")
        return df
    return None

def load_strategy(path):
    if path and Path(path).exists():
        df = pd.read_csv(path, sep="\t")
        print(f"  Loaded strategy comparison: {path}")
        return df
    return None

# ---------------------------------------------------------------------------
# Plot 1 — Concordance heatmap
# ---------------------------------------------------------------------------
def plot_heatmap(acc_df, out_dir, run_title):
    genes = _gene_cols(acc_df)
    pivot = acc_df.set_index("Tool")[genes].astype(float)
    pivot.index = [tool_label(t) for t in pivot.index]

    fig, ax = plt.subplots(figsize=FIGSIZE_HEAT)
    data = pivot.values.copy()
    mask = np.isnan(data)

    if HAS_SEABORN:
        cmap = sns.color_palette("YlOrRd", as_cmap=True)
        sns.heatmap(
            pivot, ax=ax,
            vmin=0, vmax=1,
            cmap=cmap,
            annot=True, fmt=".2f",
            linewidths=0.5, linecolor="white",
            mask=pivot.isna(),
            annot_kws={"size": 11, "weight": "bold"},
            cbar_kws={"label": "2-field concordance", "shrink": 0.7},
        )
        # Grey out NA cells
        for (r, c), val in np.ndenumerate(data):
            if mask[r, c]:
                ax.add_patch(mpatches.Rectangle(
                    (c, r), 1, 1, fill=True, color="#BDBDBD", zorder=2))
                ax.text(c + 0.5, r + 0.5, "NA", ha="center", va="center",
                        fontsize=10, color="#424242", zorder=3)
    else:
        im = ax.imshow(data, vmin=0, vmax=1, cmap="YlOrRd", aspect="auto")
        plt.colorbar(im, ax=ax, label="2-field concordance", shrink=0.7)
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        for (r, c), val in np.ndenumerate(data):
            txt = f"{val:.2f}" if not mask[r, c] else "NA"
            ax.text(c, r, txt, ha="center", va="center", fontsize=10)

    ax.set_xlabel("HLA Gene", fontsize=11)
    ax.set_ylabel("Tool", fontsize=11)
    n = int(acc_df["N_samples"].max())
    dt_label = acc_df["Data_type"].iloc[0].upper() if "Data_type" in acc_df.columns else "WGS"
    title = f"HLA Typing Tool Concordance Rates — {dt_label} (N≤{n} samples)"
    if run_title:
        title = f"{run_title}\n{title}"
    save(fig, out_dir / "01_concordance_heatmap.png", title)
    return pivot

# ---------------------------------------------------------------------------
# Plot 2 — Concordance grouped bar chart
# ---------------------------------------------------------------------------
def plot_concordance_bars(acc_df, out_dir, run_title):
    genes = _gene_cols(acc_df)
    tools = acc_df["Tool"].tolist()
    n_tools = len(tools)
    n_genes = len(genes)
    x = np.arange(n_genes)
    width = 0.8 / n_tools

    fig, ax = plt.subplots(figsize=FIGSIZE_BAR)

    for i, (_, row) in enumerate(acc_df.iterrows()):
        vals = [float(row[g]) if g in row and pd.notna(row[g]) else np.nan
                for g in genes]
        offsets = x + (i - n_tools / 2 + 0.5) * width
        bars = ax.bar(offsets, [v if not np.isnan(v) else 0 for v in vals],
                      width=width * 0.9,
                      color=tool_color(row["Tool"]),
                      label=tool_label(row["Tool"]),
                      zorder=3)
        for bar, val in zip(bars, vals):
            if not np.isnan(val) and val > 0.03:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.01,
                        f"{val:.2f}", ha="center", va="bottom",
                        fontsize=7.5, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=12)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("2-Field Concordance", fontsize=11)
    ax.set_xlabel("HLA Gene", fontsize=11)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax.axhline(0.5, color="grey", lw=1, ls="--", alpha=0.5, label="_50%")
    ax.axhline(0.8, color="grey", lw=1, ls=":", alpha=0.5, label="_80%")
    ax.text(n_genes - 0.5, 0.51, "50%", fontsize=8, color="grey")
    ax.text(n_genes - 0.5, 0.81, "80%", fontsize=8, color="grey")
    ax.grid(axis="y", alpha=0.4, zorder=0)

    title = "Per-Gene 2-Field Concordance by Tool"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "02_concordance_bars.png", title)

# ---------------------------------------------------------------------------
# Plot 3 — Tool weights stacked horizontal bar
# ---------------------------------------------------------------------------
def plot_weights(weights_data, out_dir, run_title):
    if not weights_data:
        print("  [SKIP] No weights data — skipping Plot 3")
        return
    genes_dict = weights_data.get("genes", {})
    # Filter genes to plot (skip genes where all tools have equal placeholder weights)
    plot_genes = [g for g in genes_dict if g in GENE_ORDER or g in
                  ("DQA1", "DPA1", "DPB1")]
    if not plot_genes:
        print("  [SKIP] No gene weight data found")
        return

    # Collect all tool names
    all_tools = sorted({t for gd in genes_dict.values() for t in gd.keys()}
                       - EXCLUDE_TOOLS)
    fig, ax = plt.subplots(figsize=FIGSIZE_STACK)
    y = np.arange(len(plot_genes))
    lefts = np.zeros(len(plot_genes))

    for tool in all_tools:
        vals = np.array([genes_dict.get(g, {}).get(tool, 0.0)
                         for g in plot_genes])
        bars = ax.barh(y, vals, left=lefts, color=tool_color(tool),
                       label=tool_label(tool), height=0.7, edgecolor="white",
                       linewidth=0.5)
        for bar, val in zip(bars, vals):
            if val >= 0.07:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_y() + bar.get_height() / 2,
                        f"{val:.2f}", ha="center", va="center",
                        fontsize=8.5, color="white", fontweight="bold")
        lefts += vals

    ax.set_yticks(y)
    ax.set_yticklabels(plot_genes, fontsize=11)
    ax.set_xlim(0, 1.01)
    ax.set_xlabel("Cumulative Weight (sums to 1.0 per gene)", fontsize=10)
    ax.set_ylabel("HLA Gene", fontsize=11)
    ax.legend(loc="lower right", framealpha=0.9, fontsize=9,
              bbox_to_anchor=(1.0, 0.0))
    ax.axvline(1.0, color="grey", lw=1, ls="--", alpha=0.5)

    n = weights_data.get("n_samples", "?")
    title = f"Calibrated Consensus Voting Weights per HLA Gene (N={n} samples)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "03_tool_weights.png", title)

# ---------------------------------------------------------------------------
# Plot 4 — Population-stratified accuracy
# ---------------------------------------------------------------------------
def plot_population(pop_df, out_dir, run_title):
    if pop_df is None:
        print("  [SKIP] No population data — skipping Plot 4")
        return

    # Detect populations from column names (e.g. AFR_A, EAS_B, …)
    col_parts = [(c.split("_")[0], c.split("_")[1])
                 for c in pop_df.columns if "_" in c
                 and c.split("_")[1] in GENE_ORDER + ["DQA1", "DPA1", "DPB1"]]
    populations = sorted({p for p, _ in col_parts})
    if not populations:
        print("  [SKIP] No population columns detected")
        return

    genes = sorted({g for _, g in col_parts if g in GENE_ORDER},
                   key=lambda g: GENE_ORDER.index(g) if g in GENE_ORDER else 99)

    fig, axes = plt.subplots(1, len(populations), figsize=FIGSIZE_POP,
                             sharey=True)
    if len(populations) == 1:
        axes = [axes]

    pop_df = pop_df[~pop_df["Tool"].isin(EXCLUDE_TOOLS)].copy()
    tools = pop_df["Tool"].tolist()
    n_tools = len(tools)
    n_genes = len(genes)
    width = 0.8 / n_tools
    x = np.arange(n_genes)

    for ax, pop in zip(axes, populations):
        pop_cols = {g: f"{pop}_{g}" for g in genes if f"{pop}_{g}" in pop_df.columns}
        avail_genes = [g for g in genes if g in pop_cols]
        x_ = np.arange(len(avail_genes))

        for i, (_, row) in enumerate(pop_df.iterrows()):
            vals = [float(row[pop_cols[g]])
                    if pop_cols.get(g) in row and pd.notna(row.get(pop_cols.get(g, ""), np.nan))
                    else np.nan
                    for g in avail_genes]
            offs = x_ + (i - n_tools / 2 + 0.5) * width
            ax.bar(offs, [v if not np.isnan(v) else 0 for v in vals],
                   width=width * 0.9, color=tool_color(row["Tool"]),
                   label=tool_label(row["Tool"]) if pop == populations[0] else "_",
                   zorder=3)

        ax.set_xticks(x_)
        ax.set_xticklabels(avail_genes, fontsize=11)
        ax.set_ylim(0, 1.15)
        ax.set_title(pop, fontsize=13, fontweight="bold")
        ax.grid(axis="y", alpha=0.4, zorder=0)
        ax.axhline(0.5, color="grey", lw=1, ls="--", alpha=0.4)

    axes[0].set_ylabel("2-Field Concordance", fontsize=11)
    fig.legend(*axes[0].get_legend_handles_labels(),
               loc="upper right", framealpha=0.9, fontsize=9,
               bbox_to_anchor=(1.0, 1.0))

    title = "Population-Stratified Tool Accuracy (AFR / EAS / EUR)"
    if run_title:
        title = f"{run_title} — {title}"
    fig.text(0.5, -0.02,
             "Note: EAS and EUR panels may have small sample sizes (N<5)",
             ha="center", fontsize=9, color="grey")
    save(fig, out_dir / "04_population_accuracy.png", title)

# ---------------------------------------------------------------------------
# Plot 5 — Overall tool ranking
# ---------------------------------------------------------------------------
def plot_ranking(acc_df, out_dir, run_title):
    genes = _gene_cols(acc_df)
    df = acc_df.copy()
    df["Mean_num"] = pd.to_numeric(df["Mean_concordance"], errors="coerce")
    df = df.sort_values("Mean_num", ascending=True)

    fig, (ax_main, ax_detail) = plt.subplots(1, 2, figsize=FIGSIZE_RANK,
                                              gridspec_kw={"width_ratios": [2, 3]})

    # Left — horizontal ranking bar
    colors = [tool_color(t) for t in df["Tool"]]
    bars = ax_main.barh(
        [tool_label(t) for t in df["Tool"]],
        df["Mean_num"].fillna(0),
        color=colors, height=0.6, edgecolor="white")
    for bar, (_, row) in zip(bars, df.iterrows()):
        v = row["Mean_num"]
        n = int(row["N_samples"]) if pd.notna(row.get("N_samples")) else "?"
        ax_main.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                     f"{v:.3f}  (N={n})" if pd.notna(v) else "NA",
                     va="center", fontsize=9)
    ax_main.set_xlim(0, 1.2)
    ax_main.set_xlabel("Mean 2-Field Concordance", fontsize=10)
    ax_main.set_title("Overall Ranking", fontsize=11)
    ax_main.grid(axis="x", alpha=0.4)

    # Right — per-gene detail
    n_genes = len(genes)
    n_tools = len(df)
    width = 0.8 / n_genes
    x = np.arange(n_tools)

    gene_palette = plt.cm.tab10(np.linspace(0, 1, n_genes))
    for gi, (gene, color) in enumerate(zip(genes, gene_palette)):
        if gene not in df.columns:
            continue
        vals = pd.to_numeric(df[gene], errors="coerce").fillna(0).values
        offs = x + (gi - n_genes / 2 + 0.5) * width
        ax_detail.bar(offs, vals, width=width * 0.9, color=color,
                      label=gene, zorder=3)

    ax_detail.set_xticks(x)
    ax_detail.set_xticklabels([tool_label(t) for t in df["Tool"]],
                               fontsize=9, rotation=20, ha="right")
    ax_detail.set_ylim(0, 1.1)
    ax_detail.set_ylabel("2-Field Concordance", fontsize=10)
    ax_detail.set_title("Per-Gene Breakdown", fontsize=11)
    ax_detail.legend(title="Gene", loc="upper left", fontsize=8,
                     framealpha=0.9, ncol=2)
    ax_detail.grid(axis="y", alpha=0.4, zorder=0)

    title = "HLA Tool Performance Ranking (Mean 2-Field Concordance)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "05_tool_ranking.png", title)

# ---------------------------------------------------------------------------
# Plot 6 — Voting strategy comparison (optional)
# ---------------------------------------------------------------------------
def plot_strategy(strat_df, out_dir, run_title):
    if strat_df is None:
        print("  [SKIP] No strategy comparison data — skipping Plot 6")
        return

    # Detect mode columns
    mode_cols = [c for c in strat_df.columns
                 if c.endswith("_mean") and c in
                 ("equal_mean", "read_confidence_mean", "calibrated_mean")]
    if not mode_cols:
        # Fallback: any *_mean column
        mode_cols = [c for c in strat_df.columns if c.endswith("_mean")]
    if not mode_cols:
        print("  [SKIP] Strategy file has no *_mean columns")
        return

    gene_col = "Gene" if "Gene" in strat_df.columns else strat_df.columns[0]
    sig_col  = "Significant_p05" if "Significant_p05" in strat_df.columns else None
    p_col    = "Wilcoxon_p" if "Wilcoxon_p" in strat_df.columns else None

    genes = strat_df[gene_col].tolist()
    n_genes = len(genes)
    n_modes = len(mode_cols)
    x = np.arange(n_genes)
    width = 0.8 / n_modes

    fig, ax = plt.subplots(figsize=FIGSIZE_STRAT)

    for mi, col in enumerate(mode_cols):
        mode_name = col.replace("_mean", "")
        vals = pd.to_numeric(strat_df[col], errors="coerce").fillna(0).values
        offs = x + (mi - n_modes / 2 + 0.5) * width
        color = STRATEGY_COLORS.get(mode_name, "#607D8B")
        bars = ax.bar(offs, vals, width=width * 0.9, color=color,
                      label=mode_name.replace("_", " ").title(), zorder=3,
                      edgecolor="white")
        for bar, val in zip(bars, vals):
            if val > 0.02:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.01,
                        f"{val:.2f}", ha="center", va="bottom", fontsize=8,
                        rotation=45)

    # Significance markers on calibrated bars (last mode)
    if sig_col and p_col and "calibrated_mean" in mode_cols:
        ci = mode_cols.index("calibrated_mean")
        for gi, (_, row) in enumerate(strat_df.iterrows()):
            sig = str(row.get(sig_col, "no")).strip().lower()
            if sig in ("yes", "true", "1"):
                offs_cal = gi + (ci - n_modes / 2 + 0.5) * width
                best_val = float(row.get("calibrated_mean", 0) or 0)
                ax.text(offs_cal, best_val + 0.07, "*", ha="center",
                        fontsize=14, color="black", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=11)
    ax.set_ylim(0, 1.2)
    ax.set_ylabel("Mean 2-Field Concordance", fontsize=11)
    ax.set_xlabel("HLA Gene", fontsize=11)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.9)
    ax.grid(axis="y", alpha=0.4, zorder=0)
    ax.text(0, -0.13, "* p < 0.05 (Wilcoxon signed-rank test, calibrated vs equal)",
            fontsize=8, color="grey", transform=ax.transAxes)

    title = "Voting Strategy Comparison (Equal / Read-Confidence / Calibrated)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "06_strategy_comparison.png", title)

# ---------------------------------------------------------------------------
# Trace file parsing (resource usage)
# ---------------------------------------------------------------------------

def _parse_duration(s):
    """Convert Nextflow duration string to seconds. e.g. '8m 19s' → 499.0"""
    if not s or s.strip() in ("-", ""):
        return 0.0
    total = 0.0
    for val, unit in re.findall(r"([\d.]+)\s*(d|h|m|s|ms|us)", s):
        v = float(val)
        if unit == "d":   total += v * 86400
        elif unit == "h": total += v * 3600
        elif unit == "m": total += v * 60
        elif unit == "s": total += v
        elif unit == "ms": total += v / 1000
        elif unit == "us": total += v / 1_000_000
    return total


def _parse_memory(s):
    """Convert Nextflow memory string to GB. e.g. '2.6 GB' → 2.6"""
    if not s or s.strip() in ("-", ""):
        return 0.0
    m = re.match(r"([\d.]+)\s*(B|KB|MB|GB|TB)", s.strip(), re.IGNORECASE)
    if not m:
        return 0.0
    val, unit = float(m.group(1)), m.group(2).upper()
    return {"B": val/1e9, "KB": val/1e6, "MB": val/1e3, "GB": val, "TB": val*1e3}.get(unit, 0.0)


def _parse_cpu(s):
    """Parse CPU percentage. e.g. '659.0%' → 659.0"""
    if not s or s.strip() in ("-", ""):
        return 0.0
    return float(s.strip().rstrip("%"))


def load_traces(trace_dir):
    """
    Load all trace_*.txt files from *trace_dir*, keep only COMPLETED
    HLA typing process rows.  Returns list of dicts with keys:
        tool, realtime_s, peak_rss_gb, cpu_pct, sample
    """
    trace_dir = Path(trace_dir)
    trace_files = sorted(trace_dir.glob("trace_*.txt"))
    if not trace_files:
        # Also try plain trace.txt
        trace_files = sorted(trace_dir.glob("trace*.txt"))
    if not trace_files:
        print(f"  [WARN] No trace_*.txt files found in {trace_dir}")
        return []

    records = []
    for tf in trace_files:
        try:
            with open(tf) as fh:
                header = fh.readline().strip()
                sep = "\t" if "\t" in header else None
                cols = header.split("\t") if sep else re.split(r"\s{2,}", header)
                cols = [c.strip() for c in cols]
                idx = {c: i for i, c in enumerate(cols)}

                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t") if sep else re.split(r"\s{2,}", line)
                    parts = [p.strip() for p in parts]

                    def _get(col, default=""):
                        i = idx.get(col, -1)
                        return parts[i] if 0 <= i < len(parts) else default

                    status = _get("status")
                    if status.lower() not in ("completed", "ok", "cached"):
                        continue

                    name_field = _get("name", "")
                    process = name_field.split(" ")[0].split(":")[0].upper()
                    tool = PROCESS_TO_TOOL.get(process)
                    if tool is None:
                        continue

                    # Sample name from parentheses: "HLAHD_FASTQ (NA18526)"
                    sm = re.search(r"\(([^)]+)\)", name_field)
                    sample = sm.group(1) if sm else ""

                    realtime_s = _parse_duration(_get("realtime"))
                    peak_rss_gb = _parse_memory(_get("peak_rss"))
                    cpu_pct = _parse_cpu(_get("%cpu"))

                    if realtime_s > 0:
                        records.append({
                            "tool":        tool,
                            "sample":      sample,
                            "realtime_s":  realtime_s,
                            "peak_rss_gb": peak_rss_gb,
                            "cpu_pct":     cpu_pct,
                        })
        except Exception as e:
            print(f"  [WARN] Could not parse {tf.name}: {e}")

    print(f"  Loaded {len(records)} completed typing task records from {len(trace_files)} trace file(s)")
    return records


def _boxplot_h(ax, data_by_tool, metric_fn, tools_sorted, color_fn):
    """Draw horizontal boxplot with overlaid jittered points."""
    vals = [metric_fn(data_by_tool[t]) for t in tools_sorted]
    rng = np.random.default_rng(42)

    bp = ax.boxplot(
        vals,
        vert=False,
        patch_artist=True,
        widths=0.5,
        flierprops=dict(marker="", linestyle="none"),
        medianprops=dict(color="white", linewidth=2),
        whiskerprops=dict(color="#555", linewidth=1),
        capprops=dict(color="#555", linewidth=1),
        boxprops=dict(linewidth=0.5),
    )
    for patch, tool in zip(bp["boxes"], tools_sorted):
        patch.set_facecolor(color_fn(tool))
        patch.set_alpha(0.75)

    # Jitter overlay
    for i, (tool, v) in enumerate(zip(tools_sorted, vals), start=1):
        if not v:
            continue
        jitter = rng.normal(0, 0.08, len(v))
        ax.scatter(v, [i + j for j in jitter],
                   color=color_fn(tool), alpha=0.6, s=18, zorder=5,
                   edgecolors="white", linewidths=0.3)
        # Median annotation
        med = np.median(v)
        ax.text(med, i + 0.38, f"{med:.1f}", ha="center", va="bottom",
                fontsize=7.5, color="#222")

    ax.set_yticks(range(1, len(tools_sorted) + 1))
    ax.set_yticklabels([tool_label(t) for t in tools_sorted], fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")


# ---------------------------------------------------------------------------
# Plot 7 — Tool runtime distribution
# ---------------------------------------------------------------------------
def plot_resource_runtime(records, out_dir, run_title):
    if not records:
        print("  [SKIP] No trace data — skipping Plot 7")
        return

    from collections import defaultdict
    by_tool = defaultdict(list)
    for r in records:
        if r["realtime_s"] > 0:
            by_tool[r["tool"]].append(r["realtime_s"] / 60.0)  # → minutes

    tools = sorted(by_tool, key=lambda t: np.median(by_tool[t]), reverse=True)
    if not tools:
        print("  [SKIP] No runtime data")
        return

    n_total = sum(len(by_tool[t]) for t in tools)
    fig, ax = plt.subplots(figsize=FIGSIZE_BOX)
    _boxplot_h(ax, by_tool, lambda v: v, tools, tool_color)

    # N label per tool
    for i, t in enumerate(tools, start=1):
        n = len(by_tool[t])
        ax.text(ax.get_xlim()[1] * 0.97, i, f"n={n}",
                va="center", ha="right", fontsize=8, color="#555")

    # Log scale if range > 10×
    all_vals = [v for t in tools for v in by_tool[t]]
    if max(all_vals) / max(min(all_vals), 0.01) > 10:
        ax.set_xscale("log")
        ax.set_xlabel("Runtime (minutes, log scale)", fontsize=10)
    else:
        ax.set_xlabel("Runtime (minutes)", fontsize=10)

    title = f"HLA Typing Tool Runtime Distribution (N={n_total} tasks)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "07_tool_runtime.png", title)


# ---------------------------------------------------------------------------
# Plot 8 — Peak RAM distribution
# ---------------------------------------------------------------------------
def plot_resource_ram(records, out_dir, run_title, system_ram_gb=32.0):
    if not records:
        print("  [SKIP] No trace data — skipping Plot 8")
        return

    from collections import defaultdict
    by_tool = defaultdict(list)
    for r in records:
        if r["peak_rss_gb"] > 0:
            by_tool[r["tool"]].append(r["peak_rss_gb"])

    tools = sorted(by_tool, key=lambda t: np.median(by_tool[t]), reverse=True)
    if not tools:
        print("  [SKIP] No RAM data")
        return

    fig, ax = plt.subplots(figsize=FIGSIZE_BOX)
    _boxplot_h(ax, by_tool, lambda v: v, tools, tool_color)

    if system_ram_gb and system_ram_gb > 0:
        ax.axvline(system_ram_gb, color="#E63946", linestyle="--", linewidth=1.5,
                   alpha=0.7, label=f"Node RAM ({system_ram_gb:.0f} GB)")
        ax.legend(fontsize=9, loc="lower right")

    ax.set_xlabel("Peak RSS (GB)", fontsize=10)

    n_total = sum(len(by_tool[t]) for t in tools)
    title = f"Peak RAM Usage per HLA Typing Tool (N={n_total} tasks)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "08_tool_ram.png", title)


# ---------------------------------------------------------------------------
# Plot 9 — CPU utilisation distribution
# ---------------------------------------------------------------------------
def plot_resource_cpu(records, out_dir, run_title):
    if not records:
        print("  [SKIP] No trace data — skipping Plot 9")
        return

    from collections import defaultdict
    by_tool = defaultdict(list)
    for r in records:
        if r["cpu_pct"] > 0:
            by_tool[r["tool"]].append(r["cpu_pct"])

    tools = sorted(by_tool, key=lambda t: np.median(by_tool[t]), reverse=True)
    if not tools:
        print("  [SKIP] No CPU data")
        return

    fig, ax = plt.subplots(figsize=FIGSIZE_BOX)
    _boxplot_h(ax, by_tool, lambda v: v, tools, tool_color)

    # Reference lines at common core multiples
    xmax = max(v for t in tools for v in by_tool[t]) * 1.1
    for pct, label in [(100, "1 core"), (200, "2"), (400, "4"), (800, "8")]:
        if pct < xmax:
            ax.axvline(pct, color="grey", linestyle=":", linewidth=0.8, alpha=0.6)
            ax.text(pct + 2, len(tools) + 0.5, label, fontsize=7.5, color="grey", va="top")

    ax.set_xlabel("CPU Usage (%)", fontsize=10)

    n_total = sum(len(by_tool[t]) for t in tools)
    title = f"CPU Utilisation per HLA Typing Tool (N={n_total} tasks)"
    if run_title:
        title = f"{run_title} — {title}"
    save(fig, out_dir / "09_tool_cpu.png", title)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
def _b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def write_html(out_dir, acc_df, weights_data, run_title, plots_made, trace_records=None, data_type="wgs"):
    captions = {
        "01_concordance_heatmap.png":
            "Tool concordance rates across HLA genes. Color intensity indicates 2-field "
            "concordance (0–1); grey cells indicate genes not typed by that tool.",
        "02_concordance_bars.png":
            "Per-gene grouped bar chart. Each group shows all tools for that gene. "
            "Dashed lines at 50% and 80% concordance.",
        "03_tool_weights.png":
            "Calibrated consensus voting weights. Weights per gene sum to 1.0. Tools with "
            "higher concordance receive proportionally higher weight in consensus calls.",
        "04_population_accuracy.png":
            "Population-stratified concordance by superpopulation. EAS and EUR panels "
            "may have small N and should be interpreted with caution.",
        "05_tool_ranking.png":
            "Left: tools ranked by mean concordance. Right: per-gene breakdown. "
            "Mean is computed over genes typed by that tool.",
        "06_strategy_comparison.png":
            "Comparison of three voting strategies. Asterisk (*) marks genes where the "
            "calibrated strategy is significantly better (Wilcoxon p < 0.05).",
        "07_tool_runtime.png":
            "Wall-clock runtime per HLA typing tool from Nextflow trace files. "
            "Box shows IQR; whiskers extend to 1.5×IQR; circles are individual sample runs.",
        "08_tool_ram.png":
            "Peak RSS memory usage per HLA typing tool. "
            "Red dashed line indicates available node RAM.",
        "09_tool_cpu.png":
            "CPU utilisation per tool (%). Values >100% indicate multi-core parallelism. "
            "Grey dotted lines mark 1-, 2-, 4-, 8-core reference levels.",
    }

    genes = _gene_cols(acc_df)
    n_samples = int(acc_df["N_samples"].max()) if "N_samples" in acc_df.columns else "?"

    # Accuracy table HTML
    tbl_rows = []
    for _, row in acc_df.sort_values("Mean_concordance",
                                      ascending=False, na_position="last").iterrows():
        cells = f"<td><b>{tool_label(row['Tool'])}</b></td>"
        cells += f"<td>{row.get('N_samples', '?')}</td>"
        for g in genes:
            v = row.get(g, float("nan"))
            try:
                fv = float(v)
                color = ("#c8e6c9" if fv >= 0.8 else
                         "#fff9c4" if fv >= 0.5 else
                         "#ffcdd2" if fv > 0 else "#f5f5f5")
                cells += f'<td style="background:{color};text-align:center">{fv:.2f}</td>'
            except (ValueError, TypeError):
                cells += '<td style="text-align:center;color:#999">NA</td>'
        mc = row.get("Mean_concordance", float("nan"))
        try:
            fmc = float(mc)
            cells += f'<td style="text-align:center"><b>{fmc:.3f}</b></td>'
        except (ValueError, TypeError):
            cells += '<td style="text-align:center">NA</td>'
        tbl_rows.append(f"<tr>{cells}</tr>")

    gene_headers = "".join(f"<th>{g}</th>" for g in genes)
    table_html = f"""
    <table class="tbl">
      <thead>
        <tr>
          <th>Tool</th><th>N</th>{gene_headers}<th>Mean</th>
        </tr>
      </thead>
      <tbody>{''.join(tbl_rows)}</tbody>
    </table>
    """

    # Weights table
    weights_html = ""
    if weights_data:
        gd = weights_data.get("genes", {})
        w_genes = list(gd.keys())
        w_tools = sorted({t for vals in gd.values() for t in vals if t not in EXCLUDE_TOOLS})
        hdr = "".join(f"<th>{g}</th>" for g in w_genes)
        w_rows = []
        for t in w_tools:
            cells = f"<td><b>{tool_label(t)}</b></td>"
            for g in w_genes:
                v = gd.get(g, {}).get(t, 0.0)
                color = ("#c8e6c9" if v >= 0.4 else
                         "#fff9c4" if v >= 0.2 else "#f5f5f5")
                cells += f'<td style="background:{color};text-align:center">{v:.3f}</td>'
            w_rows.append(f"<tr>{cells}</tr>")
        weights_html = f"""
        <h2>Calibrated Tool Weights</h2>
        <p>Weights per gene sum to 1.0. Used for weighted consensus voting.</p>
        <table class="tbl">
          <thead><tr><th>Tool</th>{hdr}</tr></thead>
          <tbody>{''.join(w_rows)}</tbody>
        </table>
        """

    # Accuracy / weights images (plots 01-06)
    imgs_html = ""
    for fname in [
        "01_concordance_heatmap.png", "02_concordance_bars.png",
        "03_tool_weights.png", "04_population_accuracy.png",
        "05_tool_ranking.png", "06_strategy_comparison.png",
    ]:
        p = out_dir / fname
        if not p.exists():
            continue
        cap = captions.get(fname, "")
        b64 = _b64(p)
        imgs_html += f"""
        <div class="fig">
          <img src="data:image/png;base64,{b64}" alt="{fname}" />
          <p class="cap">{cap}</p>
        </div>
        """

    # Resource usage images (plots 07-09) + summary table
    resource_html = ""
    resource_fnames = ["07_tool_runtime.png", "08_tool_ram.png", "09_tool_cpu.png"]
    resource_plots = [f for f in resource_fnames if (out_dir / f).exists()]
    if resource_plots:
        res_imgs = ""
        for fname in resource_plots:
            cap = captions.get(fname, "")
            b64 = _b64(out_dir / fname)
            res_imgs += f"""
            <div class="fig">
              <img src="data:image/png;base64,{b64}" alt="{fname}" />
              <p class="cap">{cap}</p>
            </div>
            """

        # Summary table from trace_records
        res_tbl = ""
        if trace_records:
            from collections import defaultdict
            by_tool = defaultdict(list)
            for r in trace_records:
                by_tool[r["tool"]].append(r)
            tbl_rows2 = []
            for tool in sorted(by_tool):
                recs = by_tool[tool]
                runtimes = [r["realtime_s"] / 60 for r in recs if r["realtime_s"] > 0]
                rams = [r["peak_rss_gb"] for r in recs if r["peak_rss_gb"] > 0]
                cpus = [r["cpu_pct"] for r in recs if r["cpu_pct"] > 0]
                row_cells = (
                    f"<td><b>{tool_label(tool)}</b></td>"
                    f"<td style='text-align:center'>{len(recs)}</td>"
                    f"<td style='text-align:center'>{np.median(runtimes):.1f} min</td>"
                    f"<td style='text-align:center'>{np.median(rams):.2f} GB</td>"
                    f"<td style='text-align:center'>{np.median(cpus):.0f}%</td>"
                )
                tbl_rows2.append(f"<tr>{row_cells}</tr>")
            res_tbl = f"""
            <h3>Resource Usage Summary (median per tool)</h3>
            <table class="tbl">
              <thead><tr>
                <th>Tool</th><th>N tasks</th>
                <th>Median Runtime</th><th>Median Peak RAM</th><th>Median CPU</th>
              </tr></thead>
              <tbody>{''.join(tbl_rows2)}</tbody>
            </table>
            """

        resource_html = f"""
        <h2>Resource Usage</h2>
        <p>CPU, RAM, and runtime metrics extracted from Nextflow trace files
           (COMPLETED tasks only).</p>
        {res_tbl}
        {res_imgs}
        """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>HLA Calibration Report</title>
<style>
  body {{ font-family: 'Helvetica Neue', Arial, sans-serif; margin: 40px auto;
         max-width: 1400px; color: #333; }}
  h1 {{ color: #1a237e; border-bottom: 3px solid #1a237e; padding-bottom: 8px; }}
  h2 {{ color: #283593; margin-top: 40px; }}
  .meta {{ background: #e8eaf6; border-radius: 8px; padding: 14px 20px;
           display: flex; gap: 30px; flex-wrap: wrap; font-size: 14px; }}
  .meta span {{ font-weight: bold; color: #1a237e; }}
  .tbl {{ border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 13px; }}
  .tbl th {{ background: #283593; color: white; padding: 8px 10px;
             text-align: center; }}
  .tbl td {{ border: 1px solid #e0e0e0; padding: 6px 10px; }}
  .tbl tr:hover {{ background: #f5f5f5; }}
  .fig {{ margin: 30px 0; border: 1px solid #e0e0e0; border-radius: 8px;
          overflow: hidden; }}
  .fig img {{ width: 100%; display: block; }}
  .cap {{ margin: 0; padding: 10px 16px; background: #fafafa;
          font-size: 12px; color: #555; border-top: 1px solid #e0e0e0; }}
  .legend {{ display: flex; gap: 16px; flex-wrap: wrap; margin: 10px 0; }}
  .dot {{ width: 14px; height: 14px; border-radius: 3px; display: inline-block;
          vertical-align: middle; margin-right: 4px; }}
</style>
</head>
<body>
<h1>HLA Tool Calibration Report</h1>
<div class="meta">
  <div>Title: <span>{run_title or "HLA Calibration"}</span></div>
  <div>Generated: <span>{datetime.now().strftime("%Y-%m-%d %H:%M")}</span></div>
  <div>Samples: <span>≤{n_samples}</span></div>
  <div>Resolution: <span>2-field</span></div>
  <div>Data type: <span>{data_type.upper()}</span></div>
</div>

<h2>Concordance Summary</h2>
<p>2-field concordance per tool per gene compared to 1KGP Phase 1 ground truth.
   Color coding: <span style="background:#c8e6c9;padding:2px 6px">≥0.8</span>
   <span style="background:#fff9c4;padding:2px 6px">0.5–0.8</span>
   <span style="background:#ffcdd2;padding:2px 6px">&lt;0.5</span></p>
{table_html}

{weights_html}

<h2>Figures</h2>
{imgs_html}

{resource_html}

<hr/>
<p style="font-size:11px;color:#999">
  Generated by plot_calibration_results.py — HLA Typing Pipeline v1.4
</p>
</body>
</html>"""

    out = out_dir / "calibration_report.html"
    out.write_text(html, encoding="utf-8")
    print(f"  Saved: {out}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--conf-dir",        required=True,
                    help="Directory containing tool_accuracy_*.tsv and tool_weights_*.json")
    ap.add_argument("--output-dir",      required=True,
                    help="Directory to write PNG plots and HTML report")
    ap.add_argument("--data-type",       choices=["wgs", "wes", "rna"], default=None,
                    help="Analysis type (auto-detected from accuracy filename if omitted)")
    ap.add_argument("--accuracy-file",   default=None,
                    help="Explicit path to tool_accuracy_*.tsv (overrides --data-type lookup)")
    ap.add_argument("--weights-file",    default=None,
                    help="Explicit path to tool_weights_*.json (overrides --data-type lookup)")
    ap.add_argument("--population-file", default=None,
                    help="Explicit path to tool_accuracy_*_by_population.tsv")
    ap.add_argument("--strategy-file",   default=None,
                    help="Path to strategy_comparison_calibrated.tsv (optional)")
    ap.add_argument("--title",           default=None,
                    help="Optional run title for plot suptitles")
    ap.add_argument("--trace-dir",       default=None,
                    help="Directory containing Nextflow trace_*.txt files for resource plots")
    ap.add_argument("--system-ram",      type=float, default=32.0,
                    help="Available node RAM in GB for reference line in RAM plot (default: 32)")
    args = ap.parse_args()

    conf_dir = Path(args.conf_dir)
    out_dir  = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Auto-detect data type from explicit accuracy filename if not specified
    data_type = args.data_type
    if data_type is None and args.accuracy_file:
        fname = Path(args.accuracy_file).name
        if "rna" in fname:   data_type = "rna"
        elif "wes" in fname: data_type = "wes"
        else:                data_type = "wgs"
    data_type = data_type or "wgs"

    _set_style()

    print("Loading data...")
    acc_df       = load_accuracy(conf_dir, data_type, args.accuracy_file)
    weights_data = load_weights(conf_dir, data_type, args.weights_file)
    pop_df       = load_population(conf_dir, data_type, args.population_file)
    strat_df     = load_strategy(args.strategy_file)

    run_title = args.title or ""
    plots_made = []

    print("Generating plots...")

    print("  Plot 1: Concordance heatmap")
    plot_heatmap(acc_df, out_dir, run_title)
    plots_made.append("01_concordance_heatmap.png")

    print("  Plot 2: Concordance bar chart")
    plot_concordance_bars(acc_df, out_dir, run_title)
    plots_made.append("02_concordance_bars.png")

    print("  Plot 3: Tool weights")
    plot_weights(weights_data, out_dir, run_title)
    plots_made.append("03_tool_weights.png")

    print("  Plot 4: Population-stratified accuracy")
    plot_population(pop_df, out_dir, run_title)
    plots_made.append("04_population_accuracy.png")

    print("  Plot 5: Tool ranking")
    plot_ranking(acc_df, out_dir, run_title)
    plots_made.append("05_tool_ranking.png")

    print("  Plot 6: Voting strategy comparison")
    plot_strategy(strat_df, out_dir, run_title)
    plots_made.append("06_strategy_comparison.png")

    # Resource usage plots (optional)
    trace_records = []
    if args.trace_dir:
        print("Loading trace files...")
        trace_records = load_traces(args.trace_dir)
        if trace_records:
            print("  Plot 7: Tool runtime distribution")
            plot_resource_runtime(trace_records, out_dir, run_title)
            plots_made.append("07_tool_runtime.png")

            print("  Plot 8: Peak RAM distribution")
            plot_resource_ram(trace_records, out_dir, run_title, args.system_ram)
            plots_made.append("08_tool_ram.png")

            print("  Plot 9: CPU utilisation distribution")
            plot_resource_cpu(trace_records, out_dir, run_title)
            plots_made.append("09_tool_cpu.png")

    print("  HTML report")
    write_html(out_dir, acc_df, weights_data, run_title, plots_made, trace_records, data_type)

    print(f"\nDone. Output in: {out_dir}")
    print(f"  Open: {out_dir / 'calibration_report.html'}")


if __name__ == "__main__":
    main()
