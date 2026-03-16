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
    "hlahd":    "#2196F3",   # blue
    "optitype": "#4CAF50",   # green
    "spechla":  "#FF9800",   # orange
    "arcashla": "#9C27B0",   # purple
    "xhla":     "#F44336",   # red
}
TOOL_LABELS = {
    "hlahd":    "HLA-HD",
    "optitype": "OptiType",
    "spechla":  "SpecHLA",
    "arcashla": "arcasHLA",
    "xhla":     "xHLA",
}
GENE_ORDER = ["A", "B", "C", "DRB1", "DQB1"]
EXCLUDE_TOOLS = {"hlala"}   # low N=1, all-zero

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
def load_accuracy(conf_dir):
    for name in ("tool_accuracy_wgs_v2.tsv", "tool_accuracy_wgs.tsv"):
        p = conf_dir / name
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            df = df[~df["Tool"].isin(EXCLUDE_TOOLS)].copy()
            df["Tool_label"] = df["Tool"].map(tool_label)
            print(f"  Loaded accuracy: {p.name}  ({len(df)} tools)")
            return df
    raise FileNotFoundError(f"No tool_accuracy_wgs*.tsv found in {conf_dir}")

def load_weights(conf_dir):
    for name in ("tool_weights_wgs_v2.json", "tool_weights_wgs.json"):
        p = conf_dir / name
        if p.exists():
            with open(p) as f:
                d = json.load(f)
            print(f"  Loaded weights: {p.name}")
            return d
    return None

def load_population(conf_dir):
    p = conf_dir / "tool_accuracy_wgs_by_population.tsv"
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
    title = f"HLA Typing Tool Concordance Rates — 30× WGS (N≤{n} samples)"
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
# HTML report
# ---------------------------------------------------------------------------
def _b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def write_html(out_dir, acc_df, weights_data, run_title, plots_made):
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

    # Images
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
  <div>Data type: <span>WGS (30×)</span></div>
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
    ap.add_argument("--conf-dir",       required=True,
                    help="Directory containing tool_accuracy_wgs*.tsv and tool_weights_wgs*.json")
    ap.add_argument("--output-dir",     required=True,
                    help="Directory to write PNG plots and HTML report")
    ap.add_argument("--strategy-file",  default=None,
                    help="Path to strategy_comparison_calibrated.tsv (optional)")
    ap.add_argument("--title",          default=None,
                    help="Optional run title for plot suptitles")
    args = ap.parse_args()

    conf_dir = Path(args.conf_dir)
    out_dir  = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _set_style()

    print("Loading data...")
    acc_df       = load_accuracy(conf_dir)
    weights_data = load_weights(conf_dir)
    pop_df       = load_population(conf_dir)
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

    print("  HTML report")
    write_html(out_dir, acc_df, weights_data, run_title, plots_made)

    print(f"\nDone. Output in: {out_dir}")
    print(f"  Open: {out_dir / 'calibration_report.html'}")


if __name__ == "__main__":
    main()
