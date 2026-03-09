#!/usr/bin/env python3
"""
HLA Typing Pipeline Execution Metrics Visualizer

Parses Nextflow trace files and generates publication-quality plots showing
CPU usage, RAM usage, runtime, parallelism (Gantt), and resource efficiency.

Usage:
    hla_pipeline_metrics.py --trace results/pipeline_info/trace_*.txt --outdir plots/
    hla_pipeline_metrics.py --trace trace.txt --outdir . --system-ram 14

Output:
    metrics_gantt.png          — Execution timeline (Gantt chart)
    metrics_cpu.png            — CPU usage per tool (violin plot)
    metrics_ram.png            — Peak RAM per tool (bar chart)
    metrics_runtime.png        — Runtime distribution (horizontal bar)
    metrics_efficiency.png     — Resource efficiency scatter
    execution_metrics_report.html  — Interactive HTML report (Plotly)
"""

import argparse
import json
import logging
import os
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Tool groups for colour coding
# ---------------------------------------------------------------------------
TOOL_GROUPS = {
    # Extraction / QC
    "CHECK_BAM":              "qc",
    "EXTRACT_HLA_READS":      "qc",
    "EXTRACT_HLA_READS_HLAHD": "qc",
    "BAM_TO_FASTQ":           "qc",
    "QC_BAM":                 "qc",
    "QC_FASTQ":               "qc",
    "FASTQC_BAM":             "qc",
    "FASTQC_FASTQ":           "qc",
    "BAMQC":                  "qc",
    "BAMQC_FASTQ":            "qc",
    "MULTIQC":                "qc",
    # HLA typing tools
    "SPECHLA":                "typing",
    "SPECHLA_FASTQ":          "typing",
    "HLAHD":                  "typing",
    "HLAHD_FASTQ":            "typing",
    "ARCASHLA":               "typing",
    "ARCASHLA_FASTQ":         "typing",
    "OPTITYPE":               "typing",
    "OPTITYPE_FASTQ":         "typing",
    "XHLA":                   "typing",
    "XHLA_FASTQ":             "typing",
    "HLALA":                  "typing",
    "HLASCAN":                "typing",
    "HLASCAN_FASTQ":          "typing",
    # Post-processing
    "CONSENSUS":              "postprocess",
    "MULTISOURCE_CONSENSUS":  "postprocess",
    "HLA_VISUALIZE":          "postprocess",
    "HLA_SUMMARY_REPORT":     "postprocess",
    "HLA_LOH":                "postprocess",
    "HLA_LOH_VISUALIZE":      "postprocess",
    "HLA_LOH_SUMMARY":        "postprocess",
    "HLA_PIPELINE_METRICS":   "postprocess",
}

GROUP_COLORS = {
    "qc":          "#4ECDC4",   # teal
    "typing":      "#FF6B6B",   # coral-red
    "postprocess": "#95E1D3",   # mint
    "unknown":     "#C9C9C9",   # grey
}

GROUP_LABELS = {
    "qc":          "QC / Extraction",
    "typing":      "HLA Typing",
    "postprocess": "Post-processing",
    "unknown":     "Other",
}


# ---------------------------------------------------------------------------
# Trace file parser
# ---------------------------------------------------------------------------

def _parse_duration(s: str) -> float:
    """Convert Nextflow duration string to seconds. e.g. '1m 23s' → 83.0"""
    if not s or s in ("-", ""):
        return 0.0
    s = s.strip()
    total = 0.0
    # Handle formats: "1h 23m 4s", "23m", "5s", "1d 2h", "100ms"
    for val, unit in re.findall(r"([\d.]+)\s*(d|h|m|s|ms|us)", s):
        v = float(val)
        if unit == "d":
            total += v * 86400
        elif unit == "h":
            total += v * 3600
        elif unit == "m":
            total += v * 60
        elif unit == "s":
            total += v
        elif unit == "ms":
            total += v / 1000
        elif unit == "us":
            total += v / 1_000_000
    return total


def _parse_memory(s: str) -> float:
    """Convert Nextflow memory string to GB. e.g. '1.2 GB' → 1.2, '512 MB' → 0.5"""
    if not s or s in ("-", ""):
        return 0.0
    s = s.strip()
    m = re.match(r"([\d.]+)\s*(KB|MB|GB|TB|B)", s, re.IGNORECASE)
    if not m:
        return 0.0
    val = float(m.group(1))
    unit = m.group(2).upper()
    if unit == "B":
        return val / 1e9
    if unit == "KB":
        return val / 1e6
    if unit == "MB":
        return val / 1e3
    if unit == "GB":
        return val
    if unit == "TB":
        return val * 1e3
    return 0.0


def _parse_cpu(s: str) -> float:
    """Parse CPU percentage. e.g. '150%' → 150.0, '1.5' → 1.5"""
    if not s or s in ("-", ""):
        return 0.0
    return float(s.strip().rstrip("%"))


def _tool_name(name_field: str) -> str:
    """Extract base tool name from Nextflow task name. e.g. 'HLAHD (NA12878)' → 'HLAHD'"""
    return name_field.strip().split(" ")[0].split(":")[0]


def load_trace(trace_path: str) -> List[dict]:
    """
    Parse a Nextflow trace file into a list of task records.

    Returns list of dicts with keys:
        task_id, name, tool, status, duration_s, realtime_s,
        cpu_pct, peak_rss_gb, peak_vmem_gb, rchar_mb, wchar_mb,
        submit_epoch, start_epoch
    """
    records = []
    with open(trace_path, "r") as f:
        header_line = f.readline().strip()
        # Handle tab or multi-space separator
        if "\t" in header_line:
            cols = header_line.split("\t")
        else:
            cols = re.split(r"\s{2,}", header_line)
        cols = [c.strip() for c in cols]

        # Build index map
        idx = {c: i for i, c in enumerate(cols)}
        required = {"name", "status", "duration", "realtime", "%cpu", "peak_rss", "peak_vmem"}
        missing = required - set(idx.keys())
        if missing:
            logger.warning(f"Trace missing columns: {missing}. Available: {cols}")

        for line in f:
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                parts = line.split("\t")
            else:
                parts = re.split(r"\s{2,}", line)
            parts = [p.strip() for p in parts]

            if len(parts) < max(idx.values(), default=0) + 1:
                # Pad short rows
                parts += [""] * (max(idx.values(), default=0) + 1 - len(parts))

            def get(col, default=""):
                i = idx.get(col, -1)
                if i < 0 or i >= len(parts):
                    return default
                return parts[i]

            name = get("name", "unknown")
            tool = _tool_name(name)
            status = get("status", "unknown")
            duration_s = _parse_duration(get("duration"))
            realtime_s = _parse_duration(get("realtime"))
            cpu_pct = _parse_cpu(get("%cpu"))
            peak_rss_gb = _parse_memory(get("peak_rss"))
            peak_vmem_gb = _parse_memory(get("peak_vmem"))
            rchar_mb = _parse_memory(get("rchar", "0 MB"))
            wchar_mb = _parse_memory(get("wchar", "0 MB"))

            # Parse timestamps for Gantt
            submit_str = get("submit", "")
            start_str = get("start", "")
            complete_str = get("complete", "")

            def to_ms(ts):
                # Nextflow timestamps: epoch ms or ISO-8601
                if not ts or ts == "-":
                    return None
                try:
                    return int(ts)
                except ValueError:
                    # ISO format: 2024-01-15 10:23:45.123
                    from datetime import datetime
                    try:
                        dt = datetime.strptime(ts[:23], "%Y-%m-%d %H:%M:%S.%f")
                        return int(dt.timestamp() * 1000)
                    except Exception:
                        return None

            records.append({
                "task_id":      get("task_id"),
                "name":         name,
                "tool":         tool,
                "status":       status,
                "duration_s":   duration_s,
                "realtime_s":   realtime_s,
                "cpu_pct":      cpu_pct,
                "peak_rss_gb":  peak_rss_gb,
                "peak_vmem_gb": peak_vmem_gb,
                "rchar_mb":     rchar_mb * 1000,    # convert GB→MB
                "wchar_mb":     wchar_mb * 1000,
                "submit_ms":    to_ms(submit_str),
                "start_ms":     to_ms(start_str),
                "complete_ms":  to_ms(complete_str),
                "group":        TOOL_GROUPS.get(tool, "unknown"),
            })

    logger.info(f"Loaded {len(records)} tasks from {trace_path}")
    return records


def summarise_by_tool(records: List[dict]) -> Dict[str, dict]:
    """
    Aggregate per-task records by tool name.

    Returns dict tool → {
        durations, cpu_pcts, ram_gbs, n_tasks, n_success, n_fail,
        mean_duration_s, max_ram_gb, mean_cpu_pct, group
    }
    """
    from collections import defaultdict
    agg = defaultdict(lambda: {
        "durations": [], "cpu_pcts": [], "ram_gbs": [],
        "n_tasks": 0, "n_success": 0, "n_fail": 0, "group": "unknown",
    })

    for r in records:
        tool = r["tool"]
        d = agg[tool]
        d["n_tasks"] += 1
        d["group"] = r["group"]
        if r["status"].lower() in ("completed", "ok", "cached"):
            d["n_success"] += 1
        else:
            d["n_fail"] += 1
        if r["duration_s"] > 0:
            d["durations"].append(r["duration_s"])
        if r["cpu_pct"] > 0:
            d["cpu_pcts"].append(r["cpu_pct"])
        if r["peak_rss_gb"] > 0:
            d["ram_gbs"].append(r["peak_rss_gb"])

    result = {}
    for tool, d in agg.items():
        result[tool] = {
            **d,
            "mean_duration_s": np.mean(d["durations"]) if d["durations"] else 0,
            "median_duration_s": np.median(d["durations"]) if d["durations"] else 0,
            "max_ram_gb":    np.max(d["ram_gbs"]) if d["ram_gbs"] else 0,
            "mean_ram_gb":   np.mean(d["ram_gbs"]) if d["ram_gbs"] else 0,
            "std_ram_gb":    np.std(d["ram_gbs"]) if len(d["ram_gbs"]) > 1 else 0,
            "mean_cpu_pct":  np.mean(d["cpu_pcts"]) if d["cpu_pcts"] else 0,
        }
    return result


# ---------------------------------------------------------------------------
# Plot 1: Gantt chart
# ---------------------------------------------------------------------------

def plot_gantt(records: List[dict], outdir: Path) -> Optional[Path]:
    """Timeline Gantt chart showing parallelism."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    timed = [r for r in records if r["start_ms"] and r["complete_ms"]]
    if not timed:
        logger.warning("No timestamp data available for Gantt chart")
        return None

    # Normalise to minutes from pipeline start
    t0 = min(r["start_ms"] for r in timed)
    tools_seen = []
    for r in timed:
        if r["tool"] not in tools_seen:
            tools_seen.append(r["tool"])

    tool_y = {tool: i for i, tool in enumerate(tools_seen)}
    n_tools = len(tools_seen)

    fig, ax = plt.subplots(figsize=(14, max(4, n_tools * 0.5 + 1.5)))

    for r in timed:
        start_m = (r["start_ms"] - t0) / 60000
        end_m   = (r["complete_ms"] - t0) / 60000
        dur_m   = max(end_m - start_m, 0.05)
        y       = tool_y[r["tool"]]
        color   = GROUP_COLORS.get(r["group"], GROUP_COLORS["unknown"])
        alpha   = 0.85 if r["status"].lower() in ("completed", "ok", "cached") else 0.4

        bar = mpatches.FancyBboxPatch(
            (start_m, y - 0.35), dur_m, 0.7,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor="white", linewidth=0.5, alpha=alpha,
        )
        ax.add_patch(bar)

    ax.set_yticks(range(n_tools))
    ax.set_yticklabels(tools_seen, fontsize=9)
    ax.set_ylim(-0.7, n_tools - 0.3)
    ax.set_xlabel("Time from pipeline start (minutes)", fontsize=10)
    ax.set_title("Pipeline Execution Timeline", fontsize=13, fontweight="bold", pad=12)

    # Legend
    legend_patches = [
        mpatches.Patch(color=c, label=GROUP_LABELS[g])
        for g, c in GROUP_COLORS.items()
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=8, framealpha=0.9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    plt.tight_layout()

    path = outdir / "metrics_gantt.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Gantt chart saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Plot 2: CPU usage violin plot
# ---------------------------------------------------------------------------

def plot_cpu(records: List[dict], outdir: Path) -> Path:
    """Violin + strip plot of CPU usage per tool."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from collections import defaultdict
    tool_cpus = defaultdict(list)
    for r in records:
        if r["cpu_pct"] > 0:
            tool_cpus[r["tool"]].append(r["cpu_pct"])

    tools = sorted(tool_cpus.keys(), key=lambda t: np.median(tool_cpus[t]), reverse=True)
    if not tools:
        logger.warning("No CPU data for violin plot")
        return None

    fig, ax = plt.subplots(figsize=(max(8, len(tools) * 1.1), 6))

    data = [tool_cpus[t] for t in tools]
    colors = [GROUP_COLORS.get(TOOL_GROUPS.get(t, "unknown"), GROUP_COLORS["unknown"]) for t in tools]

    # Draw violin for tools with >1 sample
    violin_data = [(d, t, c) for d, t, c in zip(data, tools, colors) if len(d) > 1]
    if violin_data:
        vparts = ax.violinplot(
            [d for d, _, _ in violin_data],
            positions=list(range(len(violin_data))),
            showmedians=True, showextrema=True, widths=0.6,
        )
        for patch, (_, _, color) in zip(vparts["bodies"], violin_data):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        vparts["cmedians"].set_color("white")
        vparts["cmedians"].set_linewidth(2)
        for key in ("cmins", "cmaxes", "cbars"):
            vparts[key].set_color("#555555")
            vparts[key].set_linewidth(0.8)
        vio_tools = [t for _, t, _ in violin_data]
    else:
        vio_tools = []

    # Draw scatter points for all tools
    for i, (t, vals, c) in enumerate(zip(tools, data, colors)):
        xi = i if t in vio_tools else i
        jitter = np.random.default_rng(hash(t) % 2**31).normal(0, 0.06, len(vals))
        ax.scatter([xi + j for j in jitter], vals, color=c, alpha=0.8, s=30, zorder=5, edgecolors="white", linewidths=0.3)

    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels(tools, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("CPU usage (%)", fontsize=10)
    ax.set_title("CPU Usage per Tool", fontsize=13, fontweight="bold", pad=12)
    ax.axhline(100, color="gray", linestyle="--", alpha=0.5, label="100% (1 core)")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    path = outdir / "metrics_cpu.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"CPU plot saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Plot 3: Peak RAM bar chart
# ---------------------------------------------------------------------------

def plot_ram(records: List[dict], outdir: Path, system_ram_gb: float = 14.0) -> Path:
    """Grouped bar chart of peak RAM per tool."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = summarise_by_tool(records)
    tools = sorted(
        [t for t in summary if summary[t]["max_ram_gb"] > 0],
        key=lambda t: summary[t]["max_ram_gb"], reverse=True,
    )
    if not tools:
        logger.warning("No RAM data for bar chart")
        return None

    fig, ax = plt.subplots(figsize=(max(8, len(tools) * 1.1), 6))

    x = np.arange(len(tools))
    means  = [summary[t]["mean_ram_gb"] for t in tools]
    stds   = [summary[t]["std_ram_gb"]  for t in tools]
    maxes  = [summary[t]["max_ram_gb"]  for t in tools]
    colors = [GROUP_COLORS.get(summary[t]["group"], GROUP_COLORS["unknown"]) for t in tools]

    bars = ax.bar(x, means, yerr=stds, capsize=4, color=colors, edgecolor="white",
                  linewidth=0.5, alpha=0.85, error_kw={"elinewidth": 1.5, "ecolor": "#555555"})
    # Overlay max as scatter
    ax.scatter(x, maxes, marker="^", color="black", s=50, zorder=5, label="Peak per run")

    # Label mean above bars
    for xi, mean, std, col in zip(x, means, stds, colors):
        ax.text(xi, mean + std + 0.05, f"{mean:.1f}", ha="center", va="bottom", fontsize=7)

    ax.axhline(system_ram_gb, color="#E63946", linestyle="--", linewidth=1.5,
               label=f"System RAM ({system_ram_gb:.0f} GB)")

    ax.set_xticks(x)
    ax.set_xticklabels(tools, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Peak RSS (GB)", fontsize=10)
    ax.set_title("Peak RAM Usage per Tool", fontsize=13, fontweight="bold", pad=12)
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    path = outdir / "metrics_ram.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"RAM plot saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Plot 4: Runtime distribution (horizontal bar)
# ---------------------------------------------------------------------------

def _fmt_duration(s: float) -> str:
    """Format seconds to human-readable string."""
    if s < 60:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{int(s // 60)}m {int(s % 60)}s"
    return f"{int(s // 3600)}h {int((s % 3600) // 60)}m"


def plot_runtime(records: List[dict], outdir: Path) -> Path:
    """Horizontal bar chart of mean runtime per tool (sorted)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    summary = summarise_by_tool(records)
    tools = sorted(
        [t for t in summary if summary[t]["mean_duration_s"] > 0],
        key=lambda t: summary[t]["mean_duration_s"],
    )
    if not tools:
        logger.warning("No runtime data for horizontal bar chart")
        return None

    means = [summary[t]["mean_duration_s"] for t in tools]
    max_t = max(means)

    # Color gradient: green (fast) → red (slow)
    cmap = LinearSegmentedColormap.from_list("speed", ["#2DC653", "#F9C74F", "#E63946"])
    colors = [cmap(m / max_t) for m in means]

    fig, ax = plt.subplots(figsize=(10, max(4, len(tools) * 0.5 + 2)))
    y = np.arange(len(tools))
    bars = ax.barh(y, means, color=colors, edgecolor="white", linewidth=0.5, height=0.65)

    # Duration labels at bar end
    for yi, mean in zip(y, means):
        ax.text(mean + max_t * 0.01, yi, _fmt_duration(mean),
                va="center", ha="left", fontsize=9, color="#333333")

    ax.set_yticks(y)
    ax.set_yticklabels(tools, fontsize=9)
    ax.set_xlabel("Mean runtime (seconds)", fontsize=10)
    ax.set_title("Runtime Distribution per Tool", fontsize=13, fontweight="bold", pad=12)

    # Secondary x-axis in minutes
    ax2 = ax.twiny()
    ax2.set_xlim(0, ax.get_xlim()[1] / 60)
    ax2.set_xlabel("Minutes", fontsize=9, color="#777777")
    ax2.tick_params(axis="x", colors="#777777", labelsize=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()

    path = outdir / "metrics_runtime.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Runtime plot saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Plot 5: Resource efficiency scatter
# ---------------------------------------------------------------------------

def plot_efficiency(records: List[dict], outdir: Path) -> Path:
    """Scatter plot: runtime vs RAM, with CPU% as bubble size."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = summarise_by_tool(records)
    tools = [
        t for t in summary
        if summary[t]["mean_duration_s"] > 0 and summary[t]["mean_ram_gb"] > 0
    ]
    if not tools:
        logger.warning("Insufficient data for efficiency scatter")
        return None

    x = np.array([summary[t]["mean_duration_s"] / 60 for t in tools])  # minutes
    y = np.array([summary[t]["mean_ram_gb"]            for t in tools])  # GB
    s = np.array([max(summary[t]["mean_cpu_pct"], 10)  for t in tools])  # CPU% → size
    colors = [GROUP_COLORS.get(summary[t]["group"], GROUP_COLORS["unknown"]) for t in tools]

    fig, ax = plt.subplots(figsize=(10, 7))

    scatter = ax.scatter(x, y, s=s * 2.5, c=colors, alpha=0.80,
                         edgecolors="white", linewidths=0.8, zorder=5)

    # Labels with offset to avoid overlap
    for xi, yi, t in zip(x, y, tools):
        ax.annotate(t, (xi, yi), fontsize=7.5, ha="center",
                    xytext=(0, 10), textcoords="offset points",
                    color="#333333")

    # Quadrant reference lines
    med_x, med_y = np.median(x), np.median(y)
    ax.axvline(med_x, color="#AAAAAA", linestyle=":", linewidth=1, alpha=0.7)
    ax.axhline(med_y, color="#AAAAAA", linestyle=":", linewidth=1, alpha=0.7)

    # Quadrant labels
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    quad_kw = dict(fontsize=8, color="#888888", fontstyle="italic", alpha=0.8)
    ax.text(0.02, 0.97, "Fast + Light",  transform=ax.transAxes, va="top", **quad_kw)
    ax.text(0.98, 0.97, "Slow + Light",  transform=ax.transAxes, va="top", ha="right", **quad_kw)
    ax.text(0.02, 0.03, "Fast + Heavy",  transform=ax.transAxes, va="bottom", **quad_kw)
    ax.text(0.98, 0.03, "Slow + Heavy",  transform=ax.transAxes, va="bottom", ha="right", **quad_kw)

    ax.set_xlabel("Mean runtime (minutes)", fontsize=10)
    ax.set_ylabel("Mean peak RAM (GB)", fontsize=10)
    ax.set_title("Resource Efficiency (bubble size = CPU%)",
                 fontsize=13, fontweight="bold", pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.25)
    plt.tight_layout()

    path = outdir / "metrics_efficiency.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Efficiency scatter saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Interactive HTML report (Plotly)
# ---------------------------------------------------------------------------

def generate_html_report(
    records: List[dict],
    outdir: Path,
    system_ram_gb: float = 14.0,
    trace_path: str = "",
) -> Path:
    """Generate interactive Plotly HTML report."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import plotly.io as pio
    except ImportError:
        logger.warning("Plotly not installed — skipping HTML report. Install with: pip install plotly")
        return None

    summary = summarise_by_tool(records)
    tools = sorted(summary.keys(), key=lambda t: summary[t]["mean_duration_s"], reverse=True)

    # ------- Fig 1: Runtime bar chart -------
    fig1 = go.Figure()
    for g, label in GROUP_LABELS.items():
        g_tools = [t for t in tools if summary[t]["group"] == g and summary[t]["mean_duration_s"] > 0]
        if not g_tools:
            continue
        fig1.add_trace(go.Bar(
            name=label,
            x=g_tools,
            y=[summary[t]["mean_duration_s"] / 60 for t in g_tools],
            marker_color=GROUP_COLORS[g],
            error_y=dict(
                type="data",
                array=[np.std(summary[t]["durations"]) / 60 for t in g_tools],
                visible=True,
            ),
            hovertemplate="<b>%{x}</b><br>Mean: %{y:.1f} min<extra></extra>",
        ))
    fig1.update_layout(
        title="Mean Runtime per Tool",
        xaxis_title="Tool", yaxis_title="Runtime (minutes)",
        barmode="group", template="plotly_white",
        legend=dict(orientation="h", y=-0.25),
    )

    # ------- Fig 2: RAM bar chart -------
    fig2 = go.Figure()
    for g, label in GROUP_LABELS.items():
        g_tools = [t for t in tools if summary[t]["group"] == g and summary[t]["mean_ram_gb"] > 0]
        if not g_tools:
            continue
        fig2.add_trace(go.Bar(
            name=label,
            x=g_tools,
            y=[summary[t]["mean_ram_gb"] for t in g_tools],
            marker_color=GROUP_COLORS[g],
            hovertemplate="<b>%{x}</b><br>Mean RAM: %{y:.2f} GB<extra></extra>",
        ))
    fig2.add_hline(y=system_ram_gb, line_dash="dash", line_color="#E63946",
                   annotation_text=f"System RAM ({system_ram_gb:.0f} GB)")
    fig2.update_layout(
        title="Peak RAM per Tool",
        xaxis_title="Tool", yaxis_title="Peak RSS (GB)",
        barmode="group", template="plotly_white",
        legend=dict(orientation="h", y=-0.25),
    )

    # ------- Fig 3: CPU box plot -------
    from collections import defaultdict
    tool_cpus = defaultdict(list)
    for r in records:
        if r["cpu_pct"] > 0:
            tool_cpus[r["tool"]].append(r["cpu_pct"])

    fig3 = go.Figure()
    for t in tools:
        if not tool_cpus[t]:
            continue
        g = summary[t]["group"]
        fig3.add_trace(go.Box(
            y=tool_cpus[t],
            name=t,
            marker_color=GROUP_COLORS.get(g, GROUP_COLORS["unknown"]),
            boxmean=True,
        ))
    fig3.add_hline(y=100, line_dash="dash", line_color="gray", annotation_text="100% (1 core)")
    fig3.update_layout(
        title="CPU Usage per Tool",
        xaxis_title="Tool", yaxis_title="CPU %",
        template="plotly_white", showlegend=False,
    )

    # ------- Fig 4: Efficiency scatter -------
    eff_tools = [t for t in tools if summary[t]["mean_duration_s"] > 0 and summary[t]["mean_ram_gb"] > 0]
    fig4 = go.Figure()
    for g, label in GROUP_LABELS.items():
        g_tools = [t for t in eff_tools if summary[t]["group"] == g]
        if not g_tools:
            continue
        fig4.add_trace(go.Scatter(
            x=[summary[t]["mean_duration_s"] / 60 for t in g_tools],
            y=[summary[t]["mean_ram_gb"] for t in g_tools],
            mode="markers+text",
            name=label,
            marker=dict(
                size=[max(summary[t]["mean_cpu_pct"] / 10, 8) for t in g_tools],
                color=GROUP_COLORS[g],
                opacity=0.8,
                line=dict(color="white", width=1),
            ),
            text=g_tools,
            textposition="top center",
            textfont=dict(size=10),
            hovertemplate="<b>%{text}</b><br>Runtime: %{x:.1f} min<br>RAM: %{y:.2f} GB<extra></extra>",
        ))
    fig4.update_layout(
        title="Resource Efficiency (bubble size ∝ CPU%)",
        xaxis_title="Mean runtime (minutes)",
        yaxis_title="Mean peak RAM (GB)",
        template="plotly_white",
    )

    # ------- Summary table -------
    table_rows = []
    for t in sorted(tools):
        s = summary[t]
        row = {
            "Tool": t,
            "Group": GROUP_LABELS.get(s["group"], s["group"]),
            "Tasks": s["n_tasks"],
            "Success": s["n_success"],
            "Failed": s["n_fail"],
            "Mean runtime": _fmt_duration(s["mean_duration_s"]),
            "Max RAM (GB)": f"{s['max_ram_gb']:.2f}",
            "Mean CPU %": f"{s['mean_cpu_pct']:.0f}",
        }
        table_rows.append(row)

    # Build HTML
    figs_html = []
    for fig in [fig1, fig2, fig3, fig4]:
        figs_html.append(pio.to_html(fig, include_plotlyjs="cdn" if not figs_html else False,
                                     full_html=False, div_id=f"fig{len(figs_html)+1}"))

    table_header = "".join(f"<th>{k}</th>" for k in table_rows[0].keys()) if table_rows else ""
    table_body   = "\n".join(
        "<tr>" + "".join(f"<td>{v}</td>" for v in row.values()) + "</tr>"
        for row in table_rows
    )

    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>HLA Pipeline Execution Metrics</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; color: #333; }}
    h1   {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
    h2   {{ color: #34495e; margin-top: 30px; }}
    .info {{ background: #ecf0f1; padding: 10px 15px; border-radius: 6px; font-size: 13px; margin-bottom: 20px; }}
    .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    .chart {{ background: white; border-radius: 8px; padding: 15px;
               box-shadow: 0 2px 8px rgba(0,0,0,.08); }}
    table {{ border-collapse: collapse; width: 100%; background: white;
              border-radius: 8px; overflow: hidden;
              box-shadow: 0 2px 8px rgba(0,0,0,.08); }}
    th {{ background: #3498db; color: white; padding: 10px 12px; text-align: left; font-size: 13px; }}
    td {{ padding: 8px 12px; font-size: 12px; border-bottom: 1px solid #ecf0f1; }}
    tr:hover td {{ background: #f1f8ff; }}
    tr:nth-child(even) td {{ background: #f9f9f9; }}
    .footer {{ margin-top: 30px; text-align: center; color: #aaa; font-size: 11px; }}
  </style>
</head>
<body>
  <h1>HLA Typing Pipeline — Execution Metrics Report</h1>
  <div class="info">
    <strong>Generated:</strong> {timestamp} &nbsp;|&nbsp;
    <strong>Trace file:</strong> {os.path.basename(trace_path) if trace_path else "N/A"} &nbsp;|&nbsp;
    <strong>Tasks:</strong> {len(records)} &nbsp;|&nbsp;
    <strong>Tools:</strong> {len(summary)}
  </div>

  <h2>Interactive Charts</h2>
  <div class="grid">
    <div class="chart">{figs_html[0]}</div>
    <div class="chart">{figs_html[1]}</div>
    <div class="chart">{figs_html[2]}</div>
    <div class="chart">{figs_html[3]}</div>
  </div>

  <h2>Summary Table</h2>
  <table>
    <thead><tr>{table_header}</tr></thead>
    <tbody>{table_body}</tbody>
  </table>

  <div class="footer">
    Generated by hla_pipeline_metrics.py &mdash; HLA Typing Pipeline v1.3.0
  </div>
</body>
</html>
"""

    path = outdir / "execution_metrics_report.html"
    path.write_text(html)
    logger.info(f"HTML report saved: {path}")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def create_demo_trace(path: Path) -> None:
    """Create a small synthetic trace file for testing."""
    import time
    now_ms = int(time.time() * 1000)
    rows = [
        "task_id\thash\tname\tstatus\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\tsubmit\tstart\tcomplete",
        f"1\ta1b2c3d4\tCHECK_BAM (FH_7087)\tCOMPLETED\t45s\t43s\t98.5%\t320 MB\t512 MB\t120 MB\t5 MB\t{now_ms}\t{now_ms+100}\t{now_ms+45100}",
        f"2\tb2c3d4e5\tEXTRACT_HLA_READS (FH_7087)\tCOMPLETED\t2m 15s\t2m 12s\t380.2%\t1.2 GB\t2.1 GB\t2.5 GB\t45 MB\t{now_ms+100}\t{now_ms+200}\t{now_ms+135200}",
        f"3\tc3d4e5f6\tSPECHLA (FH_7087)\tCOMPLETED\t4m 38s\t4m 35s\t420.0%\t8.5 GB\t12.1 GB\t350 MB\t120 MB\t{now_ms+135000}\t{now_ms+135100}\t{now_ms+413100}",
        f"4\td4e5f6g7\tHLAHD (FH_7087)\tCOMPLETED\t3m 52s\t3m 49s\t350.5%\t10.2 GB\t14.5 GB\t280 MB\t95 MB\t{now_ms+135000}\t{now_ms+135100}\t{now_ms+367100}",
        f"5\te5f6g7h8\tARCASHLA (FH_7087)\tCOMPLETED\t1m 45s\t1m 43s\t290.0%\t6.1 GB\t9.8 GB\t180 MB\t62 MB\t{now_ms+135000}\t{now_ms+135100}\t{now_ms+240100}",
        f"6\tf6g7h8i9\tXHLA (FH_7087)\tCOMPLETED\t6m 12s\t6m 10s\t180.0%\t5.8 GB\t8.5 GB\t1.2 GB\t85 MB\t{now_ms+135000}\t{now_ms+135100}\t{now_ms+507100}",
        f"7\tg7h8i9j0\tCONSENSUS (FH_7087)\tCOMPLETED\t12s\t11s\t95.2%\t512 MB\t768 MB\t25 MB\t8 MB\t{now_ms+510000}\t{now_ms+510100}\t{now_ms+522100}",
        f"8\th8i9j0k1\tHLA_VISUALIZE (FH_7087)\tCOMPLETED\t25s\t24s\t110.3%\t680 MB\t1.1 GB\t18 MB\t12 MB\t{now_ms+522000}\t{now_ms+522100}\t{now_ms+547100}",
        f"9\ti9j0k1l2\tHLA_SUMMARY_REPORT (FH_7087)\tCOMPLETED\t18s\t17s\t102.1%\t590 MB\t900 MB\t15 MB\t9 MB\t{now_ms+547000}\t{now_ms+547100}\t{now_ms+565100}",
        f"10\tj0k1l2m3\tQC_BAM (FH_7087)\tCOMPLETED\t38s\t37s\t195.0%\t280 MB\t450 MB\t95 MB\t4 MB\t{now_ms}\t{now_ms+50}\t{now_ms+38050}",
    ]
    path.write_text("\n".join(rows) + "\n")
    logger.info(f"Demo trace written to {path}")


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="HLA Pipeline Execution Metrics Visualizer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--trace", help="Nextflow trace file (or glob pattern)")
    parser.add_argument("--outdir", required=True, help="Output directory for plots and report")
    parser.add_argument("--system-ram", type=float, default=14.0,
                        help="System RAM limit in GB (default: 14.0, shown as reference line)")
    parser.add_argument("--demo", action="store_true",
                        help="Generate demo output with synthetic data (no trace file needed)")
    parser.add_argument("--no-plotly", action="store_true",
                        help="Skip Plotly HTML report (matplotlib plots only)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        demo_trace = outdir / "_demo_trace.txt"
        create_demo_trace(demo_trace)
        args.trace = str(demo_trace)

    if not args.trace:
        parser.error("--trace required (or use --demo for synthetic demo output)")

    # Resolve glob
    import glob as _glob
    trace_files = _glob.glob(args.trace)
    if not trace_files:
        parser.error(f"No trace files found matching: {args.trace}")

    # Load all matching trace files (multi-run support)
    all_records = []
    for tf in sorted(trace_files):
        try:
            all_records.extend(load_trace(tf))
        except Exception as e:
            logger.error(f"Failed to load {tf}: {e}")

    if not all_records:
        logger.error("No valid records loaded from trace files")
        sys.exit(1)

    logger.info(f"Total tasks: {len(all_records)} from {len(trace_files)} trace file(s)")

    # Generate plots
    plots_generated = []

    try:
        p = plot_gantt(all_records, outdir)
        if p:
            plots_generated.append(p)
    except Exception as e:
        logger.warning(f"Gantt chart failed: {e}")

    try:
        p = plot_cpu(all_records, outdir)
        if p:
            plots_generated.append(p)
    except Exception as e:
        logger.warning(f"CPU plot failed: {e}")

    try:
        p = plot_ram(all_records, outdir, system_ram_gb=args.system_ram)
        if p:
            plots_generated.append(p)
    except Exception as e:
        logger.warning(f"RAM plot failed: {e}")

    try:
        p = plot_runtime(all_records, outdir)
        if p:
            plots_generated.append(p)
    except Exception as e:
        logger.warning(f"Runtime plot failed: {e}")

    try:
        p = plot_efficiency(all_records, outdir)
        if p:
            plots_generated.append(p)
    except Exception as e:
        logger.warning(f"Efficiency scatter failed: {e}")

    if not args.no_plotly:
        try:
            trace_path = trace_files[0] if trace_files else ""
            p = generate_html_report(all_records, outdir,
                                     system_ram_gb=args.system_ram,
                                     trace_path=trace_path)
            if p:
                plots_generated.append(p)
        except Exception as e:
            logger.warning(f"HTML report failed: {e}")

    logger.info(f"\nDone! Generated {len(plots_generated)} output files in {outdir}:")
    for p in plots_generated:
        logger.info(f"  {p}")


if __name__ == "__main__":
    main()
