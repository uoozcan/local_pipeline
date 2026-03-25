#!/usr/bin/env python3
"""make_paper_figures.py — Generate publication-quality scientific figures
from HLA typing pipeline outputs.

Produces 16 figures in 4 groups:
  A: Resource benchmarking      (from Nextflow trace)
  B: HLA typing quality         (from consensus/comparison TSVs)
  C: Tool accuracy              (from tool_weights JSON)
  D: Summary panels             (composite)

Usage:
    python3 make_paper_figures.py \
        --consensus  results/NA19238/NA19238_consensus.txt \
        --comparison results/NA19238/NA19238_comparison.txt \
        --trace_file results/pipeline_info/trace_2026-01-01_00-00-00.txt \
        --weights_file conf/tool_weights_wgs_v3.json \
        --sample_id NA19238 \
        --output_dir results/NA19238/scientific_figures/ \
        --format both
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

# ─── Constants ────────────────────────────────────────────────────────────────

GENE_ORDER = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1',
              'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1']

TOOL_ORDER = ['HLAHD', 'SPECHLA', 'OPTITYPE', 'ARCASHLA',
              'KOURAMI', 'POLYSOLVER', 'SEQ2HLA', 'XHLA']

TOOL_COLORS = {
    'HLAHD':      '#4472C4',
    'SPECHLA':    '#ED7D31',
    'OPTITYPE':   '#E74C3C',
    'ARCASHLA':   '#2ECC71',
    'KOURAMI':    '#7030A0',
    'POLYSOLVER': '#00B0F0',
    'SEQ2HLA':    '#FF99CC',
    'XHLA':       '#92D050',
    'FAILED':     '#BFBFBF',
    'QC':         '#A9A9A9',
    'FASTQC':     '#A9A9A9',
    'CONSENSUS':  '#A9A9A9',
}

# Puhti configured memory limits per tool (GB)
CONFIGURED_MEM = {
    'HLAHD': 64, 'SPECHLA': 64, 'ARCASHLA': 32, 'OPTITYPE': 32,
    'KOURAMI': 16, 'POLYSOLVER': 16, 'SEQ2HLA': 8, 'XHLA': 32,
}

CONF_GREEN = '#27AE60'
CONF_AMBER = '#E67E22'
CONF_RED   = '#E74C3C'

FIGURE_DPI = 300

HLA_TOOL_RE = re.compile(
    r'^(HLAHD|SPECHLA|ARCASHLA|OPTITYPE|KOURAMI|POLYSOLVER|SEQ2HLA|XHLA'
    r'|T1K|HIFIHLA|HLALA|HLASCAN)',
    re.IGNORECASE
)


# ─── Style ────────────────────────────────────────────────────────────────────

def setup_style():
    plt.rcParams.update({
        'font.family': 'DejaVu Serif',
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'savefig.dpi': FIGURE_DPI,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'axes.spines.top': False,
        'axes.spines.right': False,
    })
    if HAS_SEABORN:
        sns.set_style('whitegrid')
        sns.set_context('paper', font_scale=1.0)


# ─── Parsing ──────────────────────────────────────────────────────────────────

def parse_duration(s):
    """'1h 17m 40s' / '7m 46s' / '1.3s' / '-'  →  seconds (float) or None."""
    if not s or s.strip() in ('-', ''):
        return None
    s = s.strip()
    total = 0.0
    for pattern, factor in [(r'(\d+)h', 3600), (r'(\d+)m', 60), (r'([\d.]+)s', 1)]:
        m = re.search(pattern, s)
        if m:
            total += float(m.group(1)) * factor
    return total if total > 0 else None


def parse_memory(s):
    """'2.3 GB' / '359.7 MB' / '-'  →  GB (float) or None."""
    if not s or s.strip() in ('-', ''):
        return None
    m = re.match(r'([\d.]+)\s*(B|KB|MB|GB|TB)', s.strip(), re.I)
    if not m:
        return None
    factors = {'B': 1e-9, 'KB': 1e-6, 'MB': 1e-3, 'GB': 1.0, 'TB': 1e3}
    return float(m.group(1)) * factors[m.group(2).upper()]


def parse_cpu(s):
    """'476.6%' / '-'  →  float or None."""
    if not s or s.strip() in ('-', ''):
        return None
    return float(s.strip().rstrip('%'))


def parse_trace(trace_file):
    """Return DataFrame with HLA tool rows and numeric columns added."""
    df = pd.read_csv(trace_file, sep='\t', dtype=str)
    # Filter to HLA typing tool rows only
    df = df[df['name'].str.match(HLA_TOOL_RE, na=False)].copy()
    if df.empty:
        return df

    def tool_name(n):
        m = re.match(r'^([A-Z][A-Z0-9_]+)', n.strip())
        return m.group(1).upper() if m else n.split('(')[0].strip().upper()

    df['tool']       = df['name'].apply(tool_name)
    df['duration_s'] = df['duration'].apply(parse_duration)
    df['peak_rss_gb']= df['peak_rss'].apply(parse_memory)
    df['cpu_pct']    = df['%cpu'].apply(parse_cpu)
    df['rchar_gb']   = df['rchar'].apply(parse_memory)
    df['wchar_gb']   = df['wchar'].apply(parse_memory)
    df['success']    = df['status'].str.upper().isin(['COMPLETED', 'CACHED'])
    df['submit_dt']  = pd.to_datetime(df['submit'], errors='coerce')
    return df


def _dedup_trace(df):
    """Keep the most recent completed run per tool (remove retries)."""
    completed = df[df['success']].copy()
    failed    = df[~df['success'] & ~df['tool'].isin(completed['tool'])].copy()
    best = pd.concat([
        completed.sort_values('duration_s', na_position='last').drop_duplicates('tool', keep='last'),
        failed.drop_duplicates('tool', keep='last'),
    ], ignore_index=True)
    return best


def parse_consensus(path):
    """Return DataFrame with columns Gene, Allele1, Allele2, Confidence, Reads1, Reads2."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts[0] == 'Gene':
                continue
            if len(parts) < 4:
                continue
            def safe_int(v):
                try: return int(v)
                except: return 0
            def safe_float(v):
                try: return float(v)
                except: return None
            rows.append({
                'Gene':       parts[0],
                'Allele1':    parts[1] if len(parts) > 1 else '-',
                'Allele2':    parts[2] if len(parts) > 2 else '-',
                'Confidence': safe_float(parts[3]) if len(parts) > 3 else None,
                'Reads1':     safe_int(parts[4]) if len(parts) > 4 else 0,
                'Reads2':     safe_int(parts[5]) if len(parts) > 5 else 0,
            })
    return pd.DataFrame(rows)


def parse_comparison(path):
    """Return (gene_data dict, tool_names list).

    gene_data = {
        'HLA-A': {
            'consensus': 'A*03:01/A*11:01',
            'tools': {'HLAHD': {'alleles': 'A*03:01', 'reads': 49}, ...}
        }, ...
    }
    """
    gene_data  = {}
    tool_names = []
    header     = None

    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts[0] == 'Gene':
                header = parts
                tool_names = [
                    re.match(r'^(.+)\(alleles\)$', h).group(1).strip().upper()
                    for h in header
                    if re.match(r'^(.+)\(alleles\)$', h)
                ]
                continue
            if header is None:
                continue
            row = dict(zip(header, parts))
            gene = row.get('Gene', '').strip()
            if not gene:
                continue
            gene_data[gene] = {
                'consensus': row.get('Consensus', '-'),
                'tools': {}
            }
            for tool in tool_names:
                alleles = '-'
                reads   = 0
                for col, val in row.items():
                    cl = col.lower()
                    tl = tool.lower()
                    if cl == f'{tl}(alleles)':
                        alleles = val.strip() or '-'
                    elif cl == f'{tl}(reads)':
                        try:
                            reads = int(val)
                        except (ValueError, TypeError):
                            reads = 0
                gene_data[gene]['tools'][tool] = {'alleles': alleles, 'reads': reads}

    return gene_data, tool_names


def load_weights(path):
    with open(path) as f:
        return json.load(f)


# ─── Save helper ──────────────────────────────────────────────────────────────

def _save(fig, output_dir, name, fmt):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if fmt in ('png', 'both'):
        fig.savefig(Path(output_dir) / f'{name}.png',
                    dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    if fmt in ('pdf', 'both'):
        fig.savefig(Path(output_dir) / f'{name}.pdf',
                    bbox_inches='tight', facecolor='white')


# ─── Group A: Resource Benchmarking ──────────────────────────────────────────

def fig_A1_runtime(trace_df, output_dir, fmt):
    """Horizontal bar chart: wall-clock time per HLA tool."""
    df = _dedup_trace(trace_df).dropna(subset=['duration_s'])
    if df.empty:
        return
    df = df.sort_values('duration_s')

    fig, ax = plt.subplots(figsize=(7, max(3, 0.55 * len(df))))
    colors = [TOOL_COLORS.get(t, '#888') if s else TOOL_COLORS['FAILED']
              for t, s in zip(df['tool'], df['success'])]
    bars = ax.barh(df['tool'], df['duration_s'] / 60,
                   color=colors, edgecolor='white', height=0.6)

    for bar, row in zip(bars, df.itertuples()):
        w = bar.get_width()
        secs = row.duration_s
        lbl = (f"{int(secs//3600)}h {int((secs%3600)//60)}m"
               if secs >= 3600 else f"{secs/60:.1f} min")
        ax.text(w * 1.01 + 0.1, bar.get_y() + bar.get_height() / 2,
                lbl, va='center', fontsize=8)

    ax.set_xlabel('Wall-clock time (minutes)')
    ax.set_title('Tool Execution Time')
    ax.set_xlim(0, df['duration_s'].max() / 60 * 1.3)
    failed = df[~df['success']]
    if not failed.empty:
        ax.text(0.99, 0.01,
                f"Grey = FAILED ({', '.join(failed['tool'])})",
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=7, color='grey', style='italic')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A1_runtime', fmt)
    plt.close(fig)


def fig_A2_peak_memory(trace_df, output_dir, fmt):
    """Grouped bar: peak RSS vs configured memory limit."""
    df = _dedup_trace(trace_df).dropna(subset=['peak_rss_gb'])
    if df.empty:
        return
    df = df.sort_values('peak_rss_gb', ascending=False)

    x = np.arange(len(df))
    w = 0.35
    fig, ax = plt.subplots(figsize=(max(5, len(df) * 0.9), 4.5))
    colors = [TOOL_COLORS.get(t, '#888') for t in df['tool']]
    ax.bar(x - w/2, df['peak_rss_gb'], w,
           label='Peak RSS (actual)', color=colors, alpha=0.85)
    limits = [CONFIGURED_MEM.get(t, 32) for t in df['tool']]
    ax.bar(x + w/2, limits, w,
           label='Configured limit', color='lightgrey',
           edgecolor='#555', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(df['tool'], rotation=30, ha='right')
    ax.set_ylabel('Memory (GB)')
    ax.set_title('Peak Memory Usage vs Configured Limit')
    ax.legend()
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A2_peak_memory', fmt)
    plt.close(fig)


def fig_A3_cpu_utilization(trace_df, output_dir, fmt):
    """Bar chart: % CPU per tool."""
    df = _dedup_trace(trace_df).dropna(subset=['cpu_pct'])
    if df.empty:
        return
    df = df.sort_values('cpu_pct', ascending=False)

    fig, ax = plt.subplots(figsize=(max(5, len(df) * 0.9), 4))
    colors = [TOOL_COLORS.get(t, '#888') for t in df['tool']]
    ax.bar(df['tool'], df['cpu_pct'], color=colors, edgecolor='white')
    ax.axhline(100, color='black', linestyle='--', linewidth=0.8,
               label='100% (1 core)')
    ax.set_ylabel('CPU utilization (%)')
    ax.set_title('CPU Utilization per Tool (multi-core)')
    ax.set_xticklabels(df['tool'], rotation=30, ha='right')
    ax.legend()
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A3_cpu_utilization', fmt)
    plt.close(fig)


def fig_A4_io_throughput(trace_df, output_dir, fmt):
    """Grouped bar: bytes read + written per tool."""
    df = _dedup_trace(trace_df).copy()
    df['rchar_gb'] = df['rchar_gb'].fillna(0)
    df['wchar_gb'] = df['wchar_gb'].fillna(0)
    df = df[(df['rchar_gb'] > 0) | (df['wchar_gb'] > 0)].sort_values(
        'rchar_gb', ascending=False)
    if df.empty:
        return

    x  = np.arange(len(df))
    w  = 0.35
    fig, ax = plt.subplots(figsize=(max(5, len(df) * 0.9), 4))
    colors = [TOOL_COLORS.get(t, '#888') for t in df['tool']]
    ax.bar(x - w/2, df['rchar_gb'], w, label='Read (GB)',
           color=colors, alpha=0.85)
    ax.bar(x + w/2, df['wchar_gb'], w, label='Written (GB)',
           color=colors, alpha=0.45, edgecolor='black', linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(df['tool'], rotation=30, ha='right')
    ax.set_ylabel('I/O throughput (GB)')
    ax.set_title('I/O Throughput per Tool')
    ax.legend()
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A4_io_throughput', fmt)
    plt.close(fig)


def fig_A5_gantt(trace_df, output_dir, fmt):
    """Gantt chart of parallel execution."""
    df = _dedup_trace(trace_df).dropna(subset=['submit_dt', 'duration_s']).copy()
    if df.empty:
        return
    df['end_dt'] = df['submit_dt'] + pd.to_timedelta(df['duration_s'], unit='s')
    df = df.sort_values('submit_dt')
    t0 = df['submit_dt'].min()

    fig, ax = plt.subplots(figsize=(9, max(3, 0.6 * len(df))))
    for _, row in df.iterrows():
        start_min = (row['submit_dt'] - t0).total_seconds() / 60
        dur_min   = row['duration_s'] / 60
        color = (TOOL_COLORS.get(row['tool'], '#888')
                 if row['success'] else TOOL_COLORS['FAILED'])
        ax.barh(row['tool'], dur_min, left=start_min,
                color=color, alpha=0.8, height=0.5, edgecolor='white')

    ax.set_xlabel('Time from start (minutes)')
    ax.set_title('Pipeline Execution Timeline (Gantt)')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A5_gantt', fmt)
    plt.close(fig)


def fig_A6_resource_efficiency(trace_df, output_dir, fmt):
    """Scatter: runtime (x) × peak RSS (y), bubble = %CPU."""
    df = _dedup_trace(trace_df).dropna(subset=['duration_s', 'peak_rss_gb']).copy()
    if df.empty:
        return
    df['cpu_pct'] = df['cpu_pct'].fillna(100)

    fig, ax = plt.subplots(figsize=(6, 5))
    for _, row in df.iterrows():
        color = (TOOL_COLORS.get(row['tool'], '#888')
                 if row['success'] else TOOL_COLORS['FAILED'])
        size = max(60, min(600, row['cpu_pct'] * 1.5))
        ax.scatter(row['duration_s'] / 60, row['peak_rss_gb'],
                   s=size, color=color, alpha=0.75,
                   edgecolors='black', linewidths=0.5, zorder=3)
        ax.annotate(row['tool'],
                    (row['duration_s'] / 60, row['peak_rss_gb']),
                    fontsize=8, xytext=(4, 3), textcoords='offset points')

    ax.set_xlabel('Wall-clock time (minutes)')
    ax.set_ylabel('Peak RSS (GB)')
    ax.set_title('Resource Efficiency (bubble size ∝ CPU%)')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_A6_resource_efficiency', fmt)
    plt.close(fig)


# ─── Group B: HLA Typing Quality ─────────────────────────────────────────────

def fig_B1_confidence_per_gene(consensus_df, sample_id, output_dir, fmt):
    """Bar chart per gene, color-coded by confidence threshold."""
    df = consensus_df.dropna(subset=['Confidence']).copy()
    df['Gene_short'] = df['Gene'].str.replace('HLA-', '', regex=False)
    order = [g.replace('HLA-', '') for g in GENE_ORDER if g in df['Gene'].values]
    df = df.set_index('Gene_short').reindex(order).reset_index()
    df = df.dropna(subset=['Confidence'])
    if df.empty:
        return

    colors = [CONF_GREEN if c >= 0.8 else CONF_AMBER if c >= 0.5 else CONF_RED
              for c in df['Confidence']]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(df['Gene_short'], df['Confidence'],
                  color=colors, edgecolor='white', width=0.6)
    ax.axhline(0.8, color=CONF_GREEN, linestyle='--', linewidth=0.8, alpha=0.7,
               label='High ≥ 0.8')
    ax.axhline(0.5, color=CONF_AMBER, linestyle='--', linewidth=0.8, alpha=0.7,
               label='Medium ≥ 0.5')
    ax.set_ylim(0, 1.15)
    ax.set_ylabel('Confidence score')
    ax.set_title(f'Consensus Confidence per HLA Gene — {sample_id}')
    ax.legend(loc='upper right', fontsize=8)
    for bar, val in zip(bars, df['Confidence']):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02, f'{val:.2f}',
                ha='center', va='bottom', fontsize=8)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_B1_confidence_per_gene', fmt)
    plt.close(fig)


def fig_B2_tool_concordance(gene_data, tool_names, sample_id, output_dir, fmt):
    """Heatmap: tools × genes, 1=match / 0=mismatch / NaN=not typed."""
    genes = [g for g in GENE_ORDER if g in gene_data]
    tools = [t for t in TOOL_ORDER if t in tool_names]
    tools += [t for t in tool_names if t not in TOOL_ORDER]

    if not genes or not tools:
        return

    labels_g = [g.replace('HLA-', '') for g in genes]
    matrix   = np.full((len(tools), len(genes)), np.nan)

    for j, gene in enumerate(genes):
        entry     = gene_data[gene]
        consensus = entry.get('consensus', '-')
        for i, tool in enumerate(tools):
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            if alleles in ('-', '', None):
                matrix[i, j] = np.nan
            elif consensus and (alleles == consensus
                                or consensus in alleles
                                or alleles in consensus):
                matrix[i, j] = 1.0
            else:
                matrix[i, j] = 0.0

    fig, ax = plt.subplots(
        figsize=(max(6, len(labels_g) * 0.9), max(3, len(tools) * 0.65)))

    mask = np.isnan(matrix)
    if HAS_SEABORN:
        annot = np.where(mask, '',
                         np.where(matrix == 1.0, '✓', '✗').astype(object))
        sns.heatmap(matrix, annot=annot, fmt='',
                    cmap='RdYlGn', vmin=0, vmax=1,
                    xticklabels=labels_g, yticklabels=tools,
                    linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Agreement with consensus'},
                    mask=mask, ax=ax)
    else:
        im = ax.imshow(matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
        ax.set_xticks(range(len(labels_g)))
        ax.set_yticks(range(len(tools)))
        ax.set_xticklabels(labels_g, rotation=30, ha='right')
        ax.set_yticklabels(tools)
        plt.colorbar(im, ax=ax, label='Agreement')

    ax.set_title(f'Tool vs Consensus Agreement — {sample_id}')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_B2_tool_concordance', fmt)
    plt.close(fig)


def fig_B3_read_support(consensus_df, sample_id, output_dir, fmt):
    """Grouped bar per gene: reads allele 1 vs allele 2."""
    df = consensus_df.copy()
    df['Gene_short'] = df['Gene'].str.replace('HLA-', '', regex=False)
    order = [g.replace('HLA-', '') for g in GENE_ORDER if g in df['Gene'].values]
    df = df.set_index('Gene_short').reindex(order).reset_index()
    if df[['Reads1', 'Reads2']].sum().sum() == 0:
        return  # No read data available

    x = np.arange(len(df))
    w = 0.35
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x - w/2, df['Reads1'].fillna(0), w,
           label='Allele 1', color='#4472C4', alpha=0.85)
    ax.bar(x + w/2, df['Reads2'].fillna(0), w,
           label='Allele 2', color='#ED7D31', alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(df['Gene_short'], rotation=30, ha='right')
    ax.set_ylabel('Read count')
    ax.set_title(f'Read Support per Allele — {sample_id}')
    ax.legend()
    fig.tight_layout()
    _save(fig, output_dir, 'fig_B3_read_support', fmt)
    plt.close(fig)


def fig_B4_confidence_vs_reads(consensus_df, sample_id, output_dir, fmt):
    """Scatter per locus: total reads vs confidence."""
    df = consensus_df.dropna(subset=['Confidence']).copy()
    df['TotalReads'] = df['Reads1'].fillna(0) + df['Reads2'].fillna(0)
    df['Gene_short'] = df['Gene'].str.replace('HLA-', '', regex=False)
    if df['TotalReads'].sum() == 0:
        return

    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(df['TotalReads'], df['Confidence'],
               c=range(len(df)), cmap='tab10',
               s=80, edgecolors='black', linewidths=0.5, zorder=3)
    for _, row in df.iterrows():
        ax.annotate(row['Gene_short'],
                    (row['TotalReads'], row['Confidence']),
                    fontsize=8, xytext=(4, 3), textcoords='offset points')
    ax.axhline(0.8, color=CONF_GREEN, linestyle='--', linewidth=0.8, alpha=0.7)
    ax.axhline(0.5, color=CONF_AMBER, linestyle='--', linewidth=0.8, alpha=0.7)
    ax.set_xlabel('Total reads (allele 1 + 2)')
    ax.set_ylabel('Confidence score')
    ax.set_title(f'Read Support vs Confidence — {sample_id}')
    ax.set_ylim(-0.05, 1.1)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_B4_confidence_vs_reads', fmt)
    plt.close(fig)


# ─── Group C: Tool Accuracy ───────────────────────────────────────────────────

def _raw_accuracy_matrix(weights_data):
    """Return (matrix, tools, genes) where matrix[i,j] is concordance * 100."""
    raw   = weights_data.get('raw_accuracy', {})
    if not raw:
        return None, [], []
    genes = [g for g in ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
             if any(g in v for v in raw.values())]
    tools = sorted(raw.keys(), key=lambda t: t.upper())
    matrix = np.full((len(tools), len(genes)), np.nan)
    for i, tool in enumerate(tools):
        for j, gene in enumerate(genes):
            v = raw[tool].get(gene)
            if v is not None:
                matrix[i, j] = float(v) * 100
    return matrix, [t.upper() for t in tools], genes


def fig_C1_accuracy_heatmap(weights_data, output_dir, fmt):
    """Heatmap: raw concordance tools × genes (%)."""
    matrix, tools, genes = _raw_accuracy_matrix(weights_data)
    if matrix is None:
        return

    n = weights_data.get('n_samples', '?')
    fig, ax = plt.subplots(
        figsize=(max(6, len(genes) * 0.9), max(3, len(tools) * 0.65)))
    mask = np.isnan(matrix)

    if HAS_SEABORN:
        sns.heatmap(matrix, annot=True, fmt='.1f',
                    cmap='YlOrRd', vmin=0, vmax=100,
                    xticklabels=genes, yticklabels=tools,
                    linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Concordance with GT (%)'},
                    mask=mask, ax=ax)
    else:
        display = np.where(mask, 0, matrix)
        im = ax.imshow(display, cmap='YlOrRd', vmin=0, vmax=100, aspect='auto')
        ax.set_xticks(range(len(genes)))
        ax.set_yticks(range(len(tools)))
        ax.set_xticklabels(genes, rotation=30, ha='right')
        ax.set_yticklabels(tools)
        plt.colorbar(im, ax=ax, label='Concordance (%)')

    ax.set_title(f'Tool Accuracy per HLA Gene (N={n} WGS samples)')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_C1_accuracy_heatmap', fmt)
    plt.close(fig)


def fig_C2_weight_heatmap(weights_data, output_dir, fmt):
    """Heatmap: normalized consensus weights tools × genes."""
    genes_dict = weights_data.get('genes', {})
    if not genes_dict:
        return

    genes = [g for g in ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
             if g in genes_dict]
    all_tools = set()
    for g in genes:
        all_tools.update(genes_dict[g].keys())
    tools = sorted(all_tools)

    matrix = np.zeros((len(tools), len(genes)))
    for j, gene in enumerate(genes):
        for i, tool in enumerate(tools):
            matrix[i, j] = float(genes_dict[gene].get(tool, 0.0))

    fig, ax = plt.subplots(
        figsize=(max(6, len(genes) * 0.9), max(3, len(tools) * 0.65)))

    if HAS_SEABORN:
        sns.heatmap(matrix, annot=True, fmt='.2f',
                    cmap='Blues', vmin=0, vmax=1,
                    xticklabels=genes,
                    yticklabels=[t.upper() for t in tools],
                    linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Normalized weight'}, ax=ax)
    else:
        im = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=1, aspect='auto')
        ax.set_xticks(range(len(genes)))
        ax.set_yticks(range(len(tools)))
        ax.set_xticklabels(genes, rotation=30, ha='right')
        ax.set_yticklabels([t.upper() for t in tools])
        plt.colorbar(im, ax=ax, label='Normalized weight')

    ax.set_title('Normalized Tool Weights for Consensus Voting')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_C2_weight_heatmap', fmt)
    plt.close(fig)


def fig_C3_tool_ranking(weights_data, output_dir, fmt):
    """Horizontal bar: mean accuracy per tool (A, B, C, DRB1, DQB1)."""
    raw = weights_data.get('raw_accuracy', {})
    if not raw:
        return

    core = ['A', 'B', 'C', 'DRB1', 'DQB1']
    rows = []
    for tool, gene_acc in raw.items():
        vals = [float(gene_acc[g]) * 100
                for g in core if g in gene_acc and gene_acc[g] is not None]
        if vals:
            rows.append({'tool': tool.upper(), 'mean_acc': np.mean(vals)})

    df = pd.DataFrame(rows).sort_values('mean_acc')
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(6, max(3, 0.55 * len(df))))
    colors = [TOOL_COLORS.get(t, '#888') for t in df['tool']]
    bars = ax.barh(df['tool'], df['mean_acc'],
                   color=colors, edgecolor='white', height=0.6)
    for bar, val in zip(bars, df['mean_acc']):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                f'{val:.1f}%', va='center', fontsize=8)
    ax.set_xlabel('Mean concordance with ground truth (%)')
    ax.set_title('Tool Ranking (A, B, C, DRB1, DQB1)')
    ax.set_xlim(0, 110)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_C3_tool_ranking', fmt)
    plt.close(fig)


def fig_C4_accuracy_vs_runtime(weights_data, trace_df, output_dir, fmt):
    """Scatter: mean concordance (y) vs median runtime from trace (x)."""
    raw = weights_data.get('raw_accuracy', {})
    if not raw or trace_df.empty:
        return

    core = ['A', 'B', 'C', 'DRB1', 'DQB1']
    means = {}
    for tool, gene_acc in raw.items():
        vals = [float(gene_acc[g]) * 100
                for g in core if g in gene_acc and gene_acc[g] is not None]
        if vals:
            means[tool.upper()] = np.mean(vals)

    runtime_med = (trace_df.dropna(subset=['duration_s'])
                   .groupby('tool')['duration_s'].median())

    plot_data = [(t, acc, runtime_med.get(t))
                 for t, acc in means.items()
                 if runtime_med.get(t) is not None]
    if not plot_data:
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    for tool, acc, rt in plot_data:
        color = TOOL_COLORS.get(tool, '#888')
        ax.scatter(rt / 60, acc, s=120, color=color,
                   edgecolors='black', linewidths=0.5, zorder=3)
        ax.annotate(tool, (rt / 60, acc),
                    fontsize=8, xytext=(4, 3), textcoords='offset points')

    ax.set_xlabel('Median runtime (minutes)')
    ax.set_ylabel('Mean concordance (%)')
    ax.set_title('Accuracy vs Runtime Trade-off per Tool')
    ax.set_ylim(0, 105)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_C4_accuracy_vs_runtime', fmt)
    plt.close(fig)


# ─── Group D: Summary Panels ──────────────────────────────────────────────────

def fig_D1_tool_calls_table(gene_data, tool_names, consensus_df, sample_id,
                             output_dir, fmt):
    """Visual table: all tool allele calls + consensus per gene."""
    genes     = [g for g in GENE_ORDER if g in gene_data]
    tools_ord = [t for t in TOOL_ORDER if t in tool_names]
    tools_ord += [t for t in tool_names if t not in TOOL_ORDER]

    if not genes or not tools_ord:
        return

    conf_map = {}
    if consensus_df is not None and not consensus_df.empty:
        for _, row in consensus_df.iterrows():
            conf_map[row['Gene']] = row.get('Confidence')

    col_labels = ['Gene', 'Consensus'] + tools_ord
    cell_text  = []
    cell_colors = []

    for gene in genes:
        entry     = gene_data[gene]
        consensus = entry.get('consensus', '-')
        conf      = conf_map.get(gene)

        if conf is not None and conf >= 0.8:
            cons_bg = '#D5E8D4'
        elif conf is not None and conf >= 0.5:
            cons_bg = '#FFE6CC'
        else:
            cons_bg = '#F8CECC'

        row_text  = [gene.replace('HLA-', ''), consensus]
        row_color = ['#DAE8FC', cons_bg]

        for tool in tools_ord:
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '—')
            if alleles in ('-', ''):
                alleles = '—'
            row_text.append(alleles)
            if alleles == '—':
                row_color.append('#F5F5F5')
            elif (consensus and
                  (alleles == consensus
                   or consensus in alleles
                   or alleles in consensus)):
                row_color.append('#D5E8D4')   # green = match
            else:
                row_color.append('#FFE6CC')    # amber = mismatch

        cell_text.append(row_text)
        cell_colors.append(row_color)

    n_cols = len(col_labels)
    n_rows = len(genes)
    fig, ax = plt.subplots(
        figsize=(max(8, n_cols * 1.6), max(3, n_rows * 0.45 + 1.2)))
    ax.axis('off')

    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellColours=cell_colors,
        cellLoc='center',
        loc='center',
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)

    header_color = '#4472C4'
    for j in range(n_cols):
        tbl[(0, j)].set_facecolor(header_color)
        tbl[(0, j)].set_text_props(color='white', fontweight='bold')

    ax.set_title(f'HLA Typing Results — {sample_id}',
                 pad=14, fontsize=11, fontweight='bold')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_D1_tool_calls_table', fmt)
    plt.close(fig)


def fig_D2_summary_panel(output_dir, fmt):
    """3×2 multi-panel figure composed from the 6 key PNG outputs."""
    selected = [
        ('fig_A1_runtime',             'A1 — Execution Time'),
        ('fig_A2_peak_memory',         'A2 — Memory Usage'),
        ('fig_B1_confidence_per_gene', 'B1 — Confidence per Gene'),
        ('fig_B2_tool_concordance',    'B2 — Tool Concordance'),
        ('fig_C1_accuracy_heatmap',    'C1 — Tool Accuracy'),
        ('fig_C3_tool_ranking',        'C3 — Tool Ranking'),
    ]
    panels = [(title, Path(output_dir) / f'{name}.png')
              for name, title in selected]
    panels = [(title, p) for title, p in panels if p.exists()]

    if len(panels) < 2:
        return

    nrows, ncols = 3, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 4.5))
    fig.suptitle('HLA Typing Pipeline — Summary Figure',
                 fontsize=14, fontweight='bold', y=1.01)

    for ax, (title, p) in zip(axes.flatten(), panels):
        try:
            img = plt.imread(str(p))
            ax.imshow(img)
        except Exception:
            ax.text(0.5, 0.5, f'[{title}]', ha='center', va='center',
                    transform=ax.transAxes, fontsize=9)
        ax.axis('off')
        ax.set_title(title, fontsize=9, pad=4)

    for ax in axes.flatten()[len(panels):]:
        ax.axis('off')

    fig.tight_layout()
    _save(fig, output_dir, 'fig_D2_summary_panel', fmt)
    plt.close(fig)


# ─── Combine all PDFs ─────────────────────────────────────────────────────────

def combine_pdfs(output_dir):
    """Re-render all fig_*.png files into a single all_figures.pdf."""
    pngs = sorted(Path(output_dir).glob('fig_*.png'))
    if not pngs:
        return
    out_path = Path(output_dir) / 'all_figures.pdf'
    with PdfPages(out_path) as pp:
        for p in pngs:
            try:
                img = plt.imread(str(p))
                h, w = img.shape[:2]
                fig, ax = plt.subplots(
                    figsize=(w / FIGURE_DPI, h / FIGURE_DPI))
                ax.imshow(img)
                ax.axis('off')
                fig.subplots_adjust(0, 0, 1, 1)
                pp.savefig(fig, bbox_inches='tight', pad_inches=0)
                plt.close(fig)
            except Exception:
                pass
    print(f"  → {out_path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='Generate publication-quality figures from HLA typing pipeline outputs')
    p.add_argument('--consensus',    required=True)
    p.add_argument('--comparison',   required=True)
    p.add_argument('--trace_file',   required=True)
    p.add_argument('--weights_file', default=None)
    p.add_argument('--sample_id',    default='sample')
    p.add_argument('--output_dir',   default='scientific_figures')
    p.add_argument('--format',       default='both',
                   choices=['pdf', 'png', 'both'])
    args = p.parse_args()

    setup_style()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f'[make_paper_figures] sample={args.sample_id}')

    trace_df = pd.DataFrame()
    try:
        trace_df = parse_trace(args.trace_file)
        print(f'  trace: {len(trace_df)} tool records')
    except Exception as e:
        print(f'  WARNING trace: {e}')

    consensus_df = pd.DataFrame()
    try:
        consensus_df = parse_consensus(args.consensus)
        print(f'  consensus: {len(consensus_df)} genes')
    except Exception as e:
        print(f'  WARNING consensus: {e}')

    gene_data, tool_names = {}, []
    try:
        gene_data, tool_names = parse_comparison(args.comparison)
        print(f'  comparison: {len(gene_data)} genes, tools={tool_names}')
    except Exception as e:
        print(f'  WARNING comparison: {e}')

    weights_data = {}
    if args.weights_file and Path(args.weights_file).exists():
        try:
            weights_data = load_weights(args.weights_file)
            n = weights_data.get('n_samples', '?')
            print(f'  weights: N={n}')
        except Exception as e:
            print(f'  WARNING weights: {e}')

    od  = args.output_dir
    fmt = args.format

    # ── Group A ───────────────────────────────────────────────────────────────
    if not trace_df.empty:
        print('  Group A: resource benchmarking ...')
        for fn in [fig_A1_runtime, fig_A2_peak_memory, fig_A3_cpu_utilization,
                   fig_A4_io_throughput, fig_A5_gantt, fig_A6_resource_efficiency]:
            try:
                fn(trace_df, od, fmt)
                print(f'    ✓ {fn.__name__}')
            except Exception as e:
                print(f'    ✗ {fn.__name__}: {e}')

    # ── Group B ───────────────────────────────────────────────────────────────
    if not consensus_df.empty:
        print('  Group B: HLA typing quality ...')
        for fn, kw in [
            (fig_B1_confidence_per_gene,
             {'consensus_df': consensus_df, 'sample_id': args.sample_id}),
            (fig_B3_read_support,
             {'consensus_df': consensus_df, 'sample_id': args.sample_id}),
            (fig_B4_confidence_vs_reads,
             {'consensus_df': consensus_df, 'sample_id': args.sample_id}),
        ]:
            try:
                fn(**kw, output_dir=od, fmt=fmt)
                print(f'    ✓ {fn.__name__}')
            except Exception as e:
                print(f'    ✗ {fn.__name__}: {e}')

        if gene_data:
            try:
                fig_B2_tool_concordance(
                    gene_data, tool_names, args.sample_id, od, fmt)
                print('    ✓ fig_B2_tool_concordance')
            except Exception as e:
                print(f'    ✗ fig_B2_tool_concordance: {e}')

    # ── Group C ───────────────────────────────────────────────────────────────
    if weights_data:
        print('  Group C: tool accuracy ...')
        for fn in [fig_C1_accuracy_heatmap, fig_C2_weight_heatmap,
                   fig_C3_tool_ranking]:
            try:
                fn(weights_data, od, fmt)
                print(f'    ✓ {fn.__name__}')
            except Exception as e:
                print(f'    ✗ {fn.__name__}: {e}')

        if not trace_df.empty:
            try:
                fig_C4_accuracy_vs_runtime(weights_data, trace_df, od, fmt)
                print('    ✓ fig_C4_accuracy_vs_runtime')
            except Exception as e:
                print(f'    ✗ fig_C4_accuracy_vs_runtime: {e}')

    # ── Group D ───────────────────────────────────────────────────────────────
    print('  Group D: summary panels ...')
    if gene_data:
        try:
            fig_D1_tool_calls_table(
                gene_data, tool_names, consensus_df, args.sample_id, od, fmt)
            print('    ✓ fig_D1_tool_calls_table')
        except Exception as e:
            print(f'    ✗ fig_D1_tool_calls_table: {e}')

    try:
        fig_D2_summary_panel(od, fmt)
        print('    ✓ fig_D2_summary_panel')
    except Exception as e:
        print(f'    ✗ fig_D2_summary_panel: {e}')

    if fmt in ('pdf', 'both'):
        try:
            combine_pdfs(od)
            print('    ✓ all_figures.pdf')
        except Exception as e:
            print(f'    ✗ combine_pdfs: {e}')

    pdfs = len(list(Path(od).glob('fig_*.pdf')))
    pngs = len(list(Path(od).glob('fig_*.png')))
    print(f'\n[make_paper_figures] Done — {pdfs} PDFs, {pngs} PNGs → {od}')


if __name__ == '__main__':
    main()
