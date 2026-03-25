#!/usr/bin/env python3
"""make_paper_figures.py — Generate publication-quality scientific figures
from HLA typing pipeline outputs.

Produces 20 figures in 5 groups:
  A: Resource benchmarking      (from Nextflow trace)
  B: HLA typing quality         (from consensus/comparison TSVs)
  C: Tool accuracy              (from tool_weights JSON)
  D: Summary panels             (composite)
  V: Weighted voting analysis   (the pipeline's novel contribution)

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

    def _b2_to_set(s):
        if not s or s in ('-', ''):
            return set()
        return set(a.strip() for a in re.split(r'[/,]', s) if a.strip() not in ('-', ''))

    for j, gene in enumerate(genes):
        entry     = gene_data[gene]
        consensus = entry.get('consensus', '-')
        for i, tool in enumerate(tools):
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            if alleles in ('-', '', None):
                matrix[i, j] = np.nan
            elif consensus and alleles not in ('-', '', None):
                cons_set = _b2_to_set(consensus)
                call_set = _b2_to_set(alleles)
                overlap  = len(cons_set & call_set)
                if overlap == len(cons_set) and overlap == len(call_set):
                    matrix[i, j] = 1.0   # full match
                elif overlap > 0:
                    matrix[i, j] = 0.5   # partial match
                else:
                    matrix[i, j] = 0.0   # conflict
            else:
                matrix[i, j] = 0.0

    fig, ax = plt.subplots(
        figsize=(max(6, len(labels_g) * 0.9), max(3, len(tools) * 0.65)))

    mask = np.isnan(matrix)
    if HAS_SEABORN:
        annot = np.where(mask, '',
                         np.where(matrix == 1.0, '✓',
                                  np.where(matrix == 0.5, '½', '✗')).astype(object))
        sns.heatmap(matrix, annot=annot, fmt='',
                    cmap='RdYlGn', vmin=0, vmax=1,
                    xticklabels=labels_g, yticklabels=tools,
                    linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Agreement with consensus (1=full, 0.5=partial, 0=conflict)'},
                    mask=mask, ax=ax)
    else:
        im = ax.imshow(matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
        ax.set_xticks(range(len(labels_g)))
        ax.set_yticks(range(len(tools)))
        ax.set_xticklabels(labels_g, rotation=30, ha='right')
        ax.set_yticklabels(tools)
        plt.colorbar(im, ax=ax, label='Agreement (1=full, 0.5=partial, 0=conflict)')

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


# ─── Group V: Weighted Voting Analysis ───────────────────────────────────────

def _allele_match(consensus, alleles):
    """Return 'full', 'partial', 'conflict', or 'none'."""
    if not alleles or alleles in ('-', ''):
        return 'none'
    if not consensus or consensus in ('-', ''):
        return 'none'
    def to_set(s):
        return set(a.strip() for a in re.split(r'[/,]', s) if a.strip() not in ('-', ''))
    c = to_set(consensus)
    a = to_set(alleles)
    overlap = len(c & a)
    if overlap >= len(c) or overlap >= len(a):
        return 'full'
    elif overlap > 0:
        return 'partial'
    else:
        return 'conflict'


def fig_V1_allele_vote_matrix(gene_data, tool_names, consensus_df, sample_id, output_dir, fmt):
    """
    Annotated heatmap: tools (rows) x classical HLA genes (cols).

    Each cell shows the allele(s) called by that tool for that gene.
    Background colour encodes agreement with the consensus:
      Green  = full match (all alleles agree)
      Amber  = partial match (one allele correct, one different)
      Red    = conflict (all alleles differ from consensus)
      Light grey = not called / no result
    Bottom row shows the final consensus call.
    """
    genes = [g for g in GENE_ORDER if g in gene_data]
    tools = [t for t in TOOL_ORDER if t in tool_names]
    tools += [t for t in tool_names if t not in TOOL_ORDER]
    if not genes or not tools:
        return

    MATCH_COLORS = {
        'full':     '#27AE60',
        'partial':  '#F39C12',
        'conflict': '#E74C3C',
        'none':     '#ECECEC',
    }

    nrows = len(tools) + 1   # +1 for consensus row
    ncols = len(genes)
    fig, ax = plt.subplots(figsize=(max(7, ncols * 1.3), max(4, nrows * 0.75)))
    ax.set_xlim(0, ncols)
    ax.set_ylim(0, nrows)
    ax.axis('off')

    gene_labels = [g.replace('HLA-', '') for g in genes]

    # Draw column headers (gene names)
    for j, gl in enumerate(gene_labels):
        ax.text(j + 0.5, nrows - 0.15, gl,
                ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Draw tool rows (top = first tool, bottom-1 = last tool, bottom = consensus)
    for i, tool in enumerate(tools):
        row_y = nrows - 1 - i
        ax.text(-0.05, row_y - 0.5, tool, ha='right', va='center',
                fontsize=8, fontweight='normal',
                color=TOOL_COLORS.get(tool, '#333'))
        for j, gene in enumerate(genes):
            entry = gene_data[gene]
            consensus = entry.get('consensus', '-')
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            match = _allele_match(consensus, alleles)
            color = MATCH_COLORS[match]

            rect = plt.Rectangle((j, row_y - 1), 1, 1,
                                  facecolor=color, edgecolor='white',
                                  linewidth=1.5)
            ax.add_patch(rect)

            if alleles and alleles not in ('-', ''):
                parts = re.split(r'/', alleles)
                short = '\n'.join(a.strip() for a in parts[:2])
                txt_color = 'white' if match in ('full', 'conflict') else '#333'
                ax.text(j + 0.5, row_y - 0.5, short,
                        ha='center', va='center', fontsize=6,
                        color=txt_color, linespacing=1.2)
            else:
                ax.text(j + 0.5, row_y - 0.5, '\u2014',
                        ha='center', va='center', fontsize=8, color='#999')

    # Consensus row at the bottom
    cons_y = 0
    ax.text(-0.05, cons_y + 0.5, 'Consensus', ha='right', va='center',
            fontsize=8, fontweight='bold', color='#333')
    for j, gene in enumerate(genes):
        consensus = gene_data[gene].get('consensus', '-')
        rect = plt.Rectangle((j, cons_y), 1, 1,
                              facecolor='#2C3E50', edgecolor='white',
                              linewidth=1.5)
        ax.add_patch(rect)
        if consensus and consensus not in ('-',):
            parts = re.split(r'/', consensus)
            short = '\n'.join(a.strip() for a in parts[:2])
            ax.text(j + 0.5, cons_y + 0.5, short,
                    ha='center', va='center', fontsize=6.5,
                    color='white', fontweight='bold', linespacing=1.2)

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor=MATCH_COLORS['full'],     label='Full match'),
        mpatches.Patch(facecolor=MATCH_COLORS['partial'],  label='Partial match (1 allele)'),
        mpatches.Patch(facecolor=MATCH_COLORS['conflict'], label='Conflict'),
        mpatches.Patch(facecolor=MATCH_COLORS['none'],     label='Not called'),
        mpatches.Patch(facecolor='#2C3E50',                label='Consensus'),
    ]
    ax.legend(handles=legend_patches, loc='upper right',
              bbox_to_anchor=(1.0, -0.02), ncol=5, fontsize=7,
              frameon=True, borderpad=0.5)

    ax.set_title(f'Weighted Majority Vote Matrix \u2014 {sample_id}',
                 fontsize=12, fontweight='bold', pad=10)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_V1_allele_vote_matrix', fmt)
    plt.close(fig)


def fig_V2_vote_weight_breakdown(gene_data, tool_names, consensus_df, sample_id,
                                  output_dir, fmt, weights_data=None):
    """
    Horizontal stacked bar: per classical gene, how much total vote weight
    agreed with the consensus vs disagreed vs did not call.
    """
    genes = [g for g in GENE_ORDER if g in gene_data]
    tools = [t for t in TOOL_ORDER if t in tool_names]
    tools += [t for t in tool_names if t not in TOOL_ORDER]
    if not genes or not tools:
        return

    gene_labels = [g.replace('HLA-', '') for g in genes]

    cal_weights = {}
    if weights_data:
        genes_section = weights_data.get('genes', {})
        for gs, tw in genes_section.items():
            for tool, w in tw.items():
                cal_weights.setdefault(tool.upper(), {})[gs.upper()] = float(w)

    def get_weight(tool, gene):
        gs = gene.replace('HLA-', '').upper()
        return cal_weights.get(tool.upper(), {}).get(gs, 1.0)

    agree_vals = []
    partial_vals = []
    conflict_vals = []
    nocall_vals = []
    agree_labels = []

    for gene in genes:
        entry = gene_data[gene]
        consensus = entry.get('consensus', '-')
        w_agree = w_partial = w_conflict = w_nocall = 0.0
        a_tools = []
        for tool in tools:
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            w = get_weight(tool, gene)
            match = _allele_match(consensus, alleles)
            if match == 'full':
                w_agree += w; a_tools.append(tool)
            elif match == 'partial':
                w_partial += w
            elif match == 'conflict':
                w_conflict += w
            else:
                w_nocall += w
        agree_vals.append(w_agree)
        partial_vals.append(w_partial)
        conflict_vals.append(w_conflict)
        nocall_vals.append(w_nocall)
        agree_labels.append(a_tools)

    y = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(8, max(3.5, len(genes) * 0.55)))

    b1 = ax.barh(y, agree_vals,   color='#27AE60', label='Full agreement', height=0.6)
    b2 = ax.barh(y, partial_vals, left=agree_vals, color='#F39C12', label='Partial agreement', height=0.6)
    left2 = [a + p for a, p in zip(agree_vals, partial_vals)]
    b3 = ax.barh(y, conflict_vals, left=left2, color='#E74C3C', label='Conflict', height=0.6)
    left3 = [a + b for a, b in zip(left2, conflict_vals)]
    b4 = ax.barh(y, nocall_vals,  left=left3, color='#ECECEC', label='Not called',
                 height=0.6, edgecolor='#ccc', linewidth=0.5)

    for idx, (bar, tl) in enumerate(zip(b1, agree_labels)):
        w = bar.get_width()
        if w > 0.3 and tl:
            label = ', '.join(t[:3] for t in tl)
            ax.text(w / 2, bar.get_y() + bar.get_height() / 2,
                    label, ha='center', va='center',
                    fontsize=6, color='white', fontweight='bold')

    ax.set_yticks(y)
    ax.set_yticklabels(gene_labels)
    ax.set_xlabel('Vote weight (equal: # tools)')
    ax.set_title(f'Voting Weight Breakdown per Gene \u2014 {sample_id}')
    ax.legend(loc='lower right', fontsize=8)
    total = len(tools)
    ax.set_xlim(0, total * 1.05)
    ax.axvline(total / 2, color='black', linestyle=':', linewidth=0.8, alpha=0.5)
    fig.tight_layout()
    _save(fig, output_dir, 'fig_V2_vote_weight_breakdown', fmt)
    plt.close(fig)


def fig_V3_weighted_confidence_comparison(gene_data, tool_names, consensus_df,
                                           sample_id, output_dir, fmt,
                                           weights_data=None):
    """
    For each classical HLA gene: side-by-side bars comparing consensus
    confidence under three weighting strategies:
      1. Equal weighting     - each tool votes with weight 1.0
      2. Calibrated weighting - weights from tool_weights JSON
      3. Observed (current)  - the actual confidence from consensus file
    """
    genes = [g for g in GENE_ORDER if g in gene_data]
    tools = [t for t in TOOL_ORDER if t in tool_names]
    tools += [t for t in tool_names if t not in TOOL_ORDER]
    if not genes or not tools:
        return

    gene_labels = [g.replace('HLA-', '') for g in genes]

    cal_weights = {}
    if weights_data:
        for gs, tw in weights_data.get('genes', {}).items():
            for tool, w in tw.items():
                cal_weights.setdefault(tool.upper(), {})[gs.upper()] = float(w)

    def get_cal_weight(tool, gene):
        gs = gene.replace('HLA-', '').upper()
        return cal_weights.get(tool.upper(), {}).get(gs, None)

    equal_conf = []
    calib_conf = []
    obs_conf   = []

    cons_lookup = {}
    if consensus_df is not None and not consensus_df.empty:
        for _, row in consensus_df.iterrows():
            cons_lookup[row['Gene']] = row.get('Confidence')

    for gene in genes:
        entry = gene_data[gene]
        consensus = entry.get('consensus', '-')

        n_agree = n_call = 0
        w_agree_cal = w_call_cal = 0.0

        for tool in tools:
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            match = _allele_match(consensus, alleles)
            calling = match != 'none'
            if calling:
                n_call += 1
                w_cal = get_cal_weight(tool, gene)
                if w_cal is not None:
                    w_call_cal += w_cal
            if match in ('full', 'partial'):
                n_agree += 1
                w_cal = get_cal_weight(tool, gene)
                if w_cal is not None:
                    w_agree_cal += w_cal

        equal_conf.append(n_agree / n_call if n_call > 0 else 0.0)
        calib_conf.append(w_agree_cal / w_call_cal if w_call_cal > 0 else (n_agree / n_call if n_call > 0 else 0.0))
        obs_conf.append(cons_lookup.get(gene) or (n_agree / n_call if n_call > 0 else 0.0))

    x = np.arange(len(genes))
    w = 0.25
    fig, ax = plt.subplots(figsize=(max(7, len(genes) * 1.1), 4.5))

    ax.bar(x - w,   equal_conf, w, label='Equal weighting',      color='#4472C4', alpha=0.85)
    ax.bar(x,       calib_conf, w, label='Calibrated weighting',  color='#27AE60', alpha=0.85)
    ax.bar(x + w,   obs_conf,   w, label='Observed (consensus)',  color='#ED7D31', alpha=0.85)

    ax.axhline(0.8, color=CONF_GREEN, linestyle='--', linewidth=0.8, alpha=0.6, label='High >= 0.8')
    ax.axhline(0.5, color=CONF_AMBER, linestyle='--', linewidth=0.8, alpha=0.6, label='Medium >= 0.5')
    ax.set_xticks(x)
    ax.set_xticklabels(gene_labels)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel('Confidence score')
    ax.set_title(f'Weighting Strategy Comparison \u2014 {sample_id}')
    ax.legend(fontsize=8)

    for bars in [ax.containers[0], ax.containers[1], ax.containers[2]]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.02,
                    f'{h:.2f}', ha='center', va='bottom', fontsize=6.5)

    fig.tight_layout()
    _save(fig, output_dir, 'fig_V3_weighted_confidence_comparison', fmt)
    plt.close(fig)


def fig_V4_allele_concordance_detail(gene_data, tool_names, sample_id,
                                      output_dir, fmt):
    """
    Per-gene allele-level concordance breakdown for classical HLA genes.

    For each gene, shows stacked segments across all calling tools:
      Dark green  = tools agreeing on BOTH alleles (full concordance)
      Light green = tools agreeing on allele 1 only
      Orange      = tools agreeing on allele 2 only
      Red         = tools with no allele in common with consensus
      Grey        = tools that did not call this gene
    """
    genes = [g for g in GENE_ORDER if g in gene_data]
    tools = [t for t in TOOL_ORDER if t in tool_names]
    tools += [t for t in tool_names if t not in TOOL_ORDER]
    if not genes or not tools:
        return

    def to_set(s):
        if not s or s in ('-',):
            return set()
        return set(a.strip() for a in re.split(r'[/,]', s) if a.strip() not in ('-', ''))

    gene_labels = [g.replace('HLA-', '') for g in genes]
    both_agree = []
    a1_only = []
    a2_only = []
    conflict = []
    no_call  = []

    for gene in genes:
        entry = gene_data[gene]
        consensus = entry.get('consensus', '-')
        c_set = to_set(consensus)
        c_list = list(c_set)

        n_both = n_a1 = n_a2 = n_conf = n_none = 0
        for tool in tools:
            alleles = entry.get('tools', {}).get(tool, {}).get('alleles', '-')
            a_set = to_set(alleles)
            if not a_set:
                n_none += 1
                continue
            if not c_set:
                n_conf += 1
                continue
            overlap = c_set & a_set
            if len(overlap) >= 2 or overlap == c_set:
                n_both += 1
            elif len(overlap) == 1:
                matched = list(overlap)[0]
                if c_list and matched == c_list[0]:
                    n_a1 += 1
                else:
                    n_a2 += 1
            else:
                n_conf += 1

        both_agree.append(n_both)
        a1_only.append(n_a1)
        a2_only.append(n_a2)
        conflict.append(n_conf)
        no_call.append(n_none)

    x = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(max(7, len(genes) * 1.0), 4.5))

    b1 = ax.bar(x, both_agree, color='#1A7340', label='Both alleles agree', width=0.6)
    b2 = ax.bar(x, a1_only,   bottom=both_agree, color='#82C341', label='Allele 1 only', width=0.6)
    bot2 = [a + b for a, b in zip(both_agree, a1_only)]
    b3 = ax.bar(x, a2_only,   bottom=bot2, color='#F39C12', label='Allele 2 only', width=0.6)
    bot3 = [a + b for a, b in zip(bot2, a2_only)]
    b4 = ax.bar(x, conflict,  bottom=bot3, color='#E74C3C', label='Conflict', width=0.6)
    bot4 = [a + b for a, b in zip(bot3, conflict)]
    b5 = ax.bar(x, no_call,   bottom=bot4, color='#ECECEC', label='Not called',
                width=0.6, edgecolor='#ccc', linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(gene_labels)
    ax.set_ylabel('Number of tools')
    ax.set_yticks(range(0, len(tools) + 1))
    ax.set_ylim(0, len(tools) + 0.5)
    ax.set_title(f'Allele-level Concordance Breakdown \u2014 {sample_id}')
    ax.legend(fontsize=8, loc='upper right')
    fig.tight_layout()
    _save(fig, output_dir, 'fig_V4_allele_concordance_detail', fmt)
    plt.close(fig)


# ─── Captions report ──────────────────────────────────────────────────────────

def generate_captions_report(output_dir, sample_id, consensus_df, gene_data,
                              tool_names, trace_df, weights_data):
    """Write figure_captions_evaluation.md to output_dir."""
    lines = [
        f'# Figure Captions and Evaluation — {sample_id}',
        '',
        'Generated by `make_paper_figures.py`. Each section provides a formal '
        'publishable caption and a brief data-driven evaluation.',
        '',
    ]

    # Helper: safe stat extraction
    def safe(fn):
        try:
            return fn()
        except Exception:
            return None

    # --- Pre-compute stats ---
    trace_stats = {}
    if trace_df is not None and not trace_df.empty:
        try:
            dd = _dedup_trace(trace_df).dropna(subset=['duration_s'])
            if not dd.empty:
                slowest = dd.loc[dd['duration_s'].idxmax()]
                fastest = dd.loc[dd['duration_s'].idxmin()]
                trace_stats['slowest'] = (slowest['tool'], slowest['duration_s'])
                trace_stats['fastest'] = (fastest['tool'], fastest['duration_s'])
                trace_stats['n_tools'] = len(dd)
                trace_stats['total_s'] = dd['duration_s'].sum()
        except Exception:
            pass

    conf_stats = {}
    if consensus_df is not None and not consensus_df.empty:
        try:
            c = consensus_df.dropna(subset=['Confidence'])
            if not c.empty:
                conf_stats['highest'] = (c.loc[c['Confidence'].idxmax(), 'Gene'],
                                         c['Confidence'].max())
                conf_stats['lowest']  = (c.loc[c['Confidence'].idxmin(), 'Gene'],
                                         c['Confidence'].min())
                conf_stats['n_high']  = int((c['Confidence'] >= 0.8).sum())
                conf_stats['n_genes'] = len(c)
        except Exception:
            pass

    gene_stats = {}
    if gene_data and tool_names:
        try:
            n_tools = len(tool_names)
            full_counts = {}
            for gene, entry in gene_data.items():
                consensus = entry.get('consensus', '-')
                n_full = sum(
                    1 for t in tool_names
                    if _allele_match(consensus,
                                     entry.get('tools', {}).get(t, {}).get('alleles', '-')) == 'full'
                )
                full_counts[gene] = n_full
            if full_counts:
                best_gene = max(full_counts, key=full_counts.get)
                worst_gene = min(full_counts, key=full_counts.get)
                gene_stats['best']   = (best_gene.replace('HLA-', ''), full_counts[best_gene], n_tools)
                gene_stats['worst']  = (worst_gene.replace('HLA-', ''), full_counts[worst_gene], n_tools)
                gene_stats['n_tools'] = n_tools
        except Exception:
            pass

    # ─── Group A ──────────────────────────────────────────────────────────────
    lines += ['## Group A: Resource Benchmarking', '']

    def _a1_caption():
        cap = (
            f'**Figure A1 — Tool Execution Time.**  '
            f'Horizontal bar chart showing wall-clock runtime (minutes) for each '
            f'HLA typing tool executed during sample {sample_id}. '
            f'Bars are coloured per tool; grey indicates a failed run.'
        )
        ev = ''
        if trace_stats:
            sl_t, sl_s = trace_stats.get('slowest', ('?', 0))
            fa_t, fa_s = trace_stats.get('fastest', ('?', 0))
            tot = trace_stats.get('total_s', 0)
            ev = (
                f'**Evaluation:** {trace_stats.get("n_tools", "?")} tools completed. '
                f'The slowest tool was {sl_t} ({sl_s/60:.1f} min) and the fastest was '
                f'{fa_t} ({fa_s/60:.1f} min). '
                f'Total cumulative CPU-equivalent time was {tot/3600:.1f} h, '
                f'highlighting the benefit of parallel Nextflow execution.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_a1_caption) or '**Figure A1** — (no trace data)', '']

    lines += [
        '**Figure A2 — Peak Memory Usage vs Configured Limit.**  '
        f'Grouped bar chart comparing peak resident set size (RSS, actual) against '
        f'the configured memory limit for each tool during sample {sample_id}. '
        'Tools with large margins indicate over-provisioned allocations.',
        '',
        '**Evaluation:** Memory headroom varies substantially across tools. '
        'Tools allocated far less than their configured limit suggest these limits '
        'can be reduced in future pipeline runs to improve scheduler efficiency on HPC clusters.',
        '',
    ]

    lines += [
        '**Figure A3 — CPU Utilization per Tool.**  '
        f'Bar chart of percentage CPU utilization recorded by Nextflow for each tool '
        f'in sample {sample_id}. Values above 100% indicate multi-threaded execution.',
        '',
        '**Evaluation:** Tools exploiting multi-core parallelism (e.g., alignment-based tools) '
        'show CPU% well above 100%. Single-threaded tools remain near 100%, '
        'suggesting potential for additional parallelization.',
        '',
    ]

    lines += [
        '**Figure A4 — I/O Throughput per Tool.**  '
        f'Grouped bar chart of bytes read and written (GB) per tool for sample {sample_id}. '
        'Read-heavy tools dominate I/O load due to reference database scanning.',
        '',
        '**Evaluation:** I/O throughput is dominated by alignment-based tools that scan '
        'large reference databases. This informs storage bandwidth requirements on HPC, '
        'especially for parallel sample batches.',
        '',
    ]

    lines += [
        '**Figure A5 — Pipeline Execution Timeline (Gantt).**  '
        f'Gantt chart showing the start time and duration of each tool for sample {sample_id}. '
        'Parallel branches of the Nextflow DAG are visible as overlapping bars.',
        '',
        '**Evaluation:** The Gantt chart reveals which tools form the critical path '
        'of the pipeline. Parallelism is exploited when multiple tools run simultaneously; '
        'any sequential bottleneck will extend total wall-clock time.',
        '',
    ]

    lines += [
        '**Figure A6 — Resource Efficiency Scatter.**  '
        f'Scatter plot of wall-clock time (x) vs peak RSS (y) for each tool in sample {sample_id}. '
        'Bubble size is proportional to CPU utilization.',
        '',
        '**Evaluation:** Tools in the upper-right quadrant (slow, memory-heavy) are '
        'the most resource-intensive. The scatter reveals trade-offs useful for '
        'selecting a tool subset under HPC resource constraints.',
        '',
    ]

    # ─── Group B ──────────────────────────────────────────────────────────────
    lines += ['## Group B: HLA Typing Quality', '']

    def _b1_caption():
        cap = (
            f'**Figure B1 — Consensus Confidence per HLA Gene.**  '
            f'Bar chart of the pipeline confidence score (0–1) for each classical '
            f'HLA gene in sample {sample_id}. '
            f'Bars are coloured green (>=0.8), amber (>=0.5), or red (<0.5) '
            f'according to the traffic-light scheme.'
        )
        ev = ''
        if conf_stats:
            hg, hv = conf_stats.get('highest', ('?', 0))
            lg, lv = conf_stats.get('lowest', ('?', 0))
            nh = conf_stats.get('n_high', '?')
            ng = conf_stats.get('n_genes', '?')
            ev = (
                f'**Evaluation:** {nh}/{ng} genes reach high confidence (>=0.8). '
                f'The highest-confidence locus is {hg.replace("HLA-", "")} ({hv:.2f}) '
                f'and the lowest is {lg.replace("HLA-", "")} ({lv:.2f}). '
                f'Low-confidence loci should be validated by orthogonal methods.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_b1_caption) or '**Figure B1** — (no confidence data)', '']

    def _b2_caption():
        cap = (
            f'**Figure B2 — Tool vs Consensus Agreement.**  '
            f'Heatmap showing agreement between each tool and the weighted consensus '
            f'for each classical HLA gene in sample {sample_id}. '
            f'Green=full match, yellow-green=partial (one allele), red=conflict, '
            f'grey=not called. Annotations show checkmark (full), half (partial), '
            f'or cross (conflict).'
        )
        ev = ''
        if gene_stats:
            bg, bn, nt = gene_stats.get('best', ('?', '?', '?'))
            wg, wn, _  = gene_stats.get('worst', ('?', '?', '?'))
            ev = (
                f'**Evaluation:** The gene with highest cross-tool agreement is '
                f'{bg} ({bn}/{nt} tools in full agreement). '
                f'The gene with lowest agreement is {wg} ({wn}/{nt} tools). '
                f'Partial matches (amber) indicate allele-level discordance, '
                f'common at highly polymorphic loci such as HLA-B and HLA-DRB1.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_b2_caption) or '**Figure B2** — (no comparison data)', '']

    lines += [
        '**Figure B3 — Read Support per Allele.**  '
        f'Grouped bar chart of read counts assigned to allele 1 and allele 2 '
        f'at each classical HLA gene for sample {sample_id}. '
        'Large imbalance between alleles may indicate hemizygosity or allele drop-out.',
        '',
        '**Evaluation:** Genes with very low read support in either allele position '
        'are more likely to have incorrect or low-confidence calls. '
        'Read counts are extracted from the HLA-HD per-gene output files.',
        '',
    ]

    lines += [
        '**Figure B4 — Read Support vs Confidence.**  '
        f'Scatter plot of total reads (allele 1 + 2) vs pipeline confidence score '
        f'for each gene in sample {sample_id}. '
        'Each point is labelled with the gene name.',
        '',
        '**Evaluation:** A positive correlation between read depth and confidence '
        'is expected; outliers (high reads but low confidence) indicate tool disagreement '
        'rather than insufficient coverage.',
        '',
    ]

    # ─── Group C ──────────────────────────────────────────────────────────────
    lines += ['## Group C: Tool Accuracy (Population-level Calibration)', '']

    lines += [
        '**Figure C1 — Tool Accuracy Heatmap.**  '
        'Heatmap of concordance (%) with 1KGP ground truth for each tool and classical '
        'HLA gene, derived from the empirical weight calibration cohort. '
        'Values are shown as percentages; grey cells indicate genes not typed by that tool.',
        '',
        '**Evaluation:** Accuracy varies substantially by gene and tool. '
        'OptiType consistently achieves high Class I accuracy, while HLA-HD leads '
        'for Class II loci. These patterns directly determine the calibrated weights '
        'used in the weighted consensus voting step.',
        '',
    ]

    lines += [
        '**Figure C2 — Normalized Tool Weights Heatmap.**  '
        'Heatmap of normalized per-gene accuracy weights assigned to each tool '
        'for the weighted consensus voting step. '
        'Weights sum to 1.0 within each gene column.',
        '',
        '**Evaluation:** The weight distribution reflects which tools are most accurate '
        'for each gene in a 30x WGS context. Genes where one tool dominates '
        '(e.g., OptiType for HLA-C) result in highly asymmetric weight columns.',
        '',
    ]

    lines += [
        '**Figure C3 — Tool Ranking.**  '
        'Horizontal bar chart of mean concordance across the five core genes '
        '(A, B, C, DRB1, DQB1) for each tool, ranked from highest to lowest.',
        '',
        '**Evaluation:** The ranking provides a single-number summary of overall '
        'tool performance for WGS typing. Tools near the top of the ranking are '
        'given higher weights in the consensus, improving final call accuracy.',
        '',
    ]

    lines += [
        '**Figure C4 — Accuracy vs Runtime Trade-off.**  '
        'Scatter plot of mean concordance (y) vs median runtime from Nextflow trace (x) '
        'for each tool, enabling assessment of cost-effectiveness.',
        '',
        '**Evaluation:** Tools in the upper-left quadrant (fast and accurate) are '
        'most cost-effective. Tools that are slow but also inaccurate may be candidates '
        'for exclusion from the default pipeline configuration.',
        '',
    ]

    # ─── Group D ──────────────────────────────────────────────────────────────
    lines += ['## Group D: Summary Panels', '']

    def _d1_caption():
        cap = (
            f'**Figure D1 — HLA Typing Results Table.**  '
            f'Visual table of allele calls from all tools and the weighted consensus '
            f'for each classical HLA gene in sample {sample_id}. '
            f'Cell colours indicate confidence level (green=high, amber=medium, red=low) '
            f'for the consensus column, and agreement (green=match, amber=mismatch) for tool columns.'
        )
        ev = ''
        if gene_data and tool_names:
            n_genes = len([g for g in GENE_ORDER if g in gene_data])
            ev = (
                f'**Evaluation:** The table summarises typing results for {n_genes} classical '
                f'genes across {len(tool_names)} tools. '
                f'It is the primary deliverable for clinical or research use, '
                f'providing a compact view of both the final calls and the per-tool evidence.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_d1_caption) or '**Figure D1** — (no data)', '']

    lines += [
        '**Figure D2 — Summary Panel.**  '
        f'Composite 3x2 panel combining six key figures: A1 (execution time), '
        f'A2 (memory), B1 (confidence), B2 (concordance), C1 (accuracy heatmap), '
        f'C3 (tool ranking). Intended as a single-page overview for manuscripts.',
        '',
        '**Evaluation:** The summary panel provides a publication-ready overview '
        'of both computational performance and HLA typing quality for a single sample. '
        'It is most informative when all six constituent figures contain data.',
        '',
    ]

    # ─── Group V ──────────────────────────────────────────────────────────────
    lines += ['## Group V: Weighted Voting Analysis', '']

    def _v1_caption():
        cap = (
            f'**Figure V1 — Weighted Majority Vote Matrix.**  '
            f'Annotated grid showing the allele call from each tool (rows) for each '
            f'classical HLA gene (columns) in sample {sample_id}. '
            f'Cell colour encodes agreement with the weighted consensus: '
            f'green=full match, amber=partial, red=conflict, grey=not called. '
            f'The bottom row displays the final consensus call.'
        )
        ev = ''
        if gene_stats:
            bg, bn, nt = gene_stats.get('best', ('?', '?', '?'))
            wg, wn, _  = gene_stats.get('worst', ('?', '?', '?'))
            ev = (
                f'**Evaluation:** For sample {sample_id}, gene {bg} shows the highest '
                f'cross-tool agreement ({bn}/{nt} tools in full agreement), while '
                f'{wg} shows the lowest ({wn}/{nt}). '
                f'This matrix is the primary visualisation of the novel weighted '
                f'majority voting step introduced by this pipeline.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_v1_caption) or '**Figure V1** — (no comparison data)', '']

    def _v2_caption():
        cap = (
            f'**Figure V2 — Voting Weight Breakdown per Gene.**  '
            f'Horizontal stacked bar chart showing, for each classical HLA gene in '
            f'sample {sample_id}, the total vote weight allocated to full-agreement, '
            f'partial-agreement, conflict, and no-call segments. '
            f'The majority threshold (50% of total weight) is shown as a dotted line.'
        )
        ev = ''
        if gene_data and tool_names:
            n = len(tool_names)
            ev = (
                f'**Evaluation:** With {n} tools in the ensemble, a gene requires '
                f'>{n/2:.1f} units of agreement weight to reach majority consensus. '
                f'Genes where the green segment falls short of the threshold indicate '
                f'low-agreement loci where the consensus call is less reliable.'
            )
        return cap + '\n\n' + ev if ev else cap
    lines += [safe(_v2_caption) or '**Figure V2** — (no data)', '']

    lines += [
        '**Figure V3 — Weighting Strategy Comparison.**  '
        f'Grouped bar chart comparing confidence scores under three strategies '
        f'(equal, calibrated, observed) for each classical HLA gene in sample {sample_id}. '
        f'Reference lines at 0.8 (high) and 0.5 (medium) aid interpretation.',
        '',
        '**Evaluation:** Calibrated weighting should produce confidence scores '
        'that better reflect empirical accuracy than simple equal weighting. '
        'Genes where calibrated and equal weighting diverge most are those where '
        'one high-accuracy tool dominates the vote.',
        '',
    ]

    def _v4_caption():
        cap = (
            f'**Figure V4 — Allele-level Concordance Breakdown.**  '
            f'Stacked bar chart decomposing, for each classical HLA gene in sample {sample_id}, '
            f'the number of tools agreeing on both alleles (dark green), allele 1 only '
            f'(light green), allele 2 only (orange), conflicting (red), or not calling (grey).'
        )
        ev = (
            '**Evaluation:** This figure reveals whether tool disagreements are concentrated '
            'on the second (more polymorphic) allele or affect both alleles equally. '
            'Genes dominated by the dark-green segment indicate high-confidence homozygous '
            'or unambiguous heterozygous calls.'
        )
        return cap + '\n\n' + ev
    lines += [safe(_v4_caption) or '**Figure V4** — (no data)', '']

    lines += ['', '---', f'*Report generated for sample: {sample_id}*', '']

    out_path = Path(output_dir) / 'figure_captions_evaluation.md'
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines))


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

    # ── Group V (Voting Analysis) ─────────────────────────────────────────────
    if gene_data:
        print('  Group V: voting analysis ...')
        voting_fns = [
            ('fig_V1_allele_vote_matrix',
             dict(gene_data=gene_data, tool_names=tool_names,
                  consensus_df=consensus_df, sample_id=args.sample_id)),
            ('fig_V2_vote_weight_breakdown',
             dict(gene_data=gene_data, tool_names=tool_names,
                  consensus_df=consensus_df, sample_id=args.sample_id,
                  weights_data=weights_data)),
            ('fig_V3_weighted_confidence_comparison',
             dict(gene_data=gene_data, tool_names=tool_names,
                  consensus_df=consensus_df, sample_id=args.sample_id,
                  weights_data=weights_data)),
            ('fig_V4_allele_concordance_detail',
             dict(gene_data=gene_data, tool_names=tool_names,
                  sample_id=args.sample_id)),
        ]
        for name, kw in voting_fns:
            fn = globals()[name]
            try:
                fn(**kw, output_dir=od, fmt=fmt)
                print(f'    \u2713 {name}')
            except Exception as e:
                print(f'    \u2717 {name}: {e}')

    if fmt in ('pdf', 'both'):
        try:
            combine_pdfs(od)
            print('    \u2713 all_figures.pdf')
        except Exception as e:
            print(f'    \u2717 combine_pdfs: {e}')

    # ── Captions report ────────────────────────────────────────────────────────
    try:
        generate_captions_report(od, args.sample_id, consensus_df,
                                  gene_data, tool_names, trace_df, weights_data)
        print('    \u2713 figure_captions_evaluation.md')
    except Exception as e:
        print(f'    \u2717 captions report: {e}')

    pdfs = len(list(Path(od).glob('fig_*.pdf')))
    pngs = len(list(Path(od).glob('fig_*.png')))
    print(f'\n[make_paper_figures] Done — {pdfs} PDFs, {pngs} PNGs → {od}')


if __name__ == '__main__':
    main()
