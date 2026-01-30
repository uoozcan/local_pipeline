#!/usr/bin/env python3
"""
HLA LOH Summary Report Generator

Aggregates LOH results from multiple samples into a comprehensive report.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Any
from collections import Counter
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def parse_args():
    parser = argparse.ArgumentParser(description='Generate HLA LOH summary report')
    parser.add_argument('--loh-files', nargs='+', required=True,
                        help='LOH result files from multiple samples')
    parser.add_argument('--output-tsv', required=True, help='Output TSV summary')
    parser.add_argument('--output-html', required=True, help='Output HTML report')
    return parser.parse_args()


def parse_loh_file(filepath: str) -> List[Dict[str, Any]]:
    """Parse a single LOH results file."""
    results = []

    try:
        with open(filepath, 'r') as f:
            header = None
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')

                if header is None:
                    header = parts
                    continue

                if len(parts) >= len(header):
                    row = dict(zip(header, parts))
                    results.append(row)

    except Exception as e:
        print(f"Warning: Error parsing {filepath}: {e}", file=sys.stderr)

    return results


def create_summary_plots(all_results: List[Dict], output_dir: str):
    """Create summary plots for LOH analysis."""
    if not HAS_MATPLOTLIB:
        return []

    plots = []

    # Aggregate data
    samples = list(set(r.get('Sample', '') for r in all_results))
    genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

    # Calculate LOH frequency per gene
    gene_loh_counts = Counter()
    gene_total_counts = Counter()

    for r in all_results:
        gene = r.get('HLA', '')
        if gene in genes:
            gene_total_counts[gene] += 1
            if r.get('LOH') == 'Y':
                gene_loh_counts[gene] += 1

    # Plot 1: LOH frequency by gene
    fig, ax = plt.subplots(figsize=(10, 6))

    x = range(len(genes))
    frequencies = [gene_loh_counts.get(g, 0) / max(gene_total_counts.get(g, 1), 1) * 100
                   for g in genes]

    bars = ax.bar(x, frequencies, color='#e74c3c', alpha=0.8)
    ax.set_xlabel('HLA Gene', fontsize=12)
    ax.set_ylabel('LOH Frequency (%)', fontsize=12)
    ax.set_title('HLA LOH Frequency Across Samples', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([f'HLA-{g}' for g in genes], rotation=45, ha='right')
    ax.set_ylim(0, 100)

    # Add count labels
    for i, (bar, freq) in enumerate(zip(bars, frequencies)):
        count = gene_loh_counts.get(genes[i], 0)
        total = gene_total_counts.get(genes[i], 0)
        ax.annotate(f'{count}/{total}', (bar.get_x() + bar.get_width()/2, bar.get_height()),
                   ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plot1_path = f"{output_dir}/loh_frequency_by_gene.png"
    plt.savefig(plot1_path, dpi=150, bbox_inches='tight')
    plt.close()
    plots.append(plot1_path)

    # Plot 2: Sample-level LOH count
    if len(samples) > 1:
        sample_loh_counts = Counter()
        for r in all_results:
            if r.get('LOH') == 'Y':
                sample_loh_counts[r.get('Sample', '')] += 1

        fig, ax = plt.subplots(figsize=(max(10, len(samples) * 0.5), 6))

        sample_names = sorted(samples)
        counts = [sample_loh_counts.get(s, 0) for s in sample_names]

        colors = ['#e74c3c' if c > 0 else '#2ecc71' for c in counts]
        ax.bar(range(len(sample_names)), counts, color=colors, alpha=0.8)
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Number of LOH Events', fontsize=12)
        ax.set_title('LOH Events per Sample', fontsize=14)
        ax.set_xticks(range(len(sample_names)))
        ax.set_xticklabels(sample_names, rotation=90, ha='center', fontsize=8)

        plt.tight_layout()
        plot2_path = f"{output_dir}/loh_by_sample.png"
        plt.savefig(plot2_path, dpi=150, bbox_inches='tight')
        plt.close()
        plots.append(plot2_path)

    return plots


def create_html_report(all_results: List[Dict], samples: List[str],
                       plots: List[str], output_path: str):
    """Create comprehensive HTML summary report."""

    # Calculate statistics
    total_samples = len(samples)
    samples_with_loh = len(set(r.get('Sample') for r in all_results if r.get('LOH') == 'Y'))
    total_loh_events = sum(1 for r in all_results if r.get('LOH') == 'Y')

    # Gene-level statistics
    genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    gene_stats = {}
    for gene in genes:
        gene_results = [r for r in all_results if r.get('HLA') == gene]
        loh_count = sum(1 for r in gene_results if r.get('LOH') == 'Y')
        total = len(gene_results)
        gene_stats[gene] = {'loh': loh_count, 'total': total,
                           'freq': loh_count / total * 100 if total > 0 else 0}

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>HLA LOH Summary Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #e74c3c;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
        }}
        .summary-card.warning {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .summary-card.success {{
            background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        }}
        .card-value {{
            font-size: 2.5em;
            font-weight: bold;
        }}
        .card-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: #34495e;
            color: white;
        }}
        tr:hover {{
            background: #f5f5f5;
        }}
        .loh-yes {{
            background: #fadbd8;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .gene-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .gene-card {{
            background: #ecf0f1;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .gene-card.has-loh {{
            background: #fadbd8;
            border: 2px solid #e74c3c;
        }}
        .footer {{
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #7f8c8d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>HLA Loss of Heterozygosity - Multi-Sample Summary</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="card-value">{total_samples}</div>
                <div class="card-label">Total Samples</div>
            </div>
            <div class="summary-card {'warning' if samples_with_loh > 0 else 'success'}">
                <div class="card-value">{samples_with_loh}</div>
                <div class="card-label">Samples with LOH</div>
            </div>
            <div class="summary-card {'warning' if total_loh_events > 0 else 'success'}">
                <div class="card-value">{total_loh_events}</div>
                <div class="card-label">Total LOH Events</div>
            </div>
            <div class="summary-card">
                <div class="card-value">{total_loh_events / max(total_samples, 1):.1f}</div>
                <div class="card-label">Avg LOH per Sample</div>
            </div>
        </div>

        <h2>LOH Frequency by HLA Gene</h2>
        <div class="gene-grid">
"""

    for gene in genes:
        stats = gene_stats[gene]
        has_loh = stats['loh'] > 0
        html += f"""
            <div class="gene-card {'has-loh' if has_loh else ''}">
                <div style="font-size: 1.5em; font-weight: bold;">HLA-{gene}</div>
                <div style="font-size: 2em; color: {'#e74c3c' if has_loh else '#27ae60'};">{stats['freq']:.0f}%</div>
                <div style="font-size: 0.9em; color: #666;">{stats['loh']}/{stats['total']} samples</div>
            </div>
"""

    html += """
        </div>
"""

    # Add plots
    for plot in plots:
        plot_name = Path(plot).name
        html += f"""
        <div class="plot-container">
            <img src="{plot_name}" alt="LOH Summary Plot">
        </div>
"""

    # Sample details table
    html += """
        <h2>Sample Details</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>HLA Gene</th>
                <th>Allele 1</th>
                <th>Allele 2</th>
                <th>Copy Ratio</th>
                <th>LOH</th>
            </tr>
"""

    for r in sorted(all_results, key=lambda x: (x.get('Sample', ''), x.get('HLA', ''))):
        row_class = 'loh-yes' if r.get('LOH') == 'Y' else ''
        html += f"""
            <tr class="{row_class}">
                <td>{r.get('Sample', '')}</td>
                <td>HLA-{r.get('HLA', '')}</td>
                <td>{r.get('Allele1', '')}</td>
                <td>{r.get('Allele2', '')}</td>
                <td>{r.get('CopyRatio', '')}</td>
                <td>{r.get('LOH', '')}</td>
            </tr>
"""

    html += """
        </table>

        <div class="footer">
            <p>Generated by HLA Typing Pipeline v1.2.0</p>
        </div>
    </div>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)


def main():
    args = parse_args()

    # Parse all LOH files
    all_results = []
    samples = set()

    for filepath in args.loh_files:
        results = parse_loh_file(filepath)
        all_results.extend(results)
        for r in results:
            samples.add(r.get('Sample', ''))

    samples = list(samples)

    if not all_results:
        print("Warning: No LOH results found", file=sys.stderr)

    # Write TSV summary
    with open(args.output_tsv, 'w') as f:
        f.write("Sample\tHLA\tAllele1\tAllele2\tCopyRatio\tKeptHLA\tLostHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH\n")
        for r in all_results:
            f.write(f"{r.get('Sample', '')}\t{r.get('HLA', '')}\t")
            f.write(f"{r.get('Allele1', '')}\t{r.get('Allele2', '')}\t")
            f.write(f"{r.get('CopyRatio', '')}\t{r.get('KeptHLA', '')}\t")
            f.write(f"{r.get('LostHLA', '')}\t{r.get('Freq1', '')}\t")
            f.write(f"{r.get('Freq2', '')}\t{r.get('Purity', '')}\t")
            f.write(f"{r.get('Het_num', '')}\t{r.get('LOH', '')}\n")

    # Create plots
    output_dir = str(Path(args.output_html).parent)
    plots = create_summary_plots(all_results, output_dir)

    # Create HTML report
    create_html_report(all_results, samples, plots, args.output_html)

    loh_count = sum(1 for r in all_results if r.get('LOH') == 'Y')
    print(f"LOH summary complete: {loh_count} LOH events in {len(samples)} samples")


if __name__ == '__main__':
    main()
