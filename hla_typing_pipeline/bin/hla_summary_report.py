#!/usr/bin/env python3
"""
HLA Summary Report Generator
Aggregates results from multiple samples into a comprehensive report

Creates:
- Summary statistics table
- Cross-sample allele frequency plot
- Sample quality comparison
- Interactive HTML summary report
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any
from collections import Counter

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


def parse_args():
    parser = argparse.ArgumentParser(description='HLA Summary Report Generator')
    parser.add_argument('--consensus-files', nargs='+', required=True, help='Consensus result files')
    parser.add_argument('--comparison-files', nargs='+', help='Comparison result files')
    parser.add_argument('--statistics-files', nargs='+', help='Statistics JSON files')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    return parser.parse_args()


def parse_statistics_file(filepath: str) -> Dict:
    """Parse statistics JSON file."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Warning: Could not parse {filepath}: {e}", file=sys.stderr)
        return {}


def parse_consensus_file(filepath: str) -> Dict:
    """Parse consensus file to extract sample info."""
    data = {'sample_id': Path(filepath).stem.replace('_consensus', ''), 'alleles': []}

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    data['alleles'].append({
                        'gene': parts[0],
                        'allele1': parts[1] if parts[1] != '-' else None,
                        'allele2': parts[2] if parts[2] != '-' else None,
                        'confidence': float(parts[3])
                    })
    except Exception as e:
        print(f"Warning: Could not parse {filepath}: {e}", file=sys.stderr)

    return data


def calculate_allele_frequencies(all_data: List[Dict]) -> Dict[str, Counter]:
    """Calculate allele frequencies across all samples."""
    gene_alleles = {}

    for sample_data in all_data:
        for allele_info in sample_data.get('alleles', []):
            gene = allele_info['gene']
            if gene not in gene_alleles:
                gene_alleles[gene] = Counter()

            if allele_info['allele1']:
                gene_alleles[gene][allele_info['allele1']] += 1
            if allele_info['allele2']:
                gene_alleles[gene][allele_info['allele2']] += 1

    return gene_alleles


def plot_allele_frequency(allele_freqs: Dict[str, Counter], output_dir: str, top_n: int = 5):
    """Create allele frequency plots for each gene."""
    if not HAS_MATPLOTLIB:
        return

    genes = sorted(allele_freqs.keys())
    n_genes = len(genes)

    if n_genes == 0:
        return

    # Create subplot grid
    cols = min(3, n_genes)
    rows = (n_genes + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_genes == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_genes > 1 else [axes]

    for idx, gene in enumerate(genes):
        ax = axes[idx]
        freq_counter = allele_freqs[gene]

        # Get top alleles
        top_alleles = freq_counter.most_common(top_n)
        if not top_alleles:
            ax.set_visible(False)
            continue

        alleles = [a[0].split('*')[1] if '*' in a[0] else a[0] for a in top_alleles]
        counts = [a[1] for a in top_alleles]

        colors = plt.cm.Set3(np.linspace(0, 1, len(alleles)))
        bars = ax.bar(alleles, counts, color=colors, edgecolor='black', linewidth=0.5)

        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('Allele')
        ax.set_ylabel('Frequency')
        ax.tick_params(axis='x', rotation=45)

        # Add count labels
        for bar, count in zip(bars, counts):
            ax.annotate(str(count), xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom',
                        fontsize=9)

    # Hide unused subplots
    for idx in range(len(genes), len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('HLA Allele Frequencies Across Samples', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/allele_frequency_summary.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_sample_quality_comparison(all_stats: List[Dict], output_dir: str):
    """Create sample quality comparison plot."""
    if not HAS_MATPLOTLIB or not all_stats:
        return

    samples = []
    confidences = []
    hla_reads = []
    qc_status = []

    for stats in all_stats:
        if 'sample_id' in stats:
            samples.append(stats['sample_id'])
            confidences.append(stats.get('summary', {}).get('avg_confidence', 0))
            hla_reads.append(stats.get('qc', {}).get('hla_reads', 0))
            status = stats.get('qc', {}).get('status', 'unknown')
            qc_status.append(status)

    if not samples:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Average confidence per sample
    ax1 = axes[0]
    colors = ['#2ecc71' if c >= 0.8 else '#f39c12' if c >= 0.5 else '#e74c3c' for c in confidences]
    bars = ax1.bar(samples, confidences, color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Sample', fontsize=12)
    ax1.set_ylabel('Average Confidence', fontsize=12)
    ax1.set_title('Average Typing Confidence per Sample', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 1.1)
    ax1.axhline(y=0.8, color='green', linestyle='--', alpha=0.5)
    ax1.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5)
    ax1.tick_params(axis='x', rotation=45)

    for bar, conf in zip(bars, confidences):
        ax1.annotate(f'{conf:.2f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                     xytext=(0, 3), textcoords='offset points', ha='center', va='bottom',
                     fontsize=9, fontweight='bold')

    # Plot 2: HLA reads per sample
    ax2 = axes[1]
    colors2 = ['#2ecc71' if s == 'true' else '#f39c12' if s == 'warning' else '#e74c3c' for s in qc_status]
    bars2 = ax2.bar(samples, hla_reads, color=colors2, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Sample', fontsize=12)
    ax2.set_ylabel('HLA Region Reads', fontsize=12)
    ax2.set_title('HLA Read Count per Sample', fontsize=12, fontweight='bold')
    ax2.tick_params(axis='x', rotation=45)

    # Format y-axis with K for thousands
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1000:.0f}K' if x >= 1000 else f'{x:.0f}'))

    plt.tight_layout()
    plt.savefig(f'{output_dir}/sample_quality_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_confidence_distribution(all_stats: List[Dict], output_dir: str):
    """Create confidence score distribution across all samples."""
    if not HAS_MATPLOTLIB or not all_stats:
        return

    all_confidences = []
    for stats in all_stats:
        for gene_data in stats.get('per_gene', {}).values():
            all_confidences.append(gene_data.get('confidence', 0))

    if not all_confidences:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    # Histogram
    n, bins, patches = ax.hist(all_confidences, bins=20, edgecolor='black', linewidth=0.5, alpha=0.7)

    # Color bins by confidence level
    for i, patch in enumerate(patches):
        bin_center = (bins[i] + bins[i+1]) / 2
        if bin_center >= 0.8:
            patch.set_facecolor('#2ecc71')
        elif bin_center >= 0.5:
            patch.set_facecolor('#f39c12')
        else:
            patch.set_facecolor('#e74c3c')

    ax.axvline(x=0.8, color='green', linestyle='--', linewidth=2, label='High confidence (0.8)')
    ax.axvline(x=0.5, color='orange', linestyle='--', linewidth=2, label='Medium confidence (0.5)')

    ax.set_xlabel('Confidence Score', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of HLA Typing Confidence Scores', fontsize=14, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{output_dir}/confidence_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()


def generate_summary_table(all_stats: List[Dict], all_data: List[Dict], output_dir: str):
    """Generate summary statistics table."""

    rows = []
    for stats, data in zip(all_stats, all_data):
        sample_id = stats.get('sample_id', data.get('sample_id', 'Unknown'))
        summary = stats.get('summary', {})
        qc = stats.get('qc', {})

        row = {
            'Sample': sample_id,
            'Genes_Typed': summary.get('total_genes', len(data.get('alleles', []))),
            'High_Conf': summary.get('high_confidence_genes', 0),
            'Med_Conf': summary.get('medium_confidence_genes', 0),
            'Low_Conf': summary.get('low_confidence_genes', 0),
            'Avg_Confidence': f"{summary.get('avg_confidence', 0):.2f}",
            'HLA_Reads': qc.get('hla_reads', 0),
            'QC_Status': qc.get('status', 'unknown')
        }
        rows.append(row)

    # Write TSV
    if rows:
        headers = list(rows[0].keys())
        with open(f'{output_dir}/hla_summary_statistics.tsv', 'w') as f:
            f.write('\t'.join(headers) + '\n')
            for row in rows:
                f.write('\t'.join(str(row[h]) for h in headers) + '\n')


def generate_html_summary(all_stats: List[Dict], all_data: List[Dict],
                          allele_freqs: Dict, output_dir: str):
    """Generate comprehensive HTML summary report."""

    n_samples = len(all_data)
    all_confidences = []
    for stats in all_stats:
        all_confidences.extend([g.get('confidence', 0) for g in stats.get('per_gene', {}).values()])

    avg_confidence = sum(all_confidences) / len(all_confidences) if all_confidences else 0
    high_conf_pct = sum(1 for c in all_confidences if c >= 0.8) / len(all_confidences) * 100 if all_confidences else 0

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HLA Typing Summary Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f6fa;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .card {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .metric {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
        }}
        .metric-label {{
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
            background-color: #3498db;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .status-true, .status-pass {{
            color: #27ae60;
            font-weight: bold;
        }}
        .status-warning {{
            color: #f39c12;
            font-weight: bold;
        }}
        .status-critical, .status-false {{
            color: #e74c3c;
            font-weight: bold;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 20px;
        }}
        .footer {{
            text-align: center;
            margin-top: 40px;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 HLA Typing Summary Report</h1>
        <p><strong>Total Samples:</strong> {n_samples}</p>
        <p><strong>Generated:</strong> <span id="date"></span></p>

        <h2>📊 Overall Statistics</h2>
        <div class="metrics-grid">
            <div class="metric">
                <div class="metric-value">{n_samples}</div>
                <div class="metric-label">Total Samples</div>
            </div>
            <div class="metric">
                <div class="metric-value">{avg_confidence:.2f}</div>
                <div class="metric-label">Average Confidence</div>
            </div>
            <div class="metric">
                <div class="metric-value">{high_conf_pct:.1f}%</div>
                <div class="metric-label">High Confidence Calls</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(allele_freqs)}</div>
                <div class="metric-label">HLA Genes Analyzed</div>
            </div>
        </div>

        <div class="card">
            <h2>📋 Sample Summary</h2>
            <table>
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Genes Typed</th>
                        <th>High Conf</th>
                        <th>Med Conf</th>
                        <th>Low Conf</th>
                        <th>Avg Confidence</th>
                        <th>HLA Reads</th>
                        <th>QC Status</th>
                    </tr>
                </thead>
                <tbody>
"""

    for stats, data in zip(all_stats, all_data):
        sample_id = stats.get('sample_id', data.get('sample_id', 'Unknown'))
        summary = stats.get('summary', {})
        qc = stats.get('qc', {})
        status = qc.get('status', 'unknown')
        status_class = f'status-{status}'

        html_content += f"""                    <tr>
                        <td><strong>{sample_id}</strong></td>
                        <td>{summary.get('total_genes', len(data.get('alleles', [])))}</td>
                        <td>{summary.get('high_confidence_genes', 0)}</td>
                        <td>{summary.get('medium_confidence_genes', 0)}</td>
                        <td>{summary.get('low_confidence_genes', 0)}</td>
                        <td>{summary.get('avg_confidence', 0):.2f}</td>
                        <td>{qc.get('hla_reads', 0):,}</td>
                        <td class="{status_class}">{status.upper()}</td>
                    </tr>
"""

    html_content += """                </tbody>
            </table>
        </div>

        <div class="card">
            <h2>📈 Visualizations</h2>
            <div class="plot-grid">
                <div class="plot-container">
                    <h3>Sample Quality Comparison</h3>
                    <img src="sample_quality_comparison.png" alt="Sample Quality">
                </div>
                <div class="plot-container">
                    <h3>Confidence Distribution</h3>
                    <img src="confidence_distribution.png" alt="Confidence Distribution">
                </div>
            </div>
            <div class="plot-container">
                <h3>Allele Frequencies</h3>
                <img src="allele_frequency_summary.png" alt="Allele Frequencies">
            </div>
        </div>

        <div class="card">
            <h2>🔬 Most Common Alleles</h2>
            <table>
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Most Common Allele</th>
                        <th>Frequency</th>
                        <th>Second Most Common</th>
                        <th>Frequency</th>
                    </tr>
                </thead>
                <tbody>
"""

    for gene in sorted(allele_freqs.keys()):
        top2 = allele_freqs[gene].most_common(2)
        allele1 = top2[0][0] if len(top2) > 0 else '-'
        freq1 = top2[0][1] if len(top2) > 0 else 0
        allele2 = top2[1][0] if len(top2) > 1 else '-'
        freq2 = top2[1][1] if len(top2) > 1 else 0

        html_content += f"""                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{allele1}</td>
                        <td>{freq1}</td>
                        <td>{allele2}</td>
                        <td>{freq2}</td>
                    </tr>
"""

    html_content += """                </tbody>
            </table>
        </div>

        <div class="footer">
            <p>Generated by HLA Typing Pipeline v1.1.0</p>
        </div>
    </div>

    <script>
        document.getElementById('date').textContent = new Date().toLocaleString();
    </script>
</body>
</html>
"""

    with open(f'{output_dir}/hla_summary_report.html', 'w') as f:
        f.write(html_content)


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse all consensus files
    all_data = []
    for cf in args.consensus_files:
        if Path(cf).exists():
            all_data.append(parse_consensus_file(cf))

    # Parse statistics files if provided
    all_stats = []
    if args.statistics_files:
        for sf in args.statistics_files:
            if Path(sf).exists():
                all_stats.append(parse_statistics_file(sf))
    else:
        # Create minimal stats from consensus data
        for data in all_data:
            stats = {
                'sample_id': data['sample_id'],
                'summary': {
                    'total_genes': len(data['alleles']),
                    'high_confidence_genes': sum(1 for a in data['alleles'] if a['confidence'] >= 0.8),
                    'medium_confidence_genes': sum(1 for a in data['alleles'] if 0.5 <= a['confidence'] < 0.8),
                    'low_confidence_genes': sum(1 for a in data['alleles'] if a['confidence'] < 0.5),
                    'avg_confidence': sum(a['confidence'] for a in data['alleles']) / len(data['alleles']) if data['alleles'] else 0
                },
                'qc': {'hla_reads': 0, 'status': 'unknown'},
                'per_gene': {a['gene']: {'confidence': a['confidence']} for a in data['alleles']}
            }
            all_stats.append(stats)

    if not all_data:
        print("Error: No valid data files found", file=sys.stderr)
        sys.exit(1)

    print(f"Generating summary report for {len(all_data)} samples...")

    # Calculate allele frequencies
    allele_freqs = calculate_allele_frequencies(all_data)

    # Generate plots
    plot_allele_frequency(allele_freqs, str(output_dir))
    plot_sample_quality_comparison(all_stats, str(output_dir))
    plot_confidence_distribution(all_stats, str(output_dir))

    # Generate summary table
    generate_summary_table(all_stats, all_data, str(output_dir))

    # Generate HTML report
    generate_html_summary(all_stats, all_data, allele_freqs, str(output_dir))

    print(f"Summary report saved to {output_dir}/")


if __name__ == '__main__':
    main()
