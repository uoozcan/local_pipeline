#!/usr/bin/env python3
"""
HLA Typing Visualization Script
Generates plots and statistics for HLA typing results

Creates:
- Confidence bar chart per gene
- Tool agreement heatmap
- Read coverage plot
- Interactive HTML report
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
except ImportError:
    print("Error: matplotlib/numpy required. Install with: pip install matplotlib numpy", file=sys.stderr)
    sys.exit(1)

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("Warning: pandas not available", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(description='HLA Typing Visualization')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--consensus', required=True, help='Consensus results file')
    parser.add_argument('--comparison', required=True, help='Comparison results file')
    parser.add_argument('--qc-report', help='QC report file')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    return parser.parse_args()


def parse_consensus_file(filepath: str) -> Dict[str, Any]:
    """Parse consensus results file."""
    results = {
        'metadata': {},
        'alleles': []
    }

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Tools used:'):
                results['metadata']['tools'] = line.split(':')[1].strip()
            elif line.startswith('# Resolution:'):
                results['metadata']['resolution'] = line.split(':')[1].strip()
            elif line.startswith('# Weighting'):
                results['metadata']['weighting'] = line.split(':')[1].strip()
            elif line.startswith('# Expected reads'):
                results['metadata']['expected_reads'] = line.split(':')[1].strip()
            elif line.startswith('#') or line.startswith('Gene'):
                continue
            elif line:
                parts = line.split('\t')
                if len(parts) >= 4:
                    allele_data = {
                        'gene': parts[0],
                        'allele1': parts[1] if parts[1] != '-' else None,
                        'allele2': parts[2] if parts[2] != '-' else None,
                        'confidence': float(parts[3]),
                        'reads1': int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                        'reads2': int(parts[5]) if len(parts) > 5 and parts[5].isdigit() else 0
                    }
                    results['alleles'].append(allele_data)

    return results


def parse_comparison_file(filepath: str) -> Dict[str, Any]:
    """Parse comparison results file."""
    results = {
        'metadata': {},
        'tools': [],
        'data': []
    }

    header = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Resolution:'):
                results['metadata']['resolution'] = line.split(':')[1].strip()
            elif line.startswith('# Weighting'):
                results['metadata']['weighting'] = line.split(':')[1].strip()
            elif line.startswith('#'):
                continue
            elif line.startswith('Gene'):
                header = line.split('\t')
                # Extract tool names from header
                for col in header:
                    if '(alleles)' in col:
                        tool = col.replace('(alleles)', '')
                        results['tools'].append(tool)
            elif line and header:
                parts = line.split('\t')
                results['data'].append(dict(zip(header, parts)))

    return results


def parse_qc_report(filepath: str) -> Dict[str, Any]:
    """Parse QC report file."""
    qc_data = {
        'total_reads': 0,
        'hla_reads': 0,
        'avg_read_length': 0,
        'avg_mapq': 0,
        'status': 'unknown',
        'warnings': []
    }

    if not filepath or not Path(filepath).exists():
        return qc_data

    with open(filepath, 'r') as f:
        content = f.read()

        # Parse key metrics
        for line in content.split('\n'):
            if 'Total reads:' in line:
                try:
                    qc_data['total_reads'] = int(line.split(':')[1].strip())
                except:
                    pass
            elif 'HLA region reads:' in line:
                try:
                    qc_data['hla_reads'] = int(line.split(':')[1].strip())
                except:
                    pass
            elif 'Average read length:' in line:
                try:
                    qc_data['avg_read_length'] = int(line.split(':')[1].strip().replace(' bp', ''))
                except:
                    pass
            elif 'Average mapping quality:' in line:
                try:
                    qc_data['avg_mapq'] = float(line.split(':')[1].strip())
                except:
                    pass
            elif 'Status:' in line:
                qc_data['status'] = line.split(':')[1].strip()
            elif 'WARNING:' in line or 'CRITICAL:' in line:
                qc_data['warnings'].append(line.strip())

    return qc_data


def plot_confidence_chart(consensus_data: Dict, sample_id: str, output_dir: str):
    """Create bar chart showing confidence scores per HLA gene."""
    genes = [a['gene'] for a in consensus_data['alleles']]
    confidences = [a['confidence'] for a in consensus_data['alleles']]

    # Color based on confidence level
    colors = []
    for conf in confidences:
        if conf >= 0.8:
            colors.append('#2ecc71')  # Green - high confidence
        elif conf >= 0.5:
            colors.append('#f39c12')  # Orange - medium confidence
        else:
            colors.append('#e74c3c')  # Red - low confidence

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(genes, confidences, color=colors, edgecolor='black', linewidth=0.5)

    # Add value labels on bars
    for bar, conf in zip(bars, confidences):
        height = bar.get_height()
        ax.annotate(f'{conf:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_xlabel('HLA Gene', fontsize=12)
    ax.set_ylabel('Confidence Score', fontsize=12)
    ax.set_title(f'HLA Typing Confidence Scores - {sample_id}', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.15)
    ax.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='High confidence threshold')
    ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.5, label='Medium confidence threshold')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#2ecc71', label='High (≥0.8)'),
        mpatches.Patch(facecolor='#f39c12', label='Medium (0.5-0.8)'),
        mpatches.Patch(facecolor='#e74c3c', label='Low (<0.5)')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_confidence_chart.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_read_coverage(consensus_data: Dict, sample_id: str, output_dir: str):
    """Create bar chart showing read counts per allele."""
    genes = []
    reads1 = []
    reads2 = []

    for a in consensus_data['alleles']:
        if a['allele1'] or a['allele2']:
            genes.append(a['gene'])
            reads1.append(a['reads1'])
            reads2.append(a['reads2'])

    if not genes:
        return

    x = np.arange(len(genes))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))
    bars1 = ax.bar(x - width/2, reads1, width, label='Allele 1', color='#3498db', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, reads2, width, label='Allele 2', color='#e74c3c', edgecolor='black', linewidth=0.5)

    ax.set_xlabel('HLA Gene', fontsize=12)
    ax.set_ylabel('Read Count', fontsize=12)
    ax.set_title(f'HLA Allele Read Coverage - {sample_id}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right')
    ax.legend()

    # Add expected reads line
    expected = int(consensus_data.get('metadata', {}).get('expected_reads', 1000))
    ax.axhline(y=expected, color='green', linestyle='--', alpha=0.7, label=f'Expected ({expected})')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_read_coverage.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_tool_agreement(comparison_data: Dict, sample_id: str, output_dir: str):
    """Create heatmap showing tool agreement across genes."""
    tools = comparison_data['tools']
    if not tools or not comparison_data['data']:
        return

    genes = [d['Gene'] for d in comparison_data['data']]

    # Build agreement matrix
    # 1 = allele called, 0 = no call, 0.5 = partial match with consensus
    agreement_matrix = []

    for row in comparison_data['data']:
        gene_agreement = []
        consensus = row.get('Consensus', '-')

        for tool in tools:
            tool_alleles = row.get(f'{tool}(alleles)', '-')
            if tool_alleles == '-':
                gene_agreement.append(0)
            elif consensus != '-' and tool_alleles == consensus:
                gene_agreement.append(1)
            elif consensus != '-':
                # Check partial match
                consensus_set = set(consensus.split('/'))
                tool_set = set(tool_alleles.split('/'))
                overlap = len(consensus_set & tool_set)
                gene_agreement.append(overlap / max(len(consensus_set), 1))
            else:
                gene_agreement.append(0.5)

        agreement_matrix.append(gene_agreement)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.imshow(agreement_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(tools)))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_xticklabels(tools, fontsize=10)
    ax.set_yticklabels(genes, fontsize=10)

    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Agreement with Consensus', rotation=-90, va='bottom', fontsize=10)

    # Add text annotations
    for i in range(len(genes)):
        for j in range(len(tools)):
            value = agreement_matrix[i][j]
            text_color = 'white' if value < 0.5 else 'black'
            ax.text(j, i, f'{value:.1f}', ha='center', va='center', color=text_color, fontsize=9)

    ax.set_title(f'Tool Agreement Heatmap - {sample_id}', fontsize=14, fontweight='bold')
    ax.set_xlabel('HLA Typing Tool', fontsize=12)
    ax.set_ylabel('HLA Gene', fontsize=12)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_tool_agreement.png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_qc_summary(qc_data: Dict, sample_id: str, output_dir: str):
    """Create QC summary visualization."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: Read distribution pie chart
    ax1 = axes[0]
    if qc_data['total_reads'] > 0 and qc_data['hla_reads'] > 0:
        other_reads = qc_data['total_reads'] - qc_data['hla_reads']
        sizes = [qc_data['hla_reads'], other_reads]
        labels = [f"HLA Region\n({qc_data['hla_reads']:,})", f"Other\n({other_reads:,})"]
        colors = ['#3498db', '#95a5a6']
        explode = (0.05, 0)

        ax1.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax1.set_title('Read Distribution', fontsize=12, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'No read data\navailable', ha='center', va='center', fontsize=14)
        ax1.set_title('Read Distribution', fontsize=12, fontweight='bold')

    # Plot 2: QC metrics bar chart
    ax2 = axes[1]
    metrics = ['Read Length', 'Map Quality']
    values = [qc_data['avg_read_length'], qc_data['avg_mapq']]
    thresholds = [50, 20]  # Minimum acceptable values

    x = np.arange(len(metrics))
    bars = ax2.bar(x, values, color=['#3498db', '#2ecc71'], edgecolor='black', linewidth=0.5)

    # Add threshold lines
    for i, (val, thresh) in enumerate(zip(values, thresholds)):
        ax2.axhline(y=thresh, xmin=i/len(metrics)-0.1, xmax=(i+1)/len(metrics)+0.1,
                    color='red', linestyle='--', alpha=0.7)

    ax2.set_xticks(x)
    ax2.set_xticklabels(metrics)
    ax2.set_ylabel('Value')
    ax2.set_title('QC Metrics', fontsize=12, fontweight='bold')

    # Add value labels
    for bar, val in zip(bars, values):
        ax2.annotate(f'{val:.1f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                     xytext=(0, 3), textcoords='offset points', ha='center', va='bottom',
                     fontsize=11, fontweight='bold')

    # Plot 3: Status indicator
    ax3 = axes[2]
    ax3.axis('off')

    status = qc_data['status']
    if status == 'true':
        status_color = '#2ecc71'
        status_text = 'PASS'
        status_icon = '✓'
    elif status == 'warning':
        status_color = '#f39c12'
        status_text = 'WARNING'
        status_icon = '⚠'
    else:
        status_color = '#e74c3c'
        status_text = 'CRITICAL'
        status_icon = '✗'

    circle = plt.Circle((0.5, 0.6), 0.25, color=status_color, ec='black', linewidth=2)
    ax3.add_patch(circle)
    ax3.text(0.5, 0.6, status_icon, ha='center', va='center', fontsize=40, color='white', fontweight='bold')
    ax3.text(0.5, 0.2, status_text, ha='center', va='center', fontsize=16, fontweight='bold')
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.set_title('QC Status', fontsize=12, fontweight='bold')

    plt.suptitle(f'QC Summary - {sample_id}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{sample_id}_qc_summary.png', dpi=150, bbox_inches='tight')
    plt.close()


def generate_html_report(sample_id: str, consensus_data: Dict, comparison_data: Dict,
                         qc_data: Dict, output_dir: str):
    """Generate interactive HTML report."""

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HLA Typing Report - {sample_id}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f6fa;
        }}
        .container {{
            max-width: 1200px;
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
        .confidence-high {{
            background-color: #d5f5e3;
            color: #1e8449;
        }}
        .confidence-medium {{
            background-color: #fef9e7;
            color: #9a7d0a;
        }}
        .confidence-low {{
            background-color: #fadbd8;
            color: #922b21;
        }}
        .status-pass {{
            color: #27ae60;
            font-weight: bold;
        }}
        .status-warning {{
            color: #f39c12;
            font-weight: bold;
        }}
        .status-critical {{
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
        .warning-box {{
            background-color: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 5px;
            padding: 15px;
            margin: 10px 0;
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
        <h1>🧬 HLA Typing Report</h1>
        <p><strong>Sample ID:</strong> {sample_id}</p>
        <p><strong>Generated:</strong> <span id="date"></span></p>

        <h2>📊 Summary Metrics</h2>
        <div class="metrics-grid">
            <div class="metric">
                <div class="metric-value">{len(consensus_data['alleles'])}</div>
                <div class="metric-label">HLA Genes Typed</div>
            </div>
            <div class="metric">
                <div class="metric-value">{sum(1 for a in consensus_data['alleles'] if a['confidence'] >= 0.8)}</div>
                <div class="metric-label">High Confidence</div>
            </div>
            <div class="metric">
                <div class="metric-value">{qc_data['hla_reads']:,}</div>
                <div class="metric-label">HLA Reads</div>
            </div>
            <div class="metric">
                <div class="metric-value">{qc_data['avg_read_length']}</div>
                <div class="metric-label">Avg Read Length</div>
            </div>
        </div>

        <div class="card">
            <h2>🎯 QC Status</h2>
            <p>Status: <span class="status-{qc_data['status']}">{qc_data['status'].upper()}</span></p>
"""

    # Add warnings if any
    if qc_data['warnings']:
        html_content += """            <div class="warning-box">
                <strong>⚠️ Warnings:</strong>
                <ul>
"""
        for warning in qc_data['warnings']:
            html_content += f"                    <li>{warning}</li>\n"
        html_content += """                </ul>
            </div>
"""

    html_content += """        </div>

        <div class="card">
            <h2>🧬 HLA Typing Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Gene</th>
                        <th>Allele 1</th>
                        <th>Allele 2</th>
                        <th>Confidence</th>
                        <th>Reads (A1/A2)</th>
                    </tr>
                </thead>
                <tbody>
"""

    for allele in consensus_data['alleles']:
        conf = allele['confidence']
        if conf >= 0.8:
            conf_class = 'confidence-high'
        elif conf >= 0.5:
            conf_class = 'confidence-medium'
        else:
            conf_class = 'confidence-low'

        html_content += f"""                    <tr class="{conf_class}">
                        <td><strong>{allele['gene']}</strong></td>
                        <td>{allele['allele1'] or '-'}</td>
                        <td>{allele['allele2'] or '-'}</td>
                        <td>{conf:.2f}</td>
                        <td>{allele['reads1']} / {allele['reads2']}</td>
                    </tr>
"""

    html_content += """                </tbody>
            </table>
        </div>

        <div class="card">
            <h2>📈 Visualizations</h2>

            <div class="plot-container">
                <h3>Confidence Scores</h3>
"""
    html_content += f'                <img src="{sample_id}_confidence_chart.png" alt="Confidence Chart">\n'
    html_content += """            </div>

            <div class="plot-container">
                <h3>Read Coverage</h3>
"""
    html_content += f'                <img src="{sample_id}_read_coverage.png" alt="Read Coverage">\n'
    html_content += """            </div>

            <div class="plot-container">
                <h3>Tool Agreement</h3>
"""
    html_content += f'                <img src="{sample_id}_tool_agreement.png" alt="Tool Agreement">\n'
    html_content += """            </div>

            <div class="plot-container">
                <h3>QC Summary</h3>
"""
    html_content += f'                <img src="{sample_id}_qc_summary.png" alt="QC Summary">\n'
    html_content += """            </div>
        </div>

        <div class="card">
            <h2>ℹ️ Analysis Information</h2>
            <table>
                <tr><td><strong>Tools Used</strong></td><td>"""
    html_content += consensus_data.get('metadata', {}).get('tools', 'N/A')
    html_content += """</td></tr>
                <tr><td><strong>Resolution</strong></td><td>"""
    html_content += consensus_data.get('metadata', {}).get('resolution', 'N/A')
    html_content += """</td></tr>
                <tr><td><strong>Weighting Method</strong></td><td>"""
    html_content += consensus_data.get('metadata', {}).get('weighting', 'N/A')
    html_content += """</td></tr>
                <tr><td><strong>Expected Reads</strong></td><td>"""
    html_content += str(consensus_data.get('metadata', {}).get('expected_reads', 'N/A'))
    html_content += """</td></tr>
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

    with open(f'{output_dir}/{sample_id}_report.html', 'w') as f:
        f.write(html_content)


def generate_statistics_json(sample_id: str, consensus_data: Dict, comparison_data: Dict,
                             qc_data: Dict, output_dir: str):
    """Generate JSON file with statistics."""

    # Calculate statistics
    confidences = [a['confidence'] for a in consensus_data['alleles']]
    reads = [a['reads1'] + a['reads2'] for a in consensus_data['alleles']]

    stats = {
        'sample_id': sample_id,
        'summary': {
            'total_genes': len(consensus_data['alleles']),
            'high_confidence_genes': sum(1 for c in confidences if c >= 0.8),
            'medium_confidence_genes': sum(1 for c in confidences if 0.5 <= c < 0.8),
            'low_confidence_genes': sum(1 for c in confidences if c < 0.5),
            'avg_confidence': sum(confidences) / len(confidences) if confidences else 0,
            'total_reads': sum(reads),
            'avg_reads_per_gene': sum(reads) / len(reads) if reads else 0
        },
        'qc': {
            'total_reads': qc_data['total_reads'],
            'hla_reads': qc_data['hla_reads'],
            'hla_percentage': (qc_data['hla_reads'] / qc_data['total_reads'] * 100) if qc_data['total_reads'] > 0 else 0,
            'avg_read_length': qc_data['avg_read_length'],
            'avg_mapq': qc_data['avg_mapq'],
            'status': qc_data['status'],
            'warnings': qc_data['warnings']
        },
        'per_gene': {},
        'metadata': consensus_data.get('metadata', {})
    }

    # Per-gene statistics
    for allele in consensus_data['alleles']:
        stats['per_gene'][allele['gene']] = {
            'allele1': allele['allele1'],
            'allele2': allele['allele2'],
            'confidence': allele['confidence'],
            'reads1': allele['reads1'],
            'reads2': allele['reads2'],
            'total_reads': allele['reads1'] + allele['reads2']
        }

    with open(f'{output_dir}/{sample_id}_statistics.json', 'w') as f:
        json.dump(stats, f, indent=2)


def main():
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse input files
    consensus_data = parse_consensus_file(args.consensus)
    comparison_data = parse_comparison_file(args.comparison)
    qc_data = parse_qc_report(args.qc_report) if args.qc_report else {
        'total_reads': 0, 'hla_reads': 0, 'avg_read_length': 0,
        'avg_mapq': 0, 'status': 'unknown', 'warnings': []
    }

    # Generate plots
    print(f"Generating visualizations for {args.sample}...")

    plot_confidence_chart(consensus_data, args.sample, str(output_dir))
    plot_read_coverage(consensus_data, args.sample, str(output_dir))
    plot_tool_agreement(comparison_data, args.sample, str(output_dir))
    plot_qc_summary(qc_data, args.sample, str(output_dir))

    # Generate HTML report
    generate_html_report(args.sample, consensus_data, comparison_data, qc_data, str(output_dir))

    # Generate statistics JSON
    generate_statistics_json(args.sample, consensus_data, comparison_data, qc_data, str(output_dir))

    print(f"Visualizations saved to {output_dir}/")


if __name__ == '__main__':
    main()
