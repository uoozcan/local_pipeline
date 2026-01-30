#!/usr/bin/env python3
"""
HLA LOH Visualization Script

Generates plots and reports for HLA Loss of Heterozygosity analysis.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def parse_args():
    parser = argparse.ArgumentParser(description='Visualize HLA LOH results')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--loh-results', required=True, help='LOH results file')
    parser.add_argument('--output-plot', required=True, help='Output plot file')
    parser.add_argument('--output-report', required=True, help='Output HTML report')
    return parser.parse_args()


def parse_loh_results(filepath: str) -> List[Dict[str, Any]]:
    """Parse LOH results file."""
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
        print(f"Warning: Error parsing LOH results: {e}", file=sys.stderr)

    return results


def create_loh_plot(results: List[Dict], sample: str, output_path: str):
    """Create LOH visualization plot."""
    if not HAS_MATPLOTLIB or not results:
        return False

    # Extract data
    genes = [r.get('HLA', '') for r in results]
    loh_status = [r.get('LOH', 'N') for r in results]

    # Parse copy ratios
    copy_a = []
    copy_b = []
    for r in results:
        ratio = r.get('CopyRatio', '1:1')
        parts = ratio.split(':')
        copy_a.append(int(parts[0]) if len(parts) > 0 and parts[0].isdigit() else 1)
        copy_b.append(int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 1)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Subplot 1: Copy number bar chart
    x = np.arange(len(genes))
    width = 0.35

    colors_a = ['#e74c3c' if loh == 'Y' else '#3498db' for loh in loh_status]
    colors_b = ['#e74c3c' if loh == 'Y' else '#2ecc71' for loh in loh_status]

    bars1 = ax1.bar(x - width/2, copy_a, width, label='Allele 1', color=colors_a, alpha=0.8)
    bars2 = ax1.bar(x + width/2, copy_b, width, label='Allele 2', color=colors_b, alpha=0.8)

    ax1.set_xlabel('HLA Gene', fontsize=12)
    ax1.set_ylabel('Copy Number', fontsize=12)
    ax1.set_title(f'HLA Allele Copy Numbers - {sample}', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f'HLA-{g}' for g in genes], rotation=45, ha='right')
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Normal (1 copy)')

    # Add LOH markers
    for i, loh in enumerate(loh_status):
        if loh == 'Y':
            ax1.annotate('LOH', (x[i], max(copy_a[i], copy_b[i]) + 0.2),
                        ha='center', fontsize=9, color='red', fontweight='bold')

    ax1.legend(loc='upper right')

    # Subplot 2: LOH status summary
    loh_counts = {'LOH': sum(1 for l in loh_status if l == 'Y'),
                  'Normal': sum(1 for l in loh_status if l == 'N')}

    colors_pie = ['#e74c3c', '#2ecc71']
    if loh_counts['LOH'] == 0:
        # All normal
        ax2.pie([1], labels=['No LOH'], colors=['#2ecc71'], autopct='',
                startangle=90)
    else:
        ax2.pie(loh_counts.values(), labels=loh_counts.keys(), colors=colors_pie,
                autopct='%1.0f%%', startangle=90, explode=(0.05, 0))

    ax2.set_title(f'LOH Status Summary\n({loh_counts["LOH"]} of {len(genes)} genes)', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return True


def create_html_report(results: List[Dict], sample: str, plot_path: str, output_path: str):
    """Create HTML report for LOH analysis."""

    loh_count = sum(1 for r in results if r.get('LOH') == 'Y')
    total_genes = len(results)
    purity = results[0].get('Purity', 'N/A') if results else 'N/A'

    # Get LOH genes
    loh_genes = [r for r in results if r.get('LOH') == 'Y']

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>HLA LOH Report - {sample}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
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
        .summary-box {{
            background: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
        }}
        .summary-item {{
            text-align: center;
        }}
        .summary-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .summary-label {{
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        .loh-detected {{
            color: #e74c3c;
        }}
        .no-loh {{
            color: #27ae60;
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
            background: #3498db;
            color: white;
        }}
        tr:hover {{
            background: #f5f5f5;
        }}
        .loh-yes {{
            background: #fadbd8;
            font-weight: bold;
        }}
        .plot-container {{
            text-align: center;
            margin: 30px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .alert {{
            padding: 15px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        .alert-warning {{
            background: #fcf3cf;
            border-left: 4px solid #f1c40f;
        }}
        .alert-danger {{
            background: #fadbd8;
            border-left: 4px solid #e74c3c;
        }}
        .footer {{
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>HLA Loss of Heterozygosity Report</h1>
        <p><strong>Sample:</strong> {sample}</p>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <div class="summary-box">
            <div class="summary-item">
                <div class="summary-value">{total_genes}</div>
                <div class="summary-label">HLA Genes Analyzed</div>
            </div>
            <div class="summary-item">
                <div class="summary-value {'loh-detected' if loh_count > 0 else 'no-loh'}">{loh_count}</div>
                <div class="summary-label">LOH Events Detected</div>
            </div>
            <div class="summary-item">
                <div class="summary-value">{purity}</div>
                <div class="summary-label">Tumor Purity</div>
            </div>
        </div>
"""

    # Add alert if LOH detected
    if loh_count > 0:
        html += f"""
        <div class="alert alert-danger">
            <strong>LOH Detected!</strong> Loss of heterozygosity was detected in {loh_count} HLA gene(s).
            This may impact immune recognition and has implications for immunotherapy response.
        </div>
"""

    # Add plot if exists
    plot_file = Path(plot_path).name
    if Path(plot_path).exists():
        html += f"""
        <h2>Copy Number Visualization</h2>
        <div class="plot-container">
            <img src="{plot_file}" alt="HLA LOH Plot">
        </div>
"""

    # Add results table
    html += """
        <h2>Detailed Results</h2>
        <table>
            <tr>
                <th>HLA Gene</th>
                <th>Allele 1</th>
                <th>Allele 2</th>
                <th>Copy Ratio</th>
                <th>Kept Allele</th>
                <th>Lost Allele</th>
                <th>Het SNPs</th>
                <th>LOH</th>
            </tr>
"""

    for r in results:
        row_class = 'loh-yes' if r.get('LOH') == 'Y' else ''
        html += f"""
            <tr class="{row_class}">
                <td>HLA-{r.get('HLA', '')}</td>
                <td>{r.get('Allele1', '')}</td>
                <td>{r.get('Allele2', '')}</td>
                <td>{r.get('CopyRatio', '')}</td>
                <td>{r.get('KeptHLA', '')}</td>
                <td>{r.get('LostHLA', '')}</td>
                <td>{r.get('Het_num', '')}</td>
                <td>{r.get('LOH', '')}</td>
            </tr>
"""

    html += """
        </table>

        <h2>Interpretation Guide</h2>
        <ul>
            <li><strong>Copy Ratio:</strong> The ratio of allele copies (e.g., "2:0" means 2 copies of one allele, 0 of the other)</li>
            <li><strong>LOH = Y:</strong> Loss of heterozygosity detected - one allele is lost or significantly reduced</li>
            <li><strong>LOH = N:</strong> Both alleles are present at expected copy numbers</li>
            <li><strong>Het SNPs:</strong> Number of heterozygous SNPs used for the analysis (higher = more confident)</li>
        </ul>

        <h2>Clinical Significance</h2>
        <p>HLA LOH is a mechanism of immune evasion in tumors. Loss of HLA alleles can:</p>
        <ul>
            <li>Reduce tumor antigen presentation to T cells</li>
            <li>Impact response to immune checkpoint inhibitors</li>
            <li>Affect eligibility for personalized cancer vaccines</li>
        </ul>

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

    # Parse results
    results = parse_loh_results(args.loh_results)

    if not results:
        print("Warning: No LOH results found", file=sys.stderr)
        # Create empty report
        create_html_report([], args.sample, args.output_plot, args.output_report)
        return

    # Create plot
    if HAS_MATPLOTLIB:
        create_loh_plot(results, args.sample, args.output_plot)
    else:
        print("Warning: matplotlib not available, skipping plot", file=sys.stderr)

    # Create HTML report
    create_html_report(results, args.sample, args.output_plot, args.output_report)

    print(f"LOH visualization complete for {args.sample}")


if __name__ == '__main__':
    main()
