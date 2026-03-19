#!/usr/bin/env python3
"""
Multi-source HLA Consensus Script (v1.0.0)

Integrates HLA typing results across multiple sequencing data sources
(WGS, WES, RNAseq, targeted) for the same patient.

Each source's per-sample consensus file (produced by consensus_voting.py)
is read and combined via weighted cross-source voting. A source's effective
weight per locus = base_source_weight × per_source_confidence.

Outputs:
  - Integrated consensus TSV  (Gene, Allele1, Allele2, Confidence, Discordant)
  - Source comparison TSV     (per-locus per-source breakdown + Integrated + Discordant)
  - Self-contained HTML report

Usage:
    multisource_consensus.py \\
        --patient   P001 \\
        --samples   P001_WGS,P001_WES,P001_RNA \\
        --seq-types WGS,WES,RNAseq \\
        --files     P001_WGS_consensus.txt,P001_WES_consensus.txt,P001_RNA_consensus.txt \\
        --source-weights 'WGS:1.0,WES:0.8,RNAseq:0.6,targeted:0.5' \\
        --resolution 2-field \\
        --output-consensus  P001_integrated_consensus.txt \\
        --output-comparison P001_source_comparison.txt \\
        --output-report     P001_integrated_report.html
"""

import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

CLASSICAL_GENES = [
    'HLA-A', 'HLA-B', 'HLA-C',
    'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1'
]

DEFAULT_WEIGHTS: Dict[str, float] = {
    'WGS': 1.0,
    'WES': 0.8,
    'RNAseq': 0.6,
    'targeted': 0.5,
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Multi-source HLA consensus integrator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument('--patient',            required=True,
                   help='Patient identifier')
    p.add_argument('--samples',            required=True,
                   help='Comma-separated sample IDs, parallel to --files')
    p.add_argument('--seq-types',          required=True,
                   help='Comma-separated seq_type labels (WGS, WES, RNAseq, targeted)')
    p.add_argument('--files',              required=True,
                   help='Comma-separated paths to per-source consensus files')
    p.add_argument('--source-weights',
                   default='WGS:1.0,WES:0.8,RNAseq:0.6,targeted:0.5',
                   help='Source weights as KEY:VALUE pairs, comma-separated')
    p.add_argument('--resolution',         default='2-field',
                   choices=['2-field', '4-field'])
    p.add_argument('--output-consensus',   required=True,
                   help='Path for integrated consensus output TSV')
    p.add_argument('--output-comparison',  required=True,
                   help='Path for source comparison output TSV')
    p.add_argument('--output-report',      required=True,
                   help='Path for HTML report output')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_source_weights(weights_str: str) -> Dict[str, float]:
    """Parse 'WGS:1.0,WES:0.8,...' into dict, filling missing keys from defaults."""
    weights = dict(DEFAULT_WEIGHTS)
    for part in weights_str.split(','):
        part = part.strip()
        if ':' in part:
            key, val = part.split(':', 1)
            try:
                weights[key.strip()] = float(val.strip())
            except ValueError:
                print(f"Warning: could not parse weight '{part}', using default",
                      file=sys.stderr)
    return weights


def parse_consensus_file(filepath: str) -> Dict[str, Dict]:
    """
    Parse a *_consensus.txt file produced by consensus_voting.py.

    Expected TSV format (comment lines start with #):
        Gene  Allele1  Allele2  Confidence  [Reads1  Reads2]

    Returns dict: {gene: {allele1, allele2, confidence, reads1, reads2}}
    """
    results: Dict[str, Dict] = {}
    try:
        with open(filepath, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#') or line.lower().startswith('gene'):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                gene   = parts[0].strip()
                a1     = parts[1].strip() if parts[1].strip() not in ('', '-', 'NA') else None
                a2     = parts[2].strip() if parts[2].strip() not in ('', '-', 'NA') else None
                try:
                    conf = float(parts[3])
                except (ValueError, IndexError):
                    conf = 0.0
                try:
                    reads1 = int(parts[4]) if len(parts) > 4 and parts[4].strip().isdigit() else 0
                except (ValueError, IndexError):
                    reads1 = 0
                try:
                    reads2 = int(parts[5]) if len(parts) > 5 and parts[5].strip().isdigit() else 0
                except (ValueError, IndexError):
                    reads2 = 0
                results[gene] = {
                    'allele1': a1, 'allele2': a2,
                    'confidence': conf,
                    'reads1': reads1, 'reads2': reads2,
                }
    except FileNotFoundError:
        print(f"Warning: consensus file not found: {filepath}", file=sys.stderr)
    except Exception as exc:
        print(f"Warning: error reading {filepath}: {exc}", file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Core integration algorithm
# ---------------------------------------------------------------------------

def integrate_across_sources(
    sources: List[Dict],
    gene: str,
) -> Tuple[Optional[str], Optional[str], float, Dict]:
    """
    Weighted cross-source vote for a single HLA locus.

    For each source with allele calls:
        effective_weight = base_source_weight × per_source_confidence

    Evidence is accumulated per allele string. The top two alleles by
    cumulative evidence become the integrated call.

    Returns (allele1, allele2, integrated_confidence, details_dict)
    """
    allele_evidence: Dict[str, float] = defaultdict(float)
    allele_source_map: Dict[str, List[str]] = defaultdict(list)
    per_source_calls: Dict[str, List] = {}

    total_possible_weight = 0.0

    for src in sources:
        sid   = src['sample_id']
        base_w = src['weight']
        locus  = src['locus_data'].get(gene, {})

        a1   = locus.get('allele1')
        a2   = locus.get('allele2')
        conf = locus.get('confidence', 0.0)

        per_source_calls[sid] = [a1, a2]

        if not a1 and not a2:
            continue  # source did not type this locus

        effective_w = base_w * conf
        total_possible_weight += base_w  # count source even if confidence is low

        for allele in [a1, a2]:
            if allele:
                allele_evidence[allele] += effective_w
                allele_source_map[allele].append(f"{sid}({src['seq_type']})")

    if not allele_evidence:
        return None, None, 0.0, {
            'per_source_calls': per_source_calls,
            'allele_evidence': {},
            'discordant': False,
        }

    sorted_alleles = sorted(allele_evidence.items(), key=lambda x: x[1], reverse=True)
    integrated1 = sorted_alleles[0][0] if len(sorted_alleles) >= 1 else None
    integrated2 = sorted_alleles[1][0] if len(sorted_alleles) >= 2 else None

    # Integrated confidence: evidence for top pair / max possible evidence
    top_evidence = sum(e for _, e in sorted_alleles[:2])
    max_possible = 2.0 * total_possible_weight if total_possible_weight > 0 else 1.0
    integrated_confidence = min(1.0, top_evidence / max_possible)

    # Discordance: any two sources that both typed this locus but called different allele sets
    typed_source_allele_sets = []
    for src in sources:
        locus = src['locus_data'].get(gene, {})
        a1, a2 = locus.get('allele1'), locus.get('allele2')
        if a1 or a2:
            typed_source_allele_sets.append(frozenset(x for x in [a1, a2] if x))

    discordant = len(set(typed_source_allele_sets)) > 1 if len(typed_source_allele_sets) > 1 else False

    details = {
        'per_source_calls': per_source_calls,
        'allele_evidence': dict(allele_evidence),
        'allele_sources': {k: v for k, v in allele_source_map.items()},
        'discordant': discordant,
    }
    return integrated1, integrated2, integrated_confidence, details


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_integrated_consensus(
    path: str,
    patient_id: str,
    sources: List[Dict],
    ordered_genes: List[str],
    integrated: Dict[str, Dict],
    resolution: str,
    source_weights_str: str,
) -> None:
    with open(path, 'w') as fh:
        src_labels = ', '.join(
            f"{s['sample_id']}({s['seq_type']}, w={s['weight']:.2f})" for s in sources
        )
        fh.write(f"# Multi-source Integrated HLA Consensus\n")
        fh.write(f"# Patient:    {patient_id}\n")
        fh.write(f"# Sources:    {src_labels}\n")
        fh.write(f"# Weights:    {source_weights_str}\n")
        fh.write(f"# Resolution: {resolution}\n")
        fh.write(f"# Algorithm:  effective_weight = source_weight x per_source_confidence\n#\n")
        fh.write("Gene\tAllele1\tAllele2\tConfidence\tDiscordant\n")
        for gene in ordered_genes:
            res = integrated[gene]
            disc = 'YES' if res['details'].get('discordant') else 'NO'
            fh.write(
                f"{gene}\t"
                f"{res['allele1'] or '-'}\t"
                f"{res['allele2'] or '-'}\t"
                f"{res['confidence']:.3f}\t"
                f"{disc}\n"
            )


def write_source_comparison(
    path: str,
    patient_id: str,
    sources: List[Dict],
    ordered_genes: List[str],
    integrated: Dict[str, Dict],
) -> None:
    with open(path, 'w') as fh:
        fh.write(f"# Multi-source Source Comparison\n")
        fh.write(f"# Patient: {patient_id}\n#\n")
        src_headers = '\t'.join(
            f"{s['sample_id']}({s['seq_type']})" for s in sources
        )
        fh.write(f"Gene\t{src_headers}\tIntegrated\tDiscordant\n")
        for gene in ordered_genes:
            res = integrated[gene]
            per_src = res['details']['per_source_calls']
            row = [gene]
            for src in sources:
                calls = per_src.get(src['sample_id'], [None, None])
                row.append(f"{calls[0] or '-'}/{calls[1] or '-'}")
            integrated_call = f"{res['allele1'] or '-'}/{res['allele2'] or '-'}"
            row.append(integrated_call)
            row.append('YES' if res['details'].get('discordant') else 'NO')
            fh.write('\t'.join(row) + '\n')


def write_html_report(
    path: str,
    patient_id: str,
    sources: List[Dict],
    ordered_genes: List[str],
    integrated: Dict[str, Dict],
) -> None:
    n_discordant = sum(1 for g in ordered_genes if integrated[g]['details'].get('discordant'))
    n_typed      = sum(1 for g in ordered_genes if integrated[g]['allele1'])
    avg_conf     = (
        sum(integrated[g]['confidence'] for g in ordered_genes if integrated[g]['allele1']) / n_typed
        if n_typed else 0.0
    )

    # Per-source typed locus count
    src_stats_html = ''
    for src in sources:
        n_src = sum(
            1 for g in ordered_genes
            if any(integrated[g]['details']['per_source_calls'].get(src['sample_id'], [None, None]))
        )
        src_stats_html += (
            f"<tr><td>{src['sample_id']}</td><td>{src['seq_type']}</td>"
            f"<td>{src['weight']:.2f}</td><td>{n_src}/{len(ordered_genes)}</td></tr>"
        )

    # Result table rows
    src_header_cells = ''.join(
        f"<th>{s['sample_id']}<br><small>{s['seq_type']}</small></th>" for s in sources
    )
    table_rows = ''
    for gene in ordered_genes:
        res  = integrated[gene]
        disc = res['details'].get('discordant', False)
        bg   = '#d4edda' if not disc else '#f8d7da'
        conf_pct = int(res['confidence'] * 100)
        conf_color = '#28a745' if conf_pct >= 80 else ('#fd7e14' if conf_pct >= 50 else '#dc3545')
        cells = ''
        for src in sources:
            calls = res['details']['per_source_calls'].get(src['sample_id'], [None, None])
            cells += f"<td>{calls[0] or '-'} / {calls[1] or '-'}</td>"
        table_rows += (
            f"<tr style='background:{bg}'>"
            f"<td><strong>{gene}</strong></td>"
            f"<td>{res['allele1'] or '-'}</td>"
            f"<td>{res['allele2'] or '-'}</td>"
            f"<td><span style='color:{conf_color};font-weight:bold'>{res['confidence']:.3f}</span></td>"
            f"{cells}"
            f"<td>{'&#9888; YES' if disc else 'NO'}</td>"
            f"</tr>"
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Multi-source HLA Report &mdash; {patient_id}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 2em; color: #333; }}
  h1 {{ color: #343a40; border-bottom: 2px solid #6c757d; padding-bottom: 0.3em; }}
  h2 {{ color: #495057; margin-top: 1.5em; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 1.5em; }}
  th {{ background: #343a40; color: white; padding: 8px; text-align: left; }}
  td {{ padding: 6px 8px; border-bottom: 1px solid #dee2e6; }}
  .card {{ background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px;
           padding: 1em; display: inline-block; margin: 0.5em; min-width: 140px; text-align: center; }}
  .card .num {{ font-size: 2em; font-weight: bold; }}
  .note {{ font-size: 0.85em; color: #6c757d; }}
</style>
</head>
<body>
<h1>Multi-source HLA Consensus Report</h1>
<p class="note">Pipeline v1.3.0 &mdash; Patient: <strong>{patient_id}</strong></p>

<h2>Summary</h2>
<div>
  <div class="card"><div class="num">{len(sources)}</div>Sources integrated</div>
  <div class="card"><div class="num">{n_typed}</div>Loci typed</div>
  <div class="card"><div class="num" style="color:{'#28a745' if avg_conf>=0.7 else '#fd7e14'}">{avg_conf:.2f}</div>Avg confidence</div>
  <div class="card"><div class="num" style="color:{'#dc3545' if n_discordant>0 else '#28a745'}">{n_discordant}</div>Discordant loci</div>
</div>

<h2>Sources</h2>
<table>
  <tr><th>Sample ID</th><th>Seq type</th><th>Base weight</th><th>Loci typed</th></tr>
  {src_stats_html}
</table>

<h2>Integrated Results</h2>
<p class="note">
  Green rows = concordant across all sources that typed the locus.<br>
  Red rows = discordant allele calls detected between sources (review recommended).
</p>
<table>
  <tr>
    <th>Gene</th><th>Allele 1</th><th>Allele 2</th><th>Confidence</th>
    {src_header_cells}
    <th>Discordant</th>
  </tr>
  {table_rows}
</table>

<p class="note">
  <strong>Confidence</strong>: min(1.0, sum of top-2 allele evidence / max possible evidence).
  Effective weight per source per locus = base_source_weight &times; per_source_confidence.
</p>
</body>
</html>
"""
    with open(path, 'w') as fh:
        fh.write(html)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    sample_ids = [s.strip() for s in args.samples.split(',')]
    seq_types  = [s.strip() for s in args.seq_types.split(',')]
    file_paths = [s.strip() for s in args.files.split(',')]

    if not (len(sample_ids) == len(seq_types) == len(file_paths)):
        print(
            f"Error: --samples ({len(sample_ids)}), --seq-types ({len(seq_types)}), "
            f"and --files ({len(file_paths)}) must have equal length.",
            file=sys.stderr
        )
        sys.exit(1)

    source_weights = parse_source_weights(args.source_weights)

    # Build source descriptors
    sources: List[Dict] = []
    for sid, stype, fpath in zip(sample_ids, seq_types, file_paths):
        locus_data = parse_consensus_file(fpath)
        if not locus_data:
            print(f"Warning: no data parsed from {fpath} (sample {sid})", file=sys.stderr)
        sources.append({
            'sample_id':  sid,
            'seq_type':   stype,
            'weight':     source_weights.get(stype, 0.5),
            'locus_data': locus_data,
        })

    # Collect all gene names across sources, preserving classical gene order
    all_genes: set = set()
    for src in sources:
        all_genes.update(src['locus_data'].keys())
    ordered_genes = [g for g in CLASSICAL_GENES if g in all_genes]
    ordered_genes += sorted(g for g in all_genes if g not in CLASSICAL_GENES)

    if not ordered_genes:
        print("Warning: no typed genes found across any source.", file=sys.stderr)

    # Integrate per gene
    integrated: Dict[str, Dict] = {}
    for gene in ordered_genes:
        a1, a2, conf, details = integrate_across_sources(sources, gene)
        integrated[gene] = {
            'allele1': a1, 'allele2': a2,
            'confidence': conf, 'details': details,
        }

    # Write outputs
    write_integrated_consensus(
        args.output_consensus, args.patient, sources,
        ordered_genes, integrated, args.resolution, args.source_weights
    )
    write_source_comparison(
        args.output_comparison, args.patient, sources, ordered_genes, integrated
    )
    write_html_report(
        args.output_report, args.patient, sources, ordered_genes, integrated
    )

    n_disc = sum(1 for g in ordered_genes if integrated[g]['details'].get('discordant'))
    print(f"Patient:            {args.patient}")
    print(f"Sources integrated: {len(sources)}")
    print(f"Loci typed:         {len(ordered_genes)}")
    print(f"Discordant loci:    {n_disc}/{len(ordered_genes)}")
    print(f"Integrated consensus -> {args.output_consensus}")
    print(f"Source comparison   -> {args.output_comparison}")
    print(f"HTML report         -> {args.output_report}")


if __name__ == '__main__':
    main()
