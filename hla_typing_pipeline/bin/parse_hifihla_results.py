#!/usr/bin/env python3
"""
parse_hifihla_results.py — Convert HiFi-HLA output to standard pipeline TSV.

HiFi-HLA output format (*_calls.tsv):
    #gene    allele1            confidence1    allele2            confidence2
    HLA-A    A*03:01:01:01      0.99           A*11:01:01:01      0.97
    HLA-B    B*07:02:01:01      0.98           B*35:01:01:02      0.95
    ...

Standard output format (matches all other tools in the pipeline):
    # HiFi-HLA results for <sample>
    Gene        Allele1            Allele2
    HLA-A       A*03:01:01:01      A*11:01:01:01
    HLA-B       B*07:02:01:01      B*35:01:01:02

Confidence scores are preserved in columns 3 and 4 so that
consensus_voting.py can use them directly as quality weights.

Usage:
    parse_hifihla_results.py <input_calls.tsv> <output.txt>
"""

import sys
from pathlib import Path


# Classical HLA genes of interest (non-classical are silently dropped)
CLASSICAL_GENES = {
    'HLA-A', 'HLA-B', 'HLA-C',
    'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
    'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1',
}


def _gene_key(gene: str) -> str:
    """Normalise gene name to HLA-prefixed form."""
    gene = gene.strip().lstrip('#')
    return gene if gene.startswith('HLA-') else f"HLA-{gene}"


def parse_hifihla_calls(input_path: str):
    """
    Parse a HiFi-HLA *_calls.tsv file.

    Expected columns (tab-separated, may have leading '#' on header line):
        gene  allele1  confidence1  allele2  confidence2

    Returns list of (gene, allele1, conf1, allele2, conf2).
    """
    results = []
    header_seen = False

    try:
        with open(input_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                # Detect header by presence of 'allele' or 'gene' (case-insensitive)
                lower = line.lower().lstrip('#').strip()
                if not header_seen and ('allele' in lower or lower.startswith('gene')):
                    header_seen = True
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                gene = _gene_key(parts[0])

                # allele1
                a1 = parts[1].strip() if len(parts) > 1 else ''
                a1 = a1 if a1 not in ('', '.', '-', 'NA', 'None') else ''

                # confidence1 (float 0–1)
                try:
                    c1 = float(parts[2]) if len(parts) > 2 else 0.9
                    c1 = min(1.0, max(0.0, c1))
                except ValueError:
                    c1 = 0.9

                # allele2 (may be absent → homozygous)
                a2 = parts[3].strip() if len(parts) > 3 else ''
                a2 = a2 if a2 not in ('', '.', '-', 'NA', 'None') else ''

                # confidence2
                try:
                    c2 = float(parts[4]) if len(parts) > 4 else c1
                    c2 = min(1.0, max(0.0, c2))
                except (ValueError, IndexError):
                    c2 = c1

                if a1:
                    results.append((gene, a1, c1, a2, c2))

    except OSError as exc:
        print(f"Error reading HiFi-HLA file {input_path}: {exc}", file=sys.stderr)
        sys.exit(1)

    return results


def write_standard_tsv(results, output_path: str, sample_name: str) -> None:
    """
    Write results in the standard pipeline TSV format.

    Extra columns (Conf1, Conf2) are preserved so that consensus_voting.py
    can read them as quality scores (analogous to read counts for other tools).
    """
    with open(output_path, 'w') as out:
        out.write(f"# HiFi-HLA results for {sample_name}\n")
        out.write("Gene\tAllele1\tAllele2\tConf1\tConf2\n")
        for gene, a1, c1, a2, c2 in results:
            a2_out = a2 if a2 else '-'
            out.write(f"{gene}\t{a1}\t{a2_out}\t{c1:.4f}\t{c2:.4f}\n")
        if not results:
            out.write("# No alleles called\n")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <calls.tsv> <output.txt>", file=sys.stderr)
        sys.exit(1)

    input_path  = sys.argv[1]
    output_path = sys.argv[2]
    sample_name = Path(input_path).stem.replace('_calls', '').replace('_hla_calls', '')

    results = parse_hifihla_calls(input_path)

    # Keep only classical genes
    results = [(g, a1, c1, a2, c2) for g, a1, c1, a2, c2 in results
               if g in CLASSICAL_GENES]

    write_standard_tsv(results, output_path, sample_name)
    print(f"Parsed {len(results)} loci from {Path(input_path).name} → {output_path}")


if __name__ == '__main__':
    main()
