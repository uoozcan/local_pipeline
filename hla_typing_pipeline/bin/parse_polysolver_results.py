#!/usr/bin/env python3
"""
Parse POLYSOLVER winners.hla.nofreq.txt to standard pipeline TSV format.

POLYSOLVER output format (tab-separated, no header):
    HLA-A   hla_a_03_01_01_01   hla_a_03_01_01_01
    HLA-B   hla_b_35_01         hla_b_40_01
    HLA-C   hla_c_04_01         hla_c_04_01

Standard pipeline output format:
    # POLYSOLVER results for {sample}
    Gene    Allele1    Allele2    Reads1    Reads2
    A       A*03:01    A*03:01    NA        NA
    B       B*35:01    B*40:01    NA        NA
    C       C*04:01    C*04:01    NA        NA
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description='Parse POLYSOLVER output to standard format')
    parser.add_argument('--input',  required=True, help='POLYSOLVER winners.hla.nofreq.txt file')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output file in standard pipeline format')
    return parser.parse_args()


def polysolver_to_standard(raw_allele: str) -> str:
    """
    Convert POLYSOLVER allele notation to standard HLA format.

    hla_a_03_01_01_01  →  A*03:01:01:01
    hla_b_35_01        →  B*35:01
    hla_drb1_09_01     →  DRB1*09:01
    """
    # Strip prefix and convert to upper case
    a = raw_allele.strip()
    if a.lower().startswith('hla_'):
        a = a[4:]  # remove 'hla_'
    a = a.upper()

    # Split on underscores: first token is gene, rest are numeric fields
    parts = a.split('_')
    if len(parts) < 2:
        return a  # cannot parse — return as-is

    gene = parts[0]
    fields = parts[1:]
    return f"{gene}*{':'.join(fields)}"


def main():
    args = parse_args()

    entries = []
    try:
        with open(args.input) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 3:
                    print(f"Warning: skipping malformed line: {line!r}", file=sys.stderr)
                    continue

                # Gene column: 'HLA-A' → 'A'
                gene = cols[0].replace('HLA-', '').strip()
                allele1 = polysolver_to_standard(cols[1])
                allele2 = polysolver_to_standard(cols[2])
                entries.append((gene, allele1, allele2))

    except FileNotFoundError:
        print(f"Warning: POLYSOLVER output file not found: {args.input}", file=sys.stderr)

    with open(args.output, 'w') as out:
        out.write(f"# POLYSOLVER results for {args.sample}\n")
        out.write("# Tool: POLYSOLVER (Broad Institute) — Class I, hg19\n")
        out.write("#\n")
        out.write("Gene\tAllele1\tAllele2\tReads1\tReads2\n")

        if not entries:
            print(f"Warning: no alleles parsed from {args.input}", file=sys.stderr)
        else:
            for gene, a1, a2 in entries:
                out.write(f"{gene}\t{a1}\t{a2}\tNA\tNA\n")

    total = len(entries)
    print(f"POLYSOLVER parsing complete: {total} genes typed for {args.sample}")


if __name__ == '__main__':
    main()
