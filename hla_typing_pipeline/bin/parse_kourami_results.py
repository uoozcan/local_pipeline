#!/usr/bin/env python3
"""
Parse Kourami .result file to standard pipeline TSV format.

Kourami .result format (tab-separated, one allele per line, no header):
    A*11:01:01G    546    1.0    546    546    29.0    8.0    14.0
    A*03:01:01G    532    1.0    532    532    28.5    7.5    13.5
    B*35:01:01G    498    0.99   498    495    25.0    6.0    12.0
    ...

Columns: allele, matched_bases, identity, assembled_length, matched_length,
         bottleneck_sum, bottleneck_1, bottleneck_2

The 'G' suffix denotes a G-group allele designation — stripped before output.
Lines are ordered by gene then descending matched_bases; the top 2 per gene
form the diploid call.

Standard pipeline output format:
    # Kourami results for {sample}
    Gene    Allele1    Allele2    Reads1    Reads2
    A       A*11:01    A*03:01    546       532
"""

import argparse
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Parse Kourami output to standard format')
    parser.add_argument('--input',  required=True, help='Kourami .result file')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output file in standard pipeline format')
    return parser.parse_args()


def strip_g_group(allele: str) -> str:
    """Strip trailing G-group suffix (e.g. A*11:01:01G → A*11:01:01)."""
    return allele.rstrip('G') if allele.endswith('G') else allele


def main():
    args = parse_args()

    # gene → list of (allele, matched_bases) sorted by original file order
    # Kourami already outputs alleles sorted by confidence; keep insertion order.
    gene_alleles: dict = defaultdict(list)

    try:
        with open(args.input) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                cols = line.split('\t')
                if len(cols) < 1:
                    continue

                raw_allele = cols[0].strip()
                matched_bases = 0
                if len(cols) >= 2:
                    try:
                        matched_bases = int(cols[1])
                    except ValueError:
                        pass

                allele = strip_g_group(raw_allele)

                # Extract gene name (everything before the '*')
                if '*' not in allele:
                    print(f"Warning: cannot parse allele '{allele}', skipping", file=sys.stderr)
                    continue
                gene = allele.split('*')[0].upper()

                gene_alleles[gene].append((allele, matched_bases))

    except FileNotFoundError:
        print(f"Warning: Kourami result file not found: {args.input}", file=sys.stderr)

    with open(args.output, 'w') as out:
        out.write(f"# Kourami results for {args.sample}\n")
        out.write("# Tool: Kourami v0.9.6 -- Class I, hg38/hs38NoAltDH\n")
        out.write("#\n")
        out.write("Gene\tAllele1\tAllele2\tReads1\tReads2\n")

        if not gene_alleles:
            print(f"Warning: no alleles parsed from {args.input}", file=sys.stderr)
        else:
            for gene in sorted(gene_alleles.keys()):
                alleles = gene_alleles[gene]
                a1, r1 = alleles[0] if len(alleles) > 0 else ('NA', 0)
                a2, r2 = alleles[1] if len(alleles) > 1 else (a1, r1)   # homozygous fallback
                out.write(f"{gene}\t{a1}\t{a2}\t{r1}\t{r2}\n")

    total = len(gene_alleles)
    print(f"Kourami parsing complete: {total} genes typed for {args.sample}")


if __name__ == '__main__':
    main()
