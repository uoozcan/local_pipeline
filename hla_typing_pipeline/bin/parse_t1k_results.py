#!/usr/bin/env python3
"""
parse_t1k_results.py — Convert T1K genotype TSV to standard pipeline format.

T1K genotype TSV format (tab-separated):
    #gene   allele1     abundance1  allele2     abundance2
    HLA-A   A*01:01     42.00       A*24:02     38.00
    HLA-B   B*07:02     35.00       B*44:02     30.00

Standard output format (matches other tools in the pipeline):
    # T1K results for <sample>
    Gene    Allele1     Allele2
    HLA-A   A*01:01     A*24:02
    HLA-B   B*07:02     B*44:02
"""

import sys
from pathlib import Path


def parse_t1k_genotype(input_path: str, output_path: str) -> None:
    results = []
    try:
        with open(input_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"
                allele1 = parts[1].strip() if len(parts) > 1 else '-'
                allele2 = parts[3].strip() if len(parts) > 3 else '-'
                # Normalize blanks
                allele1 = allele1 if allele1 not in ('', '.', 'NA', 'None') else '-'
                allele2 = allele2 if allele2 not in ('', '.', 'NA', 'None') else '-'
                results.append((gene, allele1, allele2))
    except Exception as e:
        print(f"Warning: Error parsing T1K file {input_path}: {e}", file=sys.stderr)

    sample_name = Path(input_path).stem.replace('_genotype', '')
    with open(output_path, 'w') as out:
        out.write(f"# T1K results for {sample_name}\n")
        out.write("Gene\tAllele1\tAllele2\n")
        for gene, a1, a2 in results:
            out.write(f"{gene}\t{a1}\t{a2}\n")
        if not results:
            out.write("# No alleles called\n")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <t1k_genotype.tsv> <output.txt>", file=sys.stderr)
        sys.exit(1)
    parse_t1k_genotype(sys.argv[1], sys.argv[2])
