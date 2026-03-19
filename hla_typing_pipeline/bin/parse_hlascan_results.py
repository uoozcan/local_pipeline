#!/usr/bin/env python3
"""
Parse HLAscan per-gene output files into standard pipeline TSV format.

HLAscan writes one file per gene (e.g. HLA-A.txt) with output like:
    [Allele]    [Score & Coverage]
    # of considered types : N
    [Type 1] HLA-A*02:01:01:01   score=XXX  coverage=YYY
    [Type 2] HLA-A*01:01:01:01   score=ZZZ  coverage=WWW
    ...
    Low Frq Candidate type: ...

The top-2 [Type X] entries are the two alleles called.
"""

import argparse
import os
import re
import sys


GENES_OF_INTEREST = ['HLA-A', 'HLA-B', 'HLA-C',
                     'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1',
                     'HLA-DPA1', 'HLA-DPB1']

# Map HLA-X gene names to short names used in output (strip HLA- prefix)
GENE_SHORT = {g: g.replace('HLA-', '') for g in GENES_OF_INTEREST}


def normalize_allele(allele):
    """Reduce allele to 2-field resolution (e.g. HLA-A*02:01:01:01 -> A*02:01)."""
    allele = allele.strip()
    # Remove HLA- prefix if present
    allele = re.sub(r'^HLA-', '', allele)
    # Keep only first two fields (gene*field1:field2)
    m = re.match(r'([A-Z0-9]+\*[0-9]+:[0-9]+)', allele)
    if m:
        return m.group(1)
    return allele


def parse_gene_file(filepath):
    """
    Parse a single HLAscan gene output file.
    Returns (allele1, allele2) at 2-field resolution, or ('NA', 'NA') on failure.
    """
    if not os.path.isfile(filepath):
        return 'NA', 'NA'

    alleles = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                # Match lines like: [Type 1] HLA-A*02:01:01:01   score=...
                m = re.match(r'\[Type\s+\d+\]\s+(HLA-[A-Z0-9*:]+)', line)
                if m:
                    allele = normalize_allele(m.group(1))
                    if allele and allele != 'NA':
                        alleles.append(allele)
    except Exception as e:
        print("Warning: could not parse {}: {}".format(filepath, e), file=sys.stderr)
        return 'NA', 'NA'

    if not alleles:
        return 'NA', 'NA'
    allele1 = alleles[0]
    allele2 = alleles[1] if len(alleles) > 1 else alleles[0]
    return allele1, allele2


def main():
    parser = argparse.ArgumentParser(
        description='Parse HLAscan per-gene output files into standard TSV')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--input-dir', required=True,
                        help='Directory containing HLA-*.txt per-gene files')
    parser.add_argument('--output', required=True, help='Output TSV file path')
    args = parser.parse_args()

    results = []
    found_any = False

    for gene in GENES_OF_INTEREST:
        short = GENE_SHORT[gene]
        gene_file = os.path.join(args.input_dir, '{}.txt'.format(gene))
        a1, a2 = parse_gene_file(gene_file)
        if a1 != 'NA':
            found_any = True
        results.append((short, a1, a2))

    with open(args.output, 'w') as out:
        out.write('# HLAscan results for {}\n'.format(args.sample))
        out.write('# Tool: HLAscan v2.1 (SyntekaBio)\n')
        if not found_any:
            out.write('# HLAscan typing failed - no valid allele calls\n')
        out.write('#\n')
        out.write('Gene\tAllele1\tAllele2\tReads1\tReads2\n')
        for gene_short, a1, a2 in results:
            out.write('{}\t{}\t{}\tNA\tNA\n'.format(gene_short, a1, a2))

    n_genes = sum(1 for _, a1, _ in results if a1 != 'NA')
    n_alleles = sum(
        (1 if a1 != 'NA' else 0) + (1 if a2 not in ('NA', a1) else 0)
        for _, a1, a2 in results
    )
    print('HLAscan parsing complete: {} genes, {} alleles for {}'.format(
        n_genes, n_alleles, args.sample))


if __name__ == '__main__':
    main()
