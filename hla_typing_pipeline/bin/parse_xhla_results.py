#!/usr/bin/env python3
"""
Parse xHLA JSON output to standard HLA typing format.

xHLA outputs results in JSON format like:
{
    "sample_id": "sample",
    "hla": {
        "alleles": [
            "A*02:01",
            "A*03:01",
            "B*07:02",
            ...
        ]
    }
}

This script converts to the standard tab-separated format:
Gene    Allele1    Allele2    Reads1    Reads2
"""

import argparse
import json
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Parse xHLA results to standard format')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--input', required=True, help='Input JSON file from xHLA')
    parser.add_argument('--output', required=True, help='Output file in standard format')
    return parser.parse_args()


def parse_xhla_json(filepath):
    """Parse xHLA JSON output file."""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        print(f"Warning: xHLA output file not found: {filepath}", file=sys.stderr)
        return None
    except json.JSONDecodeError as e:
        print(f"Warning: Error parsing xHLA JSON: {e}", file=sys.stderr)
        return None


def extract_alleles(data):
    """Extract alleles from xHLA JSON data and organize by gene."""
    if not data:
        return {}

    alleles_by_gene = defaultdict(list)

    # Get alleles from the JSON structure
    hla_data = data.get('hla', {})
    alleles = hla_data.get('alleles', [])

    # Also check alternative JSON structures
    if not alleles and 'alleles' in data:
        alleles = data['alleles']

    for allele in alleles:
        if not allele or allele == 'NA':
            continue

        # Parse allele format (e.g., "A*02:01" or "HLA-A*02:01")
        allele = allele.replace('HLA-', '')

        if '*' in allele:
            gene, rest = allele.split('*', 1)
            gene = gene.upper()

            # Normalize gene name
            if not gene.startswith('HLA-'):
                gene_normalized = f"HLA-{gene}"
            else:
                gene_normalized = gene

            alleles_by_gene[gene_normalized].append(allele)

    return alleles_by_gene


def write_standard_format(sample, alleles_by_gene, output_path):
    """Write results in standard pipeline format."""

    # Standard HLA genes in order
    standard_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1']

    with open(output_path, 'w') as f:
        # Write header
        f.write(f"# xHLA results for {sample}\n")
        f.write("# Tool: xHLA (Human Longevity Inc.)\n")
        f.write("#\n")
        f.write("Gene\tAllele1\tAllele2\tReads1\tReads2\n")

        # Write results for each gene
        for gene in standard_genes:
            alleles = alleles_by_gene.get(gene, [])

            if len(alleles) >= 2:
                allele1 = alleles[0]
                allele2 = alleles[1]
            elif len(alleles) == 1:
                allele1 = alleles[0]
                allele2 = alleles[0]  # Homozygous
            else:
                allele1 = 'NA'
                allele2 = 'NA'

            # xHLA doesn't provide read counts, use placeholder
            reads1 = 'NA'
            reads2 = 'NA'

            # Extract short gene name for output
            gene_short = gene.replace('HLA-', '')

            f.write(f"{gene_short}\t{allele1}\t{allele2}\t{reads1}\t{reads2}\n")

        # Also write any additional genes found that aren't in standard list
        for gene, alleles in alleles_by_gene.items():
            if gene not in standard_genes and alleles:
                gene_short = gene.replace('HLA-', '')
                allele1 = alleles[0] if alleles else 'NA'
                allele2 = alleles[1] if len(alleles) > 1 else allele1
                f.write(f"{gene_short}\t{allele1}\t{allele2}\tNA\tNA\n")


def main():
    args = parse_args()

    # Parse xHLA JSON
    data = parse_xhla_json(args.input)

    if data is None:
        # Create empty output file
        with open(args.output, 'w') as f:
            f.write(f"# xHLA results for {args.sample}\n")
            f.write("# ERROR: Could not parse xHLA output\n")
            f.write("Gene\tAllele1\tAllele2\tReads1\tReads2\n")
        print(f"Warning: Created empty output for {args.sample}", file=sys.stderr)
        return

    # Extract alleles by gene
    alleles_by_gene = extract_alleles(data)

    # Write standard format
    write_standard_format(args.sample, alleles_by_gene, args.output)

    # Print summary
    total_alleles = sum(len(v) for v in alleles_by_gene.values())
    genes_typed = len([g for g, a in alleles_by_gene.items() if a])
    print(f"xHLA parsing complete: {genes_typed} genes, {total_alleles} alleles for {args.sample}")


if __name__ == '__main__':
    main()
