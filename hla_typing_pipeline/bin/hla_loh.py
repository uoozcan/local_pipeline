#!/usr/bin/env python3
"""
HLA Loss of Heterozygosity (LOH) Calculator

Calculates HLA allele copy numbers and detects LOH events
using tumor purity and ploidy estimates.

Based on SpecHLA's cal.hla.copy.pl algorithm.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class AlleleFrequency:
    """Store allele frequency data."""
    allele1_freq: float
    allele2_freq: float
    het_snp_count: int


@dataclass
class LOHResult:
    """Store LOH analysis result for one HLA gene."""
    sample: str
    hla_gene: str
    allele1: str
    allele2: str
    copy_ratio: str
    kept_hla: str
    lost_hla: str
    freq1: float
    freq2: float
    purity: float
    het_num: int
    loh: str  # 'Y' or 'N'
    copy_a: int
    copy_b: int


def parse_args():
    parser = argparse.ArgumentParser(
        description='Calculate HLA copy numbers and detect LOH'
    )
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--purity', type=float, required=True,
                        help='Tumor purity (0-1)')
    parser.add_argument('--ploidy', type=float, required=True,
                        help='Tumor ploidy')
    parser.add_argument('--freq-list', required=True,
                        help='File containing list of frequency files')
    parser.add_argument('--typing-file', required=True,
                        help='HLA typing result file')
    parser.add_argument('--het-cutoff', type=int, default=5,
                        help='Minimum heterozygous SNPs for LOH call (default: 5)')
    parser.add_argument('--output', required=True,
                        help='Output LOH results file')
    parser.add_argument('--summary', required=True,
                        help='Output summary file')
    return parser.parse_args()


def parse_typing_file(filepath: str) -> Dict[str, Tuple[str, str]]:
    """
    Parse SpecHLA typing result file.

    Returns dict mapping HLA gene to (allele1, allele2) tuple.
    """
    alleles = {}

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            return alleles

        # Parse header and allele lines
        header = lines[0].strip().split('\t')
        allele_line = lines[1].strip().split('\t') if len(lines) > 1 else []

        for i, col in enumerate(header):
            if i == 0:
                continue  # Skip sample column

            # Extract HLA gene name
            if '_' in col:
                gene = col.split('_')[1]  # e.g., "sample_A" -> "A"
            else:
                gene = col.replace('HLA-', '').replace('HLA_', '')

            if i < len(allele_line):
                allele = allele_line[i]
                if gene not in alleles:
                    alleles[gene] = (allele, '-')
                else:
                    alleles[gene] = (alleles[gene][0], allele)

    except Exception as e:
        print(f"Warning: Error parsing typing file: {e}", file=sys.stderr)

    return alleles


def parse_frequency_file(filepath: str) -> Optional[AlleleFrequency]:
    """
    Parse SpecHLA frequency file.

    Format:
    Allele  Frequency
    A*02:01 0.65
    A*03:01 0.35
    Het_SNP_num: 12
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        if len(lines) < 4:
            return None

        # Skip header, get first two allele frequencies
        freq1 = float(lines[1].strip().split()[1]) if len(lines) > 1 else 0.0
        freq2 = float(lines[2].strip().split()[1]) if len(lines) > 2 else 0.0

        # Parse heterozygous SNP count
        het_line = lines[3].strip() if len(lines) > 3 else "Het_SNP_num: 0"
        het_count = int(het_line.split(':')[-1].strip()) if ':' in het_line else 0

        return AlleleFrequency(
            allele1_freq=freq1,
            allele2_freq=freq2,
            het_snp_count=het_count
        )

    except Exception as e:
        print(f"Warning: Error parsing frequency file {filepath}: {e}", file=sys.stderr)
        return None


def calculate_copy_numbers(
    freq1: float,
    freq2: float,
    purity: float,
    ploidy: float
) -> Tuple[int, int, int]:
    """
    Calculate allele copy numbers from frequencies.

    Returns (copy_a, copy_b, tag) where:
    - copy_a, copy_b are the copy numbers
    - tag indicates which allele has higher frequency (0=equal, 1=freq2 higher, 2=freq1 higher)
    """
    if freq1 == 0:
        return (0, round(ploidy), 0)

    if freq2 == 0:
        return (round(ploidy), 0, 0)

    # Calculate frequency ratio
    if freq1 > freq2:
        freq_ratio = freq1 / freq2
        tag = 2
    else:
        freq_ratio = freq2 / freq1
        tag = 1

    # Calculate copy numbers using purity/ploidy model
    # Formula derived from tumor content modeling
    copy_a = (freq_ratio * ploidy + (freq_ratio - 1) * (1/purity - 1)) / (1 + freq_ratio)
    copy_b = ploidy - copy_a

    if copy_b < 0:
        copy_b = 0
        copy_a = ploidy

    copy_a = round(copy_a)
    copy_b = round(copy_b)

    return (copy_a, copy_b, tag)


def detect_loh(
    sample: str,
    hla_gene: str,
    allele1: str,
    allele2: str,
    freq_data: AlleleFrequency,
    purity: float,
    ploidy: float,
    het_cutoff: int
) -> LOHResult:
    """
    Detect LOH for a single HLA gene.
    """
    freq1 = freq_data.allele1_freq
    freq2 = freq_data.allele2_freq
    het_num = freq_data.het_snp_count

    copy_a, copy_b, tag = calculate_copy_numbers(freq1, freq2, purity, ploidy)

    # Determine copy ratio string and which allele is kept/lost
    if tag == 0:
        if freq1 == 0:
            copy_ratio = f"0:{round(ploidy)}"
            kept_hla = allele1
            lost_hla = "homozygous"
        else:
            copy_ratio = f"{round(ploidy)}:0"
            kept_hla = allele1
            lost_hla = "homozygous"
    elif tag == 1:
        copy_ratio = f"{copy_b}:{copy_a}"
        kept_hla = allele2
        lost_hla = allele1
    else:
        copy_ratio = f"{copy_a}:{copy_b}"
        kept_hla = allele1
        lost_hla = allele2

    # Determine LOH status
    # LOH is called when copy numbers are unequal AND sufficient heterozygous SNPs
    if copy_a == copy_b or het_num < het_cutoff:
        loh = 'N'
    else:
        loh = 'Y'

    # For homozygous cases, no LOH
    if tag == 0:
        loh = 'N'
        lost_hla = "homozygous"

    return LOHResult(
        sample=sample,
        hla_gene=hla_gene,
        allele1=allele1,
        allele2=allele2,
        copy_ratio=copy_ratio,
        kept_hla=kept_hla,
        lost_hla=lost_hla,
        freq1=freq1,
        freq2=freq2,
        purity=purity,
        het_num=het_num,
        loh=loh,
        copy_a=copy_a,
        copy_b=copy_b
    )


def main():
    args = parse_args()

    # Validate inputs
    if not 0 < args.purity <= 1:
        print(f"Error: Purity must be between 0 and 1, got {args.purity}", file=sys.stderr)
        sys.exit(1)

    if args.ploidy <= 0:
        print(f"Error: Ploidy must be positive, got {args.ploidy}", file=sys.stderr)
        sys.exit(1)

    # Parse typing results
    alleles = parse_typing_file(args.typing_file)

    if not alleles:
        print("Warning: No HLA alleles found in typing file", file=sys.stderr)

    # Parse frequency file list
    freq_files = {}
    try:
        with open(args.freq_list, 'r') as f:
            for line in f:
                filepath = line.strip()
                if filepath:
                    # Extract HLA gene from filename (e.g., HLA_A_freq.txt -> A)
                    filename = Path(filepath).name
                    gene = filename.replace('HLA_', '').replace('_freq.txt', '')
                    freq_files[gene] = filepath
    except Exception as e:
        print(f"Error reading frequency list: {e}", file=sys.stderr)
        sys.exit(1)

    # Process each HLA gene
    results: List[LOHResult] = []
    hla_genes = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']

    for gene in hla_genes:
        if gene not in freq_files:
            continue

        freq_data = parse_frequency_file(freq_files[gene])
        if freq_data is None:
            continue

        allele1, allele2 = alleles.get(gene, ('-', '-'))

        result = detect_loh(
            sample=args.sample,
            hla_gene=gene,
            allele1=allele1,
            allele2=allele2,
            freq_data=freq_data,
            purity=args.purity,
            ploidy=args.ploidy,
            het_cutoff=args.het_cutoff
        )
        results.append(result)

    # Write results
    with open(args.output, 'w') as f:
        # Header
        f.write("Sample\tHLA\tAllele1\tAllele2\tCopyRatio\tKeptHLA\tLostHLA\t")
        f.write("Freq1\tFreq2\tPurity\tHet_num\tLOH\n")

        for r in results:
            f.write(f"{r.sample}\t{r.hla_gene}\t{r.allele1}\t{r.allele2}\t")
            f.write(f"{r.copy_ratio}\t{r.kept_hla}\t{r.lost_hla}\t")
            f.write(f"{r.freq1:.4f}\t{r.freq2:.4f}\t{r.purity:.2f}\t{r.het_num}\t{r.loh}\n")

    # Write summary
    loh_count = sum(1 for r in results if r.loh == 'Y')
    total_genes = len(results)

    with open(args.summary, 'w') as f:
        f.write(f"HLA LOH Analysis Summary for {args.sample}\n")
        f.write(f"{'='*50}\n\n")
        f.write(f"Parameters:\n")
        f.write(f"  Tumor Purity: {args.purity:.2f}\n")
        f.write(f"  Tumor Ploidy: {args.ploidy:.1f}\n")
        f.write(f"  Het SNP Cutoff: {args.het_cutoff}\n\n")
        f.write(f"Results:\n")
        f.write(f"  HLA genes analyzed: {total_genes}\n")
        f.write(f"  LOH events detected: {loh_count}\n\n")

        if loh_count > 0:
            f.write("LOH Events:\n")
            for r in results:
                if r.loh == 'Y':
                    f.write(f"  HLA-{r.hla_gene}: Lost {r.lost_hla} (Copy ratio {r.copy_ratio})\n")
        else:
            f.write("No LOH events detected.\n")

    print(f"LOH analysis complete: {loh_count}/{total_genes} genes with LOH")


if __name__ == '__main__':
    main()
