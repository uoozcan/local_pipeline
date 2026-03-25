#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parse seq2HLA 4-digit genotype files to standard pipeline TSV format.

seq2HLA produces three output files (prefix = value passed to -r):
    {prefix}-ClassI-class.HLAgenotype4digits    — classical Class I: A, B, C
    {prefix}-ClassI-nonclass.HLAgenotype4digits — non-classical: E, F, G, H, J, K, L, P, V (skipped)
    {prefix}-ClassII.HLAgenotype4digits          — Class II: DRB1, DQA1, DQB1, DPA1, DPB1, DRA

Each file is tab-separated with a header line:
    #Locus  Allele 1  Confidence  Allele 2  Confidence.1
    A       A*36:01   0.520369    A*11:02'  0.0
    B       B*08:01   0.0007      B*56:05'  0.0
    C       C*07:01'  0.0         C*07:01   NA

Alleles with a trailing apostrophe (') are ambiguity-flagged; stripped here.
Alleles where BOTH confidences are NA or <= CONF_THRESHOLD (0.1) are excluded.
Non-classical Class I loci (E/F/G/H/J/K/L/P/V) and non-classical Class II (DRA)
are excluded — only classical genes are reported.

Standard pipeline output format:
    # seq2HLA results for {sample}
    Gene    Allele1    Allele2    Reads1    Reads2
    A       A*36:01    A*11:02    NA        NA
"""

from __future__ import print_function
import argparse
import os
import sys


CONF_THRESHOLD = 0.1

# Only report classical HLA genes (exclude DRA and all non-classical Class I loci)
CLASSICAL_GENES = {'A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'}


def parse_args():
    parser = argparse.ArgumentParser(description='Parse seq2HLA output to standard format')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--prefix', required=True,
                        help='Output prefix used with seq2HLA -r flag (e.g. "sample_id.")')
    parser.add_argument('--output', required=True, help='Output file in standard pipeline format')
    return parser.parse_args()


def parse_class_file(filepath):
    """
    Parse one seq2HLA genotype4digits file.
    Returns list of (gene, allele1, allele2) tuples for rows passing confidence threshold.
    """
    entries = []
    if not os.path.exists(filepath):
        return entries

    try:
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                # Skip comment/header lines
                if not line or line.startswith('#'):
                    continue

                cols = line.split('\t')
                # Expected columns: Locus, Allele1, Conf1, Allele2, Conf2
                if len(cols) < 5:
                    continue

                gene  = cols[0].strip()

                # Skip non-classical genes
                if gene not in CLASSICAL_GENES:
                    continue
                a1    = cols[1].strip().rstrip("'")   # strip ambiguity apostrophe
                c1_s  = cols[2].strip()
                a2    = cols[3].strip().rstrip("'")
                c2_s  = cols[4].strip()

                # Parse confidence values.
                # 'NA' means seq2HLA could not compute a confidence (common on WGS
                # data due to low HLA-region mapping) — treat NA as unscored, NOT as 0.
                # Only filter when confidence is an explicit low numeric value.
                def parse_conf(s):
                    if s in ('NA', ''):
                        return None          # unscored — do not filter
                    try:
                        return float(s)
                    except ValueError:
                        return None

                c1 = parse_conf(c1_s)
                c2 = parse_conf(c2_s)

                # Skip only when BOTH confidences are numeric AND below threshold.
                # If either is NA (unscored) we keep the allele.
                c1_low = (c1 is not None and c1 <= CONF_THRESHOLD)
                c2_low = (c2 is not None and c2 <= CONF_THRESHOLD)
                if c1_low and c2_low:
                    print("  Skipping {}: both allele confidences below {} (c1={:.3f}, c2={:.3f})".format(
                        gene, CONF_THRESHOLD, c1, c2), file=sys.stderr)
                    continue

                # Skip loci where allele1 is 'no' (seq2HLA could not type it at all)
                if a1.lower() in ('no', 'not typed', '-', ''):
                    continue

                # Fall back to allele1 if allele2 is absent or 'not typed'
                if not a2 or a2 in ('-', 'NA', 'not typed', 'Not typed'):
                    a2 = a1

                entries.append((gene, a1, a2))

    except Exception as e:
        print("Warning: error reading {}: {}".format(filepath, e), file=sys.stderr)

    return entries


def main():
    args = parse_args()

    prefix = args.prefix
    # seq2HLA actual output filenames (observed from seq2HLA v2.3 run):
    #   classical Class I  → {prefix}-ClassI-class.HLAgenotype4digits
    #   non-classical ClassI → {prefix}-ClassI-nonclass.HLAgenotype4digits  (skipped)
    #   Class II           → {prefix}-ClassII.HLAgenotype4digits
    class1_file = "{}-ClassI-class.HLAgenotype4digits".format(prefix)
    class2_file = "{}-ClassII.HLAgenotype4digits".format(prefix)

    entries = []
    entries.extend(parse_class_file(class1_file))
    entries.extend(parse_class_file(class2_file))

    with open(args.output, 'w') as out:
        out.write("# seq2HLA results for {}\n".format(args.sample))
        out.write("# Tool: seq2HLA v2.3 - Class I + II, RNA-seq\n")
        out.write("# Confidence threshold: > {}\n".format(CONF_THRESHOLD))
        out.write("#\n")
        out.write("Gene\tAllele1\tAllele2\tReads1\tReads2\n")

        if not entries:
            print("Warning: no alleles above confidence threshold for {}".format(args.sample),
                  file=sys.stderr)
        else:
            for gene, a1, a2 in entries:
                out.write("{}\t{}\t{}\tNA\tNA\n".format(gene, a1, a2))

    total = len(entries)
    print("seq2HLA parsing complete: {} genes typed for {}".format(total, args.sample))


if __name__ == '__main__':
    main()
