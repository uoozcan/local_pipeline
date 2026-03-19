#!/usr/bin/env python3
"""
Generate HML (HLA Markup Language) v1.0.1 XML from consensus results.

HML is the interchange format used by NMDP/Be The Match transplant registries
and other clinical HLA laboratories for electronic data exchange.

Specification: http://schemas.nmdp.org/spec/hml/1.0.1

Usage:
    generate_hml.py --sample SAMPLE_ID \\
                    --consensus sample_consensus.txt \\
                    --output sample.hml.xml \\
                    [--allele-db-version 3.57.0] \\
                    [--center-id MY_LAB] \\
                    [--typing-method NextGen]
"""

import argparse
import sys
from datetime import date
from pathlib import Path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate HML v1.0.1 XML from HLA consensus results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--sample', required=True,
                        help='Sample identifier (written into <sample id="...">)')
    parser.add_argument('--consensus', required=True,
                        help='Consensus TSV produced by consensus_voting.py')
    parser.add_argument('--output', required=True,
                        help='Output HML XML file path')
    parser.add_argument('--allele-db', default='IMGT/HLA',
                        help='Allele database name')
    parser.add_argument('--allele-db-version', default='3.57.0',
                        help='IMGT/HLA release version (e.g. 3.57.0)')
    parser.add_argument('--center-id', default='HLA-PIPELINE',
                        help='Reporting centre identifier')
    parser.add_argument('--typing-method', default='NextGen',
                        help='Typing method description')
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Consensus file parsing
# ---------------------------------------------------------------------------

def parse_consensus_file(filepath: str):
    """
    Parse the TSV file produced by consensus_voting.py.

    Expected format (after comment header lines):
        Gene  Allele1  Allele2  Confidence  Reads1  Reads2

    Returns list of (gene, allele1_or_None, allele2_or_None, confidence).
    """
    results = []
    try:
        with open(filepath, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if line.startswith('Gene\t'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0]
                a1 = parts[1] if parts[1] not in ('-', '', 'None') else None
                a2 = parts[2] if parts[2] not in ('-', '', 'None') else None
                try:
                    confidence = float(parts[3]) if len(parts) > 3 else 0.0
                except ValueError:
                    confidence = 0.0
                results.append((gene, a1, a2, confidence))
    except OSError as exc:
        print(f"Error reading consensus file: {exc}", file=sys.stderr)
        sys.exit(1)
    return results


# ---------------------------------------------------------------------------
# GL String helpers
# ---------------------------------------------------------------------------

def _hla_prefix(allele: str) -> str:
    """Ensure allele has the 'HLA-' prefix required by the GL String spec."""
    return allele if allele.startswith('HLA-') else f"HLA-{allele}"


def locus_gl_string(a1, a2) -> str:
    """
    Build the GL String fragment for a single locus.

    Heterozygous: HLA-A*03:01+HLA-A*11:01
    Homozygous / single call: HLA-A*03:01
    """
    if not a1:
        return ''
    p1 = _hla_prefix(a1)
    if a2 and a2 != '-' and a2 != a1:
        return f"{p1}+{_hla_prefix(a2)}"
    return p1


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def _esc(text: str) -> str:
    """Escape the five predefined XML entities."""
    return (text
            .replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
            .replace('"', '&quot;')
            .replace("'", '&apos;'))


# ---------------------------------------------------------------------------
# HML generation
# ---------------------------------------------------------------------------

def build_hml_xml(sample_id: str,
                  consensus_results: list,
                  allele_db: str,
                  allele_db_version: str,
                  center_id: str,
                  typing_method: str) -> str:
    """
    Produce an HML v1.0.1 compliant XML string.

    Schema: http://schemas.nmdp.org/spec/hml/1.0.1

    The confidence attribute on <allele-assignment> is a pipeline extension
    (not part of the base HML schema) and records the per-locus confidence
    score (0.0–1.0) computed by the weighted consensus algorithm.
    """
    today = date.today().isoformat()
    sid   = _esc(sample_id)
    cid   = _esc(center_id)
    adb   = _esc(allele_db)
    adbv  = _esc(allele_db_version)
    meth  = _esc(typing_method)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<hml',
        '  xmlns="http://schemas.nmdp.org/spec/hml/1.0.1"',
        '  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',
        '  xsi:schemaLocation="http://schemas.nmdp.org/spec/hml/1.0.1'
        ' https://raw.githubusercontent.com/nmdp-bioinformatics/hml-fhir-app'
        '/master/src/test/resources/schemas/hml-1.0.1.xsd"',
        '  project-name="HLA-Typing-Pipeline"',
        '  version="1.0.1">',
        f'  <reporting-center',
        f'    reporting-center-id="{cid}"',
        f'    reporting-center-context="NMDP"/>',
        f'  <sample id="{sid}" center-code="{cid}">',
    ]

    typed_count = 0
    for gene, a1, a2, confidence in consensus_results:
        if not a1:
            continue
        gl_str = locus_gl_string(a1, a2)
        if not gl_str:
            continue
        typed_count += 1
        locus = gene if gene.startswith('HLA-') else f"HLA-{gene}"
        lines += [
            f'    <typing locus="{_esc(locus)}" typing-method="{meth}">',
            f'      <allele-assignment',
            f'        date="{today}"',
            f'        allele-db="{adb}"',
            f'        allele-version="{adbv}"',
            f'        confidence="{confidence:.3f}">',
            f'        <glstring>{_esc(gl_str)}</glstring>',
            f'      </allele-assignment>',
            f'    </typing>',
        ]

    lines += [
        '  </sample>',
        '</hml>',
    ]

    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not Path(args.consensus).exists():
        print(f"Error: consensus file not found: {args.consensus}", file=sys.stderr)
        sys.exit(1)

    consensus_results = parse_consensus_file(args.consensus)
    if not consensus_results:
        print(f"Warning: no consensus results found in {args.consensus}",
              file=sys.stderr)

    xml_content = build_hml_xml(
        sample_id=args.sample,
        consensus_results=consensus_results,
        allele_db=args.allele_db,
        allele_db_version=args.allele_db_version,
        center_id=args.center_id,
        typing_method=args.typing_method,
    )

    with open(args.output, 'w', encoding='utf-8') as fh:
        fh.write(xml_content)
        fh.write('\n')

    typed = sum(1 for _, a1, _, _ in consensus_results if a1)
    print(f"HML v1.0.1 written to {args.output} ({typed} loci typed)")


if __name__ == '__main__':
    main()
