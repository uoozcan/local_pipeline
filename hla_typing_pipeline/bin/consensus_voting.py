#!/usr/bin/env python3
"""
HLA Consensus Voting Script with Read-Weighted Confidence
Combines results from multiple HLA typing tools using weighted voting
based on read support confidence scores.

The confidence score is calculated as:
    confidence = actual_reads / expected_reads (capped at 1.0)

Where expected_reads is determined by the --expected-reads parameter or
calculated from the total HLA reads.

Usage:
    consensus_voting.py --sample SAMPLE --tools tool1,tool2 --files file1,file2 \
                        --resolution 2-field --min-tools 1 \
                        --expected-reads 1000 \
                        --output-consensus out.txt --output-comparison comp.txt
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any


def parse_args():
    parser = argparse.ArgumentParser(description='HLA Consensus Voting with Read-Weighted Confidence')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--tools', required=True, help='Comma-separated list of tools')
    parser.add_argument('--files', required=True, help='Comma-separated list of result files')
    parser.add_argument('--resolution', default='2-field',
                        choices=['2-field', '4-field'], help='Allele resolution')
    parser.add_argument('--min-tools', type=int, default=1,
                        help='Minimum tools required for consensus')
    parser.add_argument('--expected-reads', type=int, default=1000,
                        help='Expected reads per allele for confidence calculation')
    parser.add_argument('--weighting', default='read_confidence',
                        choices=['equal', 'read_confidence', 'tool_quality', 'calibrated'],
                        help='Weighting method for consensus voting. '
                             '"calibrated" requires --weights-file.')
    parser.add_argument('--weights-file',
                        help='JSON weights file produced by calibrate_tool_weights.py '
                             '(required when --weighting calibrated)')
    parser.add_argument('--data-type', default='wgs', choices=['wgs', 'rna'],
                        help='Sequencing data type; used for logging when weighting=calibrated')
    parser.add_argument('--output-consensus', required=True, help='Output consensus file')
    parser.add_argument('--output-comparison', required=True, help='Output comparison file')
    parser.add_argument('--output-format', default='text',
                        choices=['text', 'gl_string', 'hml', 'all'],
                        help='Output format(s). "gl_string" adds a GL String file; '
                             '"hml" adds GL String + HML v1.0.1 XML; "all" produces all formats.')
    parser.add_argument('--output-gl-string', default=None,
                        help='Path for GL String output file (auto-named when not specified)')
    return parser.parse_args()


def normalize_allele(allele: str, resolution: str = '2-field') -> Optional[str]:
    """Normalize HLA allele to specified resolution."""
    if not allele or allele in ['-', 'NA', 'None', '', 'Not typed', '?']:
        return None

    # Remove HLA- prefix if present
    allele = re.sub(r'^HLA-', '', allele)

    # Extract gene and allele fields
    match = re.match(r'^([A-Z]+\d*)\*(\d+):(\d+)(?::(\d+))?(?::(\d+))?([NLSCAQ])?$', allele)
    if not match:
        # Try alternative format
        match = re.match(r'^([A-Z]+\d*)\*(\d+):(\d+)', allele)
        if not match:
            return allele  # Return as-is if can't parse

    gene = match.group(1)
    field1 = match.group(2)
    field2 = match.group(3)

    if resolution == '2-field':
        return f"{gene}*{field1}:{field2}"
    else:  # 4-field
        field3 = match.group(4) if len(match.groups()) > 3 and match.group(4) else None
        field4 = match.group(5) if len(match.groups()) > 4 and match.group(5) else None
        if field3 and field4:
            return f"{gene}*{field1}:{field2}:{field3}:{field4}"
        elif field3:
            return f"{gene}*{field1}:{field2}:{field3}"
        else:
            return f"{gene}*{field1}:{field2}"


class AlleleCall:
    """Represents an HLA allele call with associated metrics."""
    def __init__(self, allele: str, reads: int = 0, coverage: float = 0.0,
                 quality: float = 1.0, tool: str = ''):
        self.allele = allele
        self.reads = reads
        self.coverage = coverage
        self.quality = quality  # Tool-reported quality/confidence
        self.tool = tool

    def __repr__(self):
        return f"AlleleCall({self.allele}, reads={self.reads}, quality={self.quality:.2f})"


def parse_hlahd_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse HLA-HD result file with read counts.
    HLA-HD format: Gene\tAllele1\tAllele2\t[additional columns may include reads]
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    seen_alleles = set()
                    for i, allele in enumerate(parts[1:3]):
                        allele = allele.strip()
                        if allele and allele not in ['-', 'Not typed', 'NA']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                # Try to extract read count from additional columns
                                reads = 0
                                if len(parts) > 3 + i:
                                    try:
                                        reads = int(parts[3 + i])
                                    except (ValueError, IndexError):
                                        pass

                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=reads,
                                    quality=1.0 if reads > 0 else 0.5,
                                    tool='hlahd'
                                ))
                                seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing HLA-HD file {filepath}: {e}", file=sys.stderr)
    return results


def parse_arcashla_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse arcasHLA result file.
    Also tries to parse the JSON file if available for read counts.
    """
    results = defaultdict(list)

    # Try to find and parse JSON file (published alongside txt as {sample}_arcashla.json)
    json_data = {}
    parent_dir = Path(filepath).parent
    sample_name = Path(filepath).stem.replace('_arcashla', '')
    possible_jsons = [
        # Primary: same dir as txt, _arcashla.json suffix (written by module)
        parent_dir / f"{sample_name}_arcashla.json",
        # Fallback: arcasHLA native genotype.json locations
        parent_dir / f"{sample_name}.genotype.json",
        parent_dir / sample_name / f"{sample_name}.genotype.json",
    ]

    for jp in possible_jsons:
        if jp.exists():
            try:
                with open(jp, 'r') as f:
                    json_data = json.load(f)
                break
            except Exception:
                pass

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    gene_short = gene.replace('HLA-', '')

                    for allele in parts[1:]:
                        allele = allele.strip()
                        if allele and allele != '-':
                            normalized = normalize_allele(allele, resolution)
                            if normalized:
                                # arcasHLA genotype.json values are allele lists,
                                # not dicts with read counts — reads unavailable
                                reads = 0

                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=reads,
                                    quality=1.0,
                                    tool='arcashla'
                                ))
    except Exception as e:
        print(f"Warning: Error parsing arcasHLA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_optitype_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse OptiType result file.
    OptiType provides an objective score which can be used as confidence.
    """
    results = defaultdict(list)
    objective_score = 0.0

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    # Check for objective score in comments
                    if 'Objective' in line:
                        try:
                            objective_score = float(line.split(':')[-1].strip())
                        except:
                            pass
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    # Try to extract reads/objective from additional columns
                    reads = 0
                    if len(parts) > 3:
                        try:
                            reads = int(float(parts[3]))
                        except (ValueError, IndexError):
                            pass

                    for allele in parts[1:3]:
                        allele = allele.strip()
                        if allele and allele != '-':
                            # OptiType format handling
                            if '*' not in allele and ':' in allele:
                                gene_short = gene.replace('HLA-', '')
                                allele = f"{gene_short}*{allele}"
                            normalized = normalize_allele(allele, resolution)
                            if normalized:
                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=reads,
                                    quality=min(1.0, objective_score / 100) if objective_score else 1.0,
                                    tool='optitype'
                                ))
    except Exception as e:
        print(f"Warning: Error parsing OptiType file {filepath}: {e}", file=sys.stderr)
    return results


def parse_spechla_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse SpecHLA result file.

    Supports two formats:
    1. Wide format (hla.result.txt): single sample row with columns HLA_A_1, HLA_A_2, HLA_B_1, ...
       Header: Sample  HLA_A_1  HLA_A_2  HLA_B_1  HLA_B_2  ...
       Data:   NA12878 A*01:01  A*03:01  B*07:02  ...
    2. Long format (one row per gene):
       Gene    Allele1  Allele2  [Reads1  Reads2  Quality]
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]

        if not lines:
            return results

        # Detect wide format: header starts with 'Sample' and contains 'HLA_A_1'
        first_parts = lines[0].split('\t')
        is_wide = (first_parts[0].strip().lower() == 'sample' and
                   any('HLA_' in p for p in first_parts))

        if is_wide:
            # Wide format: parse header → map col index to (gene, allele_num)
            header = first_parts
            col_map: Dict[int, tuple] = {}  # col_idx -> ('HLA-A', 1) or ('HLA-A', 2)
            for idx, col in enumerate(header):
                col = col.strip()
                # e.g. HLA_A_1, HLA_DRB1_2
                if col.startswith('HLA_'):
                    parts2 = col.split('_')
                    # HLA_A_1 → ['HLA', 'A', '1']
                    # HLA_DRB1_1 → ['HLA', 'DRB1', '1']
                    if len(parts2) >= 3:
                        gene_name = '_'.join(parts2[1:-1])  # e.g. 'A', 'DRB1'
                        allele_num = int(parts2[-1])
                        col_map[idx] = (f"HLA-{gene_name}", allele_num)

            # Parse data rows (skip header)
            for line in lines[1:]:
                parts = line.split('\t')
                seen: Dict[str, set] = defaultdict(set)
                for idx, (gene, _allele_num) in col_map.items():
                    if idx >= len(parts):
                        continue
                    allele = parts[idx].strip()
                    if allele and allele not in ['-', 'NA', '']:
                        normalized = normalize_allele(allele, resolution)
                        if normalized and normalized not in seen[gene]:
                            results[gene].append(AlleleCall(
                                allele=normalized,
                                reads=0,
                                quality=1.0,
                                tool='spechla'
                            ))
                            seen[gene].add(normalized)
        else:
            # Long format: Gene  Allele1  Allele2  [Reads1  Reads2  Quality]
            header: Optional[List[str]] = None
            for line in lines:
                parts = line.split('\t')
                if parts[0].lower() in ['gene', 'locus']:
                    header = [p.lower() for p in parts]
                    continue

                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"
                    seen_alleles: set = set()
                    for i, allele in enumerate(parts[1:3]):
                        allele = allele.strip()
                        if allele and allele not in ['-', 'NA']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                reads = 0
                                if len(parts) > 3 + i:
                                    try:
                                        reads = int(parts[3 + i])
                                    except (ValueError, IndexError):
                                        pass
                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=reads,
                                    quality=1.0,
                                    tool='spechla'
                                ))
                                seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing SpecHLA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_hlala_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse HLA*LA result file.
    HLA*LA provides quality scores for each allele call.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            header = None
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                parts = line.split('\t')

                # Check for header
                if parts[0].lower() in ['gene', 'locus']:
                    header = [p.lower() for p in parts]
                    continue

                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    seen_alleles = set()
                    for i, allele in enumerate(parts[1:3]):
                        allele = allele.strip()
                        if allele and allele not in ['-', '?', 'NA']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                quality = 1.0
                                reads = 0

                                # Try to get quality from additional columns
                                if header and 'quality' in header:
                                    try:
                                        quality = float(parts[header.index('quality')])
                                    except (ValueError, IndexError):
                                        pass
                                elif len(parts) > 3:
                                    try:
                                        quality = float(parts[3])
                                    except (ValueError, IndexError):
                                        pass

                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=reads,
                                    quality=quality,
                                    tool='hlala'
                                ))
                                seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing HLA*LA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_xhla_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse xHLA result file.
    xHLA does not report read counts (always NA), so reads is set to 0
    and quality is set to 0.5 to reflect lower confidence in consensus voting.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                parts = line.split('\t')

                if parts[0].lower() in ['gene', 'locus']:
                    continue

                if len(parts) >= 3:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    seen_alleles = set()
                    for allele in parts[1:3]:
                        allele = allele.strip()
                        if allele and allele not in ['-', 'NA']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=0,
                                    quality=0.5,
                                    tool='xhla'
                                ))
                                seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing xHLA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_hlascan_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse HLAscan result file (standard TSV format produced by parse_hlascan_results.py).
    HLAscan does not report read counts; quality set to 0.5.
    Note: HLAscan v2.1 only supports hg19; hg38 results will be all NA.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue

                parts = line.split('\t')

                if parts[0].lower() in ['gene', 'locus']:
                    continue

                if len(parts) >= 3:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    seen_alleles = set()
                    for allele in parts[1:3]:
                        allele = allele.strip()
                        if allele and allele not in ['-', 'NA']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                results[gene].append(AlleleCall(
                                    allele=normalized,
                                    reads=0,
                                    quality=0.5,
                                    tool='hlascan'
                                ))
                                seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing HLAscan file {filepath}: {e}", file=sys.stderr)
    return results


def parse_t1k_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse T1K genotype TSV file.
    T1K format: Gene\tAllele1\tAbundance1\tAllele2\tAbundance2
    Abundance is a float (read-equivalent score); used as quality proxy.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"
                seen_alleles = set()
                for i, allele_idx in enumerate([1, 3]):
                    if allele_idx >= len(parts):
                        break
                    allele = parts[allele_idx].strip()
                    if not allele or allele in ['-', 'NA', 'None', '.']:
                        continue
                    normalized = normalize_allele(allele, resolution)
                    if normalized and normalized not in seen_alleles:
                        abundance_idx = allele_idx + 1
                        try:
                            abundance = float(parts[abundance_idx]) if abundance_idx < len(parts) else 0.0
                        except (ValueError, IndexError):
                            abundance = 0.0
                        results[gene].append(AlleleCall(
                            allele=normalized,
                            reads=int(abundance),
                            quality=min(1.0, abundance / 50.0) if abundance > 0 else 0.5,
                            tool='t1k'
                        ))
                        seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing T1K file {filepath}: {e}", file=sys.stderr)
    return results


def parse_hifihla_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse HiFi-HLA result file produced by parse_hifihla_results.py.

    Format:
        Gene  Allele1  Allele2  Conf1  Conf2

    HiFi reads have intrinsically high base accuracy (Q20+), so quality is
    set to 1.0 and the confidence score from HiFi-HLA (col 3) is used
    as a scaled read count for read_confidence weighting.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                if line.startswith('Gene\t'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue

                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"

                seen_alleles: set = set()
                for i, (allele_col, conf_col) in enumerate([(1, 3), (2, 4)]):
                    if allele_col >= len(parts):
                        break
                    allele = parts[allele_col].strip()
                    if not allele or allele in ('-', 'NA', 'None', ''):
                        continue
                    normalized = normalize_allele(allele, resolution)
                    if not normalized or normalized in seen_alleles:
                        continue

                    # Use HiFi-HLA confidence (0–1) scaled to a pseudo read count
                    # so that read_confidence weighting treats these like well-supported alleles
                    try:
                        conf = float(parts[conf_col]) if conf_col < len(parts) else 0.9
                    except (ValueError, IndexError):
                        conf = 0.9
                    conf = min(1.0, max(0.0, conf))
                    pseudo_reads = int(conf * 100)  # e.g. 0.99 → 99

                    results[gene].append(AlleleCall(
                        allele=normalized,
                        reads=pseudo_reads,
                        quality=1.0,   # long-read base accuracy is uniformly high
                        tool='hifihla'
                    ))
                    seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing HiFi-HLA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_polysolver_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse POLYSOLVER standard TSV output (produced by parse_polysolver_results.py).

    Format (after conversion):
        Gene    Allele1    Allele2    Reads1    Reads2
        A       A*03:01    A*03:01    NA        NA
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"
                seen_alleles: set = set()
                for allele in parts[1:3]:
                    allele = allele.strip()
                    if allele and allele not in ('-', 'NA', 'None'):
                        normalized = normalize_allele(allele, resolution)
                        if normalized and normalized not in seen_alleles:
                            results[gene].append(AlleleCall(
                                allele=normalized,
                                reads=0,        # POLYSOLVER does not report read counts
                                quality=1.0,
                                tool='polysolver'
                            ))
                            seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing POLYSOLVER file {filepath}: {e}", file=sys.stderr)
    return results


def parse_seq2hla_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse seq2HLA standard TSV output (produced by parse_seq2hla_results.py).

    Format (after conversion):
        Gene    Allele1      Allele2      Reads1    Reads2
        A       A*36:01      A*11:02      NA        NA
        DRB1    DRB1*09:01   DRB1*01:01   NA        NA
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"
                seen_alleles: set = set()
                for allele in parts[1:3]:
                    allele = allele.strip()
                    if allele and allele not in ('-', 'NA', 'None'):
                        normalized = normalize_allele(allele, resolution)
                        if normalized and normalized not in seen_alleles:
                            results[gene].append(AlleleCall(
                                allele=normalized,
                                reads=0,        # seq2HLA does not report read counts
                                quality=1.0,
                                tool='seq2hla'
                            ))
                            seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing seq2HLA file {filepath}: {e}", file=sys.stderr)
    return results


def parse_kourami_results(filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """
    Parse Kourami standard TSV output (produced by parse_kourami_results.py).

    Format (after conversion):
        Gene    Allele1      Allele2      Reads1    Reads2
        A       A*11:01:01   A*03:01:01   546       532

    Read counts (matched_bases) are available for read_confidence weighting.
    """
    results = defaultdict(list)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or line.startswith('Gene') or not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                gene = parts[0].strip()
                if not gene.startswith('HLA-'):
                    gene = f"HLA-{gene}"
                seen_alleles: set = set()
                for i, allele in enumerate(parts[1:3]):
                    allele = allele.strip()
                    if allele and allele not in ('-', 'NA', 'None'):
                        normalized = normalize_allele(allele, resolution)
                        if normalized and normalized not in seen_alleles:
                            reads = 0
                            if len(parts) > 3 + i:
                                try:
                                    reads = int(parts[3 + i])
                                except (ValueError, IndexError):
                                    pass
                            results[gene].append(AlleleCall(
                                allele=normalized,
                                reads=reads,
                                quality=1.0 if reads > 0 else 0.5,
                                tool='kourami'
                            ))
                            seen_alleles.add(normalized)
    except Exception as e:
        print(f"Warning: Error parsing Kourami file {filepath}: {e}", file=sys.stderr)
    return results


def parse_results(tool: str, filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """Parse results based on tool type."""
    parsers = {
        'spechla':    parse_spechla_results,
        'hlahd':      parse_hlahd_results,
        'hlala':      parse_hlala_results,
        'arcashla':   parse_arcashla_results,
        'optitype':   parse_optitype_results,
        'xhla':       parse_xhla_results,
        'hlascan':    parse_hlascan_results,
        't1k':        parse_t1k_results,
        'hifihla':    parse_hifihla_results,
        'polysolver': parse_polysolver_results,
        'seq2hla':    parse_seq2hla_results,
        'kourami':    parse_kourami_results,
    }
    parser = parsers.get(tool.lower())
    if parser:
        return parser(filepath, resolution)
    else:
        print(f"Warning: Unknown tool {tool}, attempting generic parsing", file=sys.stderr)
        return parse_arcashla_results(filepath, resolution)


def calculate_read_confidence(reads: int, expected_reads: int) -> float:
    """
    Calculate confidence score based on read support.

    confidence = min(1.0, actual_reads / expected_reads)

    This gives higher confidence to alleles with more read support.
    """
    if expected_reads <= 0:
        return 1.0 if reads > 0 else 0.0
    return min(1.0, reads / expected_reads)


def generate_gl_string(consensus_results: Dict, ordered_genes: List[str]) -> str:
    """
    Generate a GL String (GL String Consortium format) from consensus results.

    Format:  HLA-A*03:01+HLA-A*11:01^HLA-B*07:02+HLA-B*35:01^...
    Loci are separated by '^'; two alleles at a locus are separated by '+'.
    The 'HLA-' prefix is required on every allele.
    Loci with no call are omitted.

    See: https://glstring.org/
    """
    def with_hla_prefix(allele: str) -> str:
        return f"HLA-{allele}" if not allele.startswith('HLA-') else allele

    locus_parts = []
    for gene in ordered_genes:
        result = consensus_results.get(gene)
        if not result:
            continue
        a1 = result.get('allele1')
        a2 = result.get('allele2')
        if not a1 or a1 == '-':
            continue
        if a2 and a2 != '-' and a2 != a1:
            locus_parts.append(f"{with_hla_prefix(a1)}+{with_hla_prefix(a2)}")
        else:
            locus_parts.append(with_hla_prefix(a1))

    return '^'.join(locus_parts)


def load_calibrated_weights(weights_file: str) -> Dict[str, Dict[str, float]]:
    """
    Load per-gene, per-tool weights from a JSON file produced by
    calibrate_tool_weights.py.

    Returns:
        {gene_short: {tool: weight}}
        e.g. {'A': {'hlahd': 0.34, 'spechla': 0.33, ...}, 'B': {...}, ...}
    """
    with open(weights_file, 'r') as f:
        payload = json.load(f)
    genes_section = payload.get('genes', {})
    # Keys in the JSON are short gene names (A, B, C, DRB1, …).
    # consensus_voting uses HLA-prefixed names internally, so we build
    # a lookup for both forms.
    weights: Dict[str, Dict[str, float]] = {}
    for gene, tool_weights in genes_section.items():
        weights[gene] = {t: float(w) for t, w in tool_weights.items()}
        weights[f"HLA-{gene}"] = weights[gene]  # alias
    return weights


def weighted_consensus_vote(
    allele_data: Dict[str, List[AlleleCall]],
    min_tools: int = 1,
    expected_reads: int = 1000,
    weighting: str = 'read_confidence',
    calibrated_weights: Optional[Dict[str, Dict[str, float]]] = None,
    gene: str = '',
) -> Tuple[Optional[str], Optional[str], float, Dict[str, Any]]:
    """
    Perform weighted consensus voting on allele calls from different tools.

    Weighting methods:
    - 'equal':        Each tool gets equal weight (traditional majority voting)
    - 'read_confidence': Weight by read support confidence
    - 'tool_quality': Weight by tool-reported quality scores
    - 'calibrated':   Use empirically calibrated per-gene weights from a JSON
                      file (produced by calibrate_tool_weights.py).  Requires
                      calibrated_weights and gene to be supplied.

    Returns:
        (allele1, allele2, confidence, details)
    """
    if not allele_data:
        return None, None, 0.0, {}

    # Collect all alleles with their weights
    allele_weights = defaultdict(float)
    allele_reads = defaultdict(int)
    allele_tools = defaultdict(list)
    tool_alleles = {}

    for tool, calls in allele_data.items():
        tool_alleles[tool] = set()
        for call in calls:
            tool_alleles[tool].add(call.allele)

            # Calculate weight based on weighting method
            if weighting == 'equal':
                weight = 1.0
            elif weighting == 'read_confidence':
                read_conf = calculate_read_confidence(call.reads, expected_reads)
                # Combine read confidence with tool quality
                weight = (read_conf * 0.7 + call.quality * 0.3) if call.reads > 0 else call.quality * 0.5
            elif weighting == 'tool_quality':
                weight = call.quality
            elif weighting == 'calibrated':
                # Base weight from empirical calibration (per-gene, per-tool)
                gene_key = gene or call.tool
                gene_weights = (calibrated_weights or {}).get(gene_key, {})
                base = gene_weights.get(call.tool, 0.0)
                # Scale by read evidence when available to reward higher-confidence calls
                if call.reads > 0:
                    read_conf = calculate_read_confidence(call.reads, expected_reads)
                    weight = base * (0.7 * read_conf + 0.3)
                else:
                    # No read support — halve the base weight (xHLA, HLA*LA)
                    weight = base * 0.5
            else:
                weight = 1.0

            allele_weights[call.allele] += weight
            allele_reads[call.allele] += call.reads
            allele_tools[call.allele].append(tool)

    if not allele_weights:
        return None, None, 0.0, {'tool_alleles': tool_alleles}

    # Sort alleles by weighted score
    sorted_alleles = sorted(allele_weights.items(), key=lambda x: x[1], reverse=True)

    num_tools = len(allele_data)
    max_weight = num_tools  # Maximum possible weight

    consensus1 = None
    consensus2 = None
    confidence = 0.0

    # Select top alleles that meet minimum tool threshold
    if len(sorted_alleles) >= 1:
        allele1, weight1 = sorted_alleles[0]
        tools_count1 = len(allele_tools[allele1])
        if tools_count1 >= min_tools:
            consensus1 = allele1
            confidence = weight1 / max_weight if max_weight > 0 else 0.0

    if len(sorted_alleles) >= 2:
        allele2, weight2 = sorted_alleles[1]
        tools_count2 = len(allele_tools[allele2])
        if tools_count2 >= min_tools:
            consensus2 = allele2
            # Combined confidence
            total_weight = weight1 + weight2 if consensus1 else weight2
            confidence = total_weight / (2 * max_weight) if max_weight > 0 else 0.0

    details = {
        'tool_alleles': tool_alleles,
        'allele_weights': dict(allele_weights),
        'allele_reads': dict(allele_reads),
        'allele_tools': dict(allele_tools),
        'weighting_method': weighting
    }

    return consensus1, consensus2, min(1.0, confidence), details


def main():
    args = parse_args()

    tools = args.tools.split(',')
    files = args.files.split(',')

    if len(tools) != len(files):
        print("Error: Number of tools must match number of files", file=sys.stderr)
        sys.exit(1)

    # Load calibrated weights if requested
    calibrated_weights = None
    if args.weighting == 'calibrated':
        if not args.weights_file:
            print("Error: --weights-file is required when --weighting calibrated", file=sys.stderr)
            sys.exit(1)
        if not Path(args.weights_file).exists():
            print(f"Error: weights file not found: {args.weights_file}", file=sys.stderr)
            sys.exit(1)
        calibrated_weights = load_calibrated_weights(args.weights_file)
        print(f"Loaded calibrated weights from {args.weights_file}", file=sys.stderr)

    # Parse all results
    all_results = {}
    for tool, filepath in zip(tools, files):
        results = parse_results(tool, filepath, args.resolution)
        all_results[tool] = results

    # Get all genes across all tools
    all_genes = set()
    for results in all_results.values():
        all_genes.update(results.keys())

    # Classical HLA genes in order
    classical_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1',
                       'HLA-DPA1', 'HLA-DPB1']
    ordered_genes = [g for g in classical_genes if g in all_genes]
    ordered_genes.extend([g for g in sorted(all_genes) if g not in ordered_genes])

    # Perform weighted consensus voting
    consensus_results = {}
    for gene in ordered_genes:
        gene_alleles = {tool: results.get(gene, [])
                        for tool, results in all_results.items()
                        if gene in results}

        consensus1, consensus2, confidence, details = weighted_consensus_vote(
            gene_alleles,
            args.min_tools,
            args.expected_reads,
            args.weighting,
            calibrated_weights=calibrated_weights,
            gene=gene,
        )

        consensus_results[gene] = {
            'allele1': consensus1,
            'allele2': consensus2,
            'confidence': confidence,
            'details': details
        }

    # Write consensus output
    with open(args.output_consensus, 'w') as f:
        f.write(f"# HLA Consensus Results for {args.sample}\n")
        f.write(f"# Tools used: {', '.join(tools)}\n")
        f.write(f"# Resolution: {args.resolution}\n")
        f.write(f"# Weighting method: {args.weighting}\n")
        if args.weighting == 'calibrated' and args.weights_file:
            f.write(f"# Weights file: {args.weights_file}\n")
            f.write(f"# Data type: {args.data_type}\n")
        f.write(f"# Expected reads per allele: {args.expected_reads}\n")
        f.write(f"# Minimum tools for consensus: {args.min_tools}\n")
        f.write("#\n")
        f.write("Gene\tAllele1\tAllele2\tConfidence\tReads1\tReads2\n")

        for gene in ordered_genes:
            result = consensus_results[gene]
            a1 = result['allele1'] or '-'
            a2 = result['allele2'] or '-'
            conf = f"{result['confidence']:.2f}"

            # Get read counts for consensus alleles
            reads1 = result['details'].get('allele_reads', {}).get(result['allele1'], 0) if result['allele1'] else 0
            reads2 = result['details'].get('allele_reads', {}).get(result['allele2'], 0) if result['allele2'] else 0

            f.write(f"{gene}\t{a1}\t{a2}\t{conf}\t{reads1}\t{reads2}\n")

    # Write comparison output
    with open(args.output_comparison, 'w') as f:
        f.write(f"# HLA Tool Comparison for {args.sample}\n")
        f.write(f"# Resolution: {args.resolution}\n")
        f.write(f"# Weighting method: {args.weighting}\n")
        f.write("#\n")

        # Header
        header = ['Gene'] + [f"{t}(alleles)" for t in tools] + [f"{t}(reads)" for t in tools] + ['Consensus', 'Weight', 'Agreement']
        f.write('\t'.join(header) + '\n')

        for gene in ordered_genes:
            row = [gene]
            result = consensus_results[gene]
            details = result['details']
            tool_alleles = details.get('tool_alleles', {})

            # Tool allele results
            for tool in tools:
                tool_calls = all_results[tool].get(gene, [])
                if tool_calls:
                    row.append('/'.join(sorted(set(c.allele for c in tool_calls))))
                else:
                    row.append('-')

            # Tool read counts
            for tool in tools:
                tool_calls = all_results[tool].get(gene, [])
                if tool_calls:
                    total_reads = sum(c.reads for c in tool_calls)
                    row.append(str(total_reads) if total_reads > 0 else '-')
                else:
                    row.append('-')

            # Consensus
            consensus = []
            if result['allele1']:
                consensus.append(result['allele1'])
            if result['allele2']:
                consensus.append(result['allele2'])
            row.append('/'.join(consensus) if consensus else '-')

            # Weight score
            if result['allele1'] and details.get('allele_weights'):
                weight = details['allele_weights'].get(result['allele1'], 0)
                row.append(f"{weight:.2f}")
            else:
                row.append('-')

            # Agreement
            agreement = f"{result['confidence']*100:.0f}%"
            row.append(agreement)

            f.write('\t'.join(row) + '\n')

    # Write GL String output when requested
    output_format = args.output_format
    if output_format in ('gl_string', 'hml', 'all'):
        gl_string = generate_gl_string(consensus_results, ordered_genes)
        gl_file = (args.output_gl_string or
                   args.output_consensus.replace('_consensus.txt', '_gl_string.txt'))
        with open(gl_file, 'w') as f:
            f.write(f"# GL String for {args.sample}\n")
            f.write(f"# Tools: {', '.join(tools)}\n")
            f.write(f"# Resolution: {args.resolution}\n")
            f.write(f"# Weighting: {args.weighting}\n")
            f.write(f"# Format: GL String Consortium (GLSC) — https://glstring.org/\n")
            f.write(f"{args.sample}\t{gl_string}\n")
        print(f"GL String written to {gl_file}")

    print(f"Consensus results written to {args.output_consensus}")
    print(f"Comparison results written to {args.output_comparison}")
    print(f"Weighting method used: {args.weighting}")


if __name__ == '__main__':
    main()
