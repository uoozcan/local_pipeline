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
                        choices=['equal', 'read_confidence', 'tool_quality'],
                        help='Weighting method for consensus voting')
    parser.add_argument('--output-consensus', required=True, help='Output consensus file')
    parser.add_argument('--output-comparison', required=True, help='Output comparison file')
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

    # Try to find and parse JSON file for detailed read counts
    json_path = Path(filepath).with_suffix('.json')
    json_data = {}

    # Also check for genotype.json pattern
    parent_dir = Path(filepath).parent
    sample_name = Path(filepath).stem.replace('_arcashla', '')
    possible_jsons = [
        json_path,
        parent_dir / f"{sample_name}.genotype.json",
        parent_dir / sample_name / f"{sample_name}.genotype.json"
    ]

    for jp in possible_jsons:
        if jp.exists():
            try:
                with open(jp, 'r') as f:
                    json_data = json.load(f)
                break
            except:
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
                                # Try to get read count from JSON data
                                reads = 0
                                if gene_short in json_data:
                                    # arcasHLA JSON may have read counts
                                    gene_info = json_data.get(gene_short, {})
                                    if isinstance(gene_info, dict) and 'reads' in gene_info:
                                        reads = gene_info.get('reads', 0)

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
    Parse SpecHLA result file with read counts.
    SpecHLA output typically includes: Gene\tAllele1\tAllele2\tReads1\tReads2\tQuality
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

                # Check if this is a header line
                if parts[0].lower() in ['gene', 'locus']:
                    header = [p.lower() for p in parts]
                    continue

                if len(parts) >= 2:
                    gene = parts[0].strip()
                    if not gene.startswith('HLA-'):
                        gene = f"HLA-{gene}"

                    seen_alleles = set()

                    # Parse with or without header
                    if header:
                        # Use header to find columns
                        for i, allele in enumerate(parts[1:3], 1):
                            allele = allele.strip()
                            if allele and allele not in ['-', 'NA']:
                                normalized = normalize_allele(allele, resolution)
                                if normalized and normalized not in seen_alleles:
                                    reads = 0
                                    quality = 1.0

                                    # Look for reads column
                                    reads_col = f'reads{i}' if f'reads{i}' in header else 'reads'
                                    if reads_col in header:
                                        try:
                                            reads = int(parts[header.index(reads_col)])
                                        except:
                                            pass

                                    # Look for quality column
                                    if 'quality' in header or 'confidence' in header:
                                        qual_col = 'quality' if 'quality' in header else 'confidence'
                                        try:
                                            quality = float(parts[header.index(qual_col)])
                                        except:
                                            pass

                                    results[gene].append(AlleleCall(
                                        allele=normalized,
                                        reads=reads,
                                        quality=quality,
                                        tool='spechla'
                                    ))
                                    seen_alleles.add(normalized)
                    else:
                        # No header - assume standard format
                        for i, allele in enumerate(parts[1:3]):
                            allele = allele.strip()
                            if allele and allele not in ['-', 'NA']:
                                normalized = normalize_allele(allele, resolution)
                                if normalized and normalized not in seen_alleles:
                                    reads = 0
                                    # Try to get reads from columns 3,4 or 4,5
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
                        if allele and allele not in ['-', '?']:
                            normalized = normalize_allele(allele, resolution)
                            if normalized and normalized not in seen_alleles:
                                quality = 1.0
                                reads = 0

                                # Try to get quality from additional columns
                                if header and 'quality' in header:
                                    try:
                                        quality = float(parts[header.index('quality')])
                                    except:
                                        pass
                                elif len(parts) > 3:
                                    try:
                                        quality = float(parts[3])
                                    except:
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


def parse_results(tool: str, filepath: str, resolution: str = '2-field') -> Dict[str, List[AlleleCall]]:
    """Parse results based on tool type."""
    parsers = {
        'spechla': parse_spechla_results,
        'hlahd': parse_hlahd_results,
        'hlala': parse_hlala_results,
        'arcashla': parse_arcashla_results,
        'optitype': parse_optitype_results,
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


def weighted_consensus_vote(
    allele_data: Dict[str, List[AlleleCall]],
    min_tools: int = 1,
    expected_reads: int = 1000,
    weighting: str = 'read_confidence'
) -> Tuple[Optional[str], Optional[str], float, Dict[str, Any]]:
    """
    Perform weighted consensus voting on allele calls from different tools.

    Weighting methods:
    - 'equal': Each tool gets equal weight (traditional majority voting)
    - 'read_confidence': Weight by read support confidence
    - 'tool_quality': Weight by tool-reported quality scores

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
            args.weighting
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

    print(f"Consensus results written to {args.output_consensus}")
    print(f"Comparison results written to {args.output_comparison}")
    print(f"Weighting method used: {args.weighting}")


if __name__ == '__main__':
    main()
