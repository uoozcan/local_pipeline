#!/usr/bin/env python3
"""
calibrate_tool_weights.py

Empirically calibrate HLA typing tool weights using 1000 Genomes Project (1KGP)
ground-truth data with IPD-IMGT/HLA validated types.

For each tool, per-gene 2-field concordance is computed against the ground truth,
then normalized within each gene so tool weights sum to unity. Data-type-specific
weights are calculated separately for WGS and RNA-seq inputs.

Usage
-----
# Step 1: download 1KGP ground truth
python calibrate_tool_weights.py download-gt --output 1kgp_hla_gt.tsv

# Step 2: run calibration (Class I only — A, B, C)
python calibrate_tool_weights.py calibrate \\
    --ground-truth 1kgp_hla_gt.tsv \\
    --results-dir /path/to/tool_results/ \\
    --data-type wgs \\
    --genes A,B,C \\
    --resolution 2-field \\
    --output-weights ../conf/tool_weights_wgs.json \\
    --output-table ../conf/tool_accuracy_wgs.tsv

Expected results-dir layout
----------------------------
results-dir/
  hlahd/     -> {sample}_hlahd.txt   (one file per 1KGP sample)
  spechla/   -> {sample}_spechla.txt
  arcashla/  -> {sample}_arcashla.txt
  optitype/  -> {sample}_optitype.txt
  hlala/     -> {sample}_hlala.txt
  xhla/      -> {sample}_xhla.txt

Ground-truth TSV format (produced by download-gt or supplied manually)
-----------------------------------------------------------------------
Sample  A_1      A_2      B_1      B_2      C_1     C_2     DRB1_1    DRB1_2    DQB1_1    DQB1_2
NA12878 A*01:01  A*03:01  B*07:02  B*08:01  C*07:01 C*07:02 DRB1*03:01 DRB1*07:01 DQB1*02:01 DQB1*03:03

References
----------
Gourraud et al. (2014) PLoS ONE  — 1KGP Phase 1 HLA types (1,070 samples)
Abi-Rached et al. (2018)         — Extended loci (DPA1, DPB1, DQA1)
"""

import argparse
import json
import logging
import os
import re
import shlex
import sys
import urllib.request
from collections import defaultdict
from math import log
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Reuse parsers from consensus_voting.py (same directory)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(_SCRIPT_DIR))

from consensus_voting import (
    AlleleCall,
    normalize_allele,
    parse_arcashla_results,
    parse_hifihla_results,
    parse_hlahd_results,
    parse_hlala_results,
    parse_optitype_results,
    parse_spechla_results,
    parse_t1k_results,
    parse_xhla_results,
    weighted_consensus_vote,
    load_calibrated_weights,
)

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Canonical gene order for output tables
CLASSICAL_GENES = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']

# Gene coverage per tool: which genes each tool is expected to type
TOOL_GENE_COVERAGE: Dict[str, Set[str]] = {
    'hlahd':    set(CLASSICAL_GENES),
    'spechla':  set(CLASSICAL_GENES),
    'arcashla': set(CLASSICAL_GENES),
    'hlala':    {'A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'},
    'optitype': {'A', 'B', 'C'},
    'xhla':     {'A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1'},   # no DQA1/DPA1
    't1k':      set(CLASSICAL_GENES),                        # full Class I + II
    'hifihla':  set(CLASSICAL_GENES),                        # long-read, 4-field
}

# Parser dispatch (mirrors consensus_voting.py)
TOOL_PARSERS = {
    'hlahd':    parse_hlahd_results,
    'spechla':  parse_spechla_results,
    'arcashla': parse_arcashla_results,
    'optitype': parse_optitype_results,
    'hlala':    parse_hlala_results,
    'xhla':     parse_xhla_results,
    't1k':      parse_t1k_results,
    'hifihla':  parse_hifihla_results,
}

# 1KGP integrated call set panel (sample → population → superpopulation)
PANEL_URL = (
    'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
    'integrated_call_samples_v3.20130502.ALL.panel'
)

# Superpopulation codes
SUPERPOPULATIONS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

# 1KGP phase-1 HLA ground truth (Gourraud et al. 2014 supplementary data S1)
GT_URL = (
    'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/'
    '20140725_hla_genotypes/20140702_hla_diversity.txt'
)


# ---------------------------------------------------------------------------
# Ground-truth handling
# ---------------------------------------------------------------------------

def download_1kgp_ground_truth(output_path: str) -> None:
    """
    Download and reformat the 1KGP HLA ground-truth file (Gourraud et al. 2014).

    The raw 1KGP file encodes alleles as bare field pairs (e.g. '01:01') without
    the gene prefix.  This function adds the standard prefix (e.g. 'A*01:01') and
    writes a tab-separated file with one header line:

        Sample  A_1  A_2  B_1  B_2  C_1  C_2  DRB1_1  DRB1_2  DQB1_1  DQB1_2

    The 2014 file covers: A, B, C, DRB1, DQB1 for 1,070 samples.
    """
    logger.info(f"Downloading 1KGP ground truth from {GT_URL} ...")
    try:
        with urllib.request.urlopen(GT_URL, timeout=120) as resp:
            raw = resp.read().decode('utf-8')
    except Exception as e:
        logger.error(f"Download failed: {e}")
        logger.info("Please download manually and supply with --ground-truth")
        sys.exit(1)

    lines = raw.strip().split('\n')
    if not lines:
        logger.error("Downloaded file is empty")
        sys.exit(1)

    # The 1KGP Phase 1 file (20140702_hla_diversity.txt) is space-separated with
    # quoted values.  The header line is just "Sample" with no allele column names.
    # Allele columns are positional (0-indexed after shlex split):
    #   0: Sample   1: Population
    #   2: A_1  3: A_2   4: B_1  5: B_2   6: C_1  7: C_2
    #   8: DRB1_1  9: DRB1_2  10: DQB1_1  11: DQB1_2
    GT_GENE_COLS = {
        'A':    [2, 3],
        'B':    [4, 5],
        'C':    [6, 7],
        'DRB1': [8, 9],
        'DQB1': [10, 11],
    }
    out_genes = [g for g in ['A', 'B', 'C', 'DRB1', 'DQB1'] if g in CLASSICAL_GENES]
    out_header = ['Sample'] + [f'{g}_{i+1}' for g in out_genes for i in range(2)]
    logger.info(f"Using positional column mapping: { {g: GT_GENE_COLS[g] for g in out_genes} }")

    rows = ['\t'.join(out_header)]
    for line in lines[1:]:
        if not line.strip():
            continue
        try:
            parts = shlex.split(line)
        except ValueError:
            # Fall back to simple split if shlex fails (malformed line)
            parts = line.split()
        if len(parts) < 2:
            continue
        sample = parts[0].strip()
        row = [sample]
        for gene in out_genes:
            for col_idx in GT_GENE_COLS[gene]:
                raw_allele = parts[col_idx].strip() if col_idx < len(parts) else ''
                # Handle ambiguous alleles like '01:01/01:02' — take first option
                raw_allele = raw_allele.split('/')[0].strip()
                if raw_allele in ('', '-', 'NA', 'Not typed', 'NT'):
                    row.append('-')
                else:
                    # Add gene prefix if missing (raw file has bare fields like '01:01')
                    if '*' not in raw_allele:
                        row.append(f"{gene}*{raw_allele}")
                    else:
                        row.append(raw_allele)
        rows.append('\t'.join(row))

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(rows) + '\n')
    logger.info(f"Ground truth written to {output_path}  ({len(rows)-1} samples)")


def load_ground_truth(filepath: str, resolution: str = '2-field') -> Dict[str, Dict[str, List[str]]]:
    """
    Load ground-truth HLA types from a TSV or CSV file.

    Supports two header naming conventions:
    - Underscore:  Sample, A_1, A_2, B_1, B_2, DRB1_1, DRB1_2, ...  (legacy TSV)
    - No separator: sample, A1, A2, B1, B2, C1, C2, ...              (new CSV format)

    Handles duplicate sample IDs: exact duplicates are collapsed; conflicting
    duplicates keep the first occurrence and emit a warning.

    Returns:
        {sample_id: {gene_short: [allele1, allele2]}}
        e.g. {'NA12878': {'A': ['A*01:01', 'A*11:01'], 'B': [...], ...}}
    """
    # Auto-detect separator
    with open(filepath, 'r') as f:
        first_line = f.readline()
    sep = ',' if ',' in first_line and first_line.count(',') >= first_line.count('\t') else '\t'

    gt: Dict[str, Dict[str, List[str]]] = {}
    conflicts: List[str] = []

    with open(filepath, 'r') as f:
        header: Optional[List[str]] = None
        col_map: Dict[int, str] = {}   # col_index → gene_short

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(sep)]

            if header is None:
                header = parts
                # Build col_map: parse both 'A_1'/'DRB1_2' and 'A1'/'B2'/'C1' style
                for idx, col in enumerate(header[1:], start=1):
                    # Pattern 1: GENE_N  e.g. A_1, DRB1_2
                    m = re.match(r'^([A-Z0-9]+)_\d+$', col)
                    if m:
                        col_map[idx] = m.group(1)
                        continue
                    # Pattern 2: GENE+digit  e.g. A1, A2, B1, DRB11, DQB12
                    # Greedy match: longest known gene name
                    m2 = re.match(r'^(DRB1|DQA1|DQB1|DPA1|DPB1|[A-Z])(\d+)$', col)
                    if m2:
                        col_map[idx] = m2.group(1)
                continue

            sample = parts[0].strip()
            if not sample:
                continue

            row_gt: Dict[str, List[str]] = defaultdict(list)
            for col_idx, gene in col_map.items():
                raw = parts[col_idx] if col_idx < len(parts) else ''
                raw = raw.strip()
                if raw in ('-', 'NA', 'Not typed', 'NT', ''):
                    continue
                # Strip ambiguous notation e.g. A*01:01/A*01:02 → take first
                raw = raw.split('/')[0].strip()
                # Add gene prefix if missing (bare allele: '01:01' → 'A*01:01')
                if '*' not in raw:
                    raw = f"{gene}*{raw}"
                normalized = normalize_allele(raw, resolution)
                if normalized:
                    row_gt[gene].append(normalized)

            if sample in gt:
                # Duplicate: check if data is identical
                if dict(gt[sample]) != dict(row_gt):
                    conflicts.append(sample)
                # Keep first occurrence regardless
            else:
                gt[sample] = row_gt  # type: ignore[assignment]

    if conflicts:
        unique_conflicts = sorted(set(conflicts))
        logger.warning(
            f"Conflicting duplicate entries for {len(unique_conflicts)} samples "
            f"(keeping first occurrence): {unique_conflicts[:5]}{'...' if len(unique_conflicts) > 5 else ''}"
        )
    logger.info(f"Loaded ground truth for {len(gt)} samples "
                f"(sep={repr(sep)}, genes={sorted(set(g for s in gt.values() for g in s))})")
    return gt


# ---------------------------------------------------------------------------
# Allele comparison
# ---------------------------------------------------------------------------

def allele_set_match(
    tool_calls: List[AlleleCall],
    truth_alleles: List[str],
) -> Tuple[int, int]:
    """
    Compare a tool's allele calls for one gene to the ground truth.

    A locus is considered *correct* when the set of unique called alleles
    (at 2-field resolution) equals the set of ground-truth alleles (order-
    independent).  For homozygous loci (one unique truth allele) a single
    matching call is sufficient.

    Returns:
        (correct_alleles, total_alleles)
        where total_alleles = len(truth_alleles) (1 or 2)
    """
    if not truth_alleles:
        return 0, 0

    called = {c.allele for c in tool_calls}
    truth_set = set(truth_alleles)

    correct = len(called & truth_set)
    return correct, len(truth_set)


def locus_correct(tool_calls: List[AlleleCall], truth_alleles: List[str]) -> bool:
    """Return True if the tool's genotype exactly matches the ground truth (unordered)."""
    if not truth_alleles or not tool_calls:
        return False
    called_set = {c.allele for c in tool_calls}
    truth_set = set(truth_alleles)
    return called_set == truth_set


# ---------------------------------------------------------------------------
# Core calibration
# ---------------------------------------------------------------------------

def calibrate(
    ground_truth: Dict[str, Dict[str, List[str]]],
    results_dir: str,
    data_type: str,
    resolution: str = '2-field',
    genes: Optional[List[str]] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Compute per-tool, per-gene accuracy against 1KGP ground truth.

    Scans `results_dir` for subdirectories named after each tool and
    looks for files matching `{sample}_{tool}.txt`.

    Args:
        genes: List of gene names to evaluate (e.g. ['A','B','C']).
               None means all classical genes.

    Returns:
        {tool: {gene: concordance_rate}}
        where concordance_rate is correct loci / typed loci (0.0–1.0)
    """
    eval_genes: Set[str] = set(genes) if genes else set(CLASSICAL_GENES)
    logger.info(f"Evaluating genes: {sorted(eval_genes)}")

    results_dir_path = Path(results_dir)

    # Discover available tools from subdirectories
    available_tools = [
        t for t in TOOL_PARSERS
        if (results_dir_path / t).is_dir()
    ]

    if not available_tools:
        logger.error(
            f"No tool subdirectories found in {results_dir}. "
            f"Expected subdirs named: {list(TOOL_PARSERS.keys())}"
        )
        sys.exit(1)

    logger.info(f"Tools found: {available_tools}")

    # Structures: correct[tool][gene], total[tool][gene]
    correct: Dict[str, Dict[str, int]] = {t: defaultdict(int) for t in available_tools}
    total:   Dict[str, Dict[str, int]] = {t: defaultdict(int) for t in available_tools}

    for sample, truth_genes in ground_truth.items():
        for tool in available_tools:
            result_file = results_dir_path / tool / f"{sample}_{tool}.txt"
            if not result_file.exists():
                continue

            parser = TOOL_PARSERS[tool]
            tool_results = parser(str(result_file), resolution)  # {HLA-GENE: [AlleleCall]}

            for gene_short, truth_alleles in truth_genes.items():
                if not truth_alleles:
                    continue

                # Skip genes not in the requested evaluation set
                if gene_short not in eval_genes:
                    continue

                # Skip genes not supported by this tool
                if gene_short not in TOOL_GENE_COVERAGE.get(tool, set()):
                    continue

                hla_gene = f"HLA-{gene_short}"
                tool_calls = tool_results.get(hla_gene, [])

                total[tool][gene_short] += 1
                if locus_correct(tool_calls, truth_alleles):
                    correct[tool][gene_short] += 1

    # Compute concordance rates
    accuracy: Dict[str, Dict[str, float]] = {}
    for tool in available_tools:
        accuracy[tool] = {}
        for gene in CLASSICAL_GENES:
            n = total[tool].get(gene, 0)
            c = correct[tool].get(gene, 0)
            if n > 0:
                accuracy[tool][gene] = c / n
            else:
                # Gene not evaluated (not in GT or not typed by tool)
                accuracy[tool][gene] = None  # type: ignore[assignment]

        typed_genes = [g for g in CLASSICAL_GENES if accuracy[tool][g] is not None]
        if typed_genes:
            mean_acc = mean(accuracy[tool][g] for g in typed_genes)  # type: ignore[misc]
            logger.info(
                f"  {tool:12s}  mean concordance = {mean_acc:.3f}  "
                f"(n_genes={len(typed_genes)}, "
                f"n_samples={max(total[tool].values(), default=0)})"
            )

    return accuracy


def derive_weights(
    accuracy: Dict[str, Dict[str, float]],
    n_samples: Dict[str, int],
    data_type: str,
    resolution: str,
    genes: Optional[List[str]] = None,
) -> Dict[str, Dict[str, float]]:
    """
    Compute Method B (per-gene normalized) weights from per-tool accuracy.

    For each gene, weights across tools that cover it sum to 1.0.
    Tools that do not cover a gene (or have no accuracy data) get weight 0.0.
    Genes not in `genes` are assigned equal placeholder weights so that
    consensus_voting.py can still fall back gracefully.

    Returns:
        {gene: {tool: weight}}  (gene-first for easy lookup in consensus voting)
    """
    eval_genes: Set[str] = set(genes) if genes else set(CLASSICAL_GENES)
    tools = list(accuracy.keys())
    weights: Dict[str, Dict[str, float]] = {gene: {} for gene in CLASSICAL_GENES}

    for gene in CLASSICAL_GENES:
        # Genes outside the evaluation set: equal weight among covering tools
        if gene not in eval_genes:
            covering = [t for t in tools if gene in TOOL_GENE_COVERAGE.get(t, set())]
            eq = 1.0 / len(covering) if covering else 0.0
            for t in tools:
                weights[gene][t] = eq if t in covering else 0.0
            continue
        # Only include tools that (a) cover the gene and (b) have accuracy data
        covering = [
            t for t in tools
            if gene in TOOL_GENE_COVERAGE.get(t, set())
            and accuracy[t].get(gene) is not None
        ]
        if not covering:
            for t in tools:
                weights[gene][t] = 0.0
            continue

        gene_accs = {t: accuracy[t][gene] for t in covering}  # type: ignore[index]
        total_acc = sum(gene_accs.values())

        if total_acc == 0:
            # All tools failed on this gene — fall back to equal weights
            eq = 1.0 / len(covering)
            for t in tools:
                weights[gene][t] = eq if t in covering else 0.0
        else:
            for t in tools:
                if t in covering:
                    weights[gene][t] = gene_accs[t] / total_acc
                else:
                    weights[gene][t] = 0.0

    return weights


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_weights_json(
    weights: Dict[str, Dict[str, float]],
    accuracy: Dict[str, Dict[str, float]],
    n_samples: int,
    data_type: str,
    resolution: str,
    output_path: str,
) -> None:
    """Write calibrated weights to a JSON file consumed by consensus_voting.py."""
    tools = list(accuracy.keys())
    payload = {
        'method': 'normalized_concordance',
        'resolution': resolution,
        'data_type': data_type,
        'n_samples': n_samples,
        'description': (
            'Per-gene normalized 2-field concordance weights derived from '
            '1000 Genomes Project ground-truth HLA types (Gourraud et al. 2014). '
            'Weights per gene sum to 1.0 across tools that cover that gene.'
        ),
        'genes': {
            gene: {t: round(weights[gene].get(t, 0.0), 4) for t in tools}
            for gene in CLASSICAL_GENES
        },
        'raw_accuracy': {
            tool: {
                gene: round(accuracy[tool].get(gene, 0.0) or 0.0, 4)
                for gene in CLASSICAL_GENES
            }
            for tool in tools
        },
    }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(payload, f, indent=2)
    logger.info(f"Weights written to {output_path}")


def write_accuracy_table(
    accuracy: Dict[str, Dict[str, float]],
    weights: Dict[str, Dict[str, float]],
    total_counts: Dict[str, int],
    data_type: str,
    output_path: str,
    genes: Optional[List[str]] = None,
) -> None:
    """
    Write a human-readable TSV accuracy table for the manuscript supplementary.

    Columns: Tool | Data_type | N_samples | <genes...> | Mean_concordance
    Only the evaluated genes appear as columns (others are not measured).
    """
    report_genes = genes if genes else CLASSICAL_GENES
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        gene_cols = '\t'.join(report_genes)
        f.write(f"Tool\tData_type\tN_samples\t{gene_cols}\tMean_concordance\n")
        for tool, gene_acc in sorted(accuracy.items()):
            n = total_counts.get(tool, 0)
            gene_vals = []
            for gene in report_genes:
                v = gene_acc.get(gene)
                gene_vals.append(f"{v:.3f}" if v is not None else 'NA')
            typed = [gene_acc[g] for g in report_genes if gene_acc.get(g) is not None]
            mean_conc = f"{mean(typed):.3f}" if typed else 'NA'  # type: ignore[misc]
            f.write(f"{tool}\t{data_type}\t{n}\t" + '\t'.join(gene_vals) + f"\t{mean_conc}\n")
    logger.info(f"Accuracy table written to {output_path}")


def write_weights_summary(
    weights: Dict[str, Dict[str, float]],
    accuracy: Dict[str, Dict[str, float]],
    genes: Optional[List[str]] = None,
) -> None:
    """Print a compact summary table to stdout."""
    report_genes = genes if genes else CLASSICAL_GENES
    tools = list(accuracy.keys())
    print("\n=== Per-gene normalized weights (Method B) ===")
    header = f"{'Gene':<8}" + ''.join(f"{t:>12}" for t in tools)
    print(header)
    print('-' * len(header))
    for gene in report_genes:
        row = f"{gene:<8}" + ''.join(
            f"{weights[gene].get(t, 0.0):>12.3f}" for t in tools
        )
        print(row)

    print("\n=== Raw concordance rates ===")
    print(header)
    print('-' * len(header))
    for gene in report_genes:
        row = f"{gene:<8}" + ''.join(
            f"{accuracy[t].get(gene) or 0.0:>12.3f}" for t in tools
        )
        print(row)


def write_progress_summary(
    accuracy: Dict[str, Dict[str, float]],
    total_counts: Dict[str, int],
    genes: Optional[List[str]] = None,
    target_n: int = 50,
) -> None:
    """Print a progress table showing sample counts and concordance per tool toward target."""
    report_genes = genes if genes else CLASSICAL_GENES
    tools = list(accuracy.keys())

    print(f"\n=== Progress toward {target_n}-sample target ===")
    gene_header = ''.join(f"{g:>10}" for g in report_genes)
    header = f"{'Tool':<12}  {'N':>6}  {gene_header}  {'Remaining':>10}"
    print(header)
    print('-' * len(header))

    for tool in sorted(tools):
        n = total_counts.get(tool, 0)
        remaining = max(0, target_n - n)
        gene_vals = ''
        for gene in report_genes:
            v = accuracy[tool].get(gene)
            gene_vals += f"{v:>10.3f}" if v is not None else f"{'NA':>10}"
        print(f"{tool:<12}  {n:>6}  {gene_vals}  {remaining:>10}")


# ---------------------------------------------------------------------------
# Population file helpers
# ---------------------------------------------------------------------------

def download_population_file(output_path: str) -> None:
    """
    Download the 1KGP integrated call set panel file and reformat as
    a two-column TSV: sample_id → superpopulation (AFR/AMR/EAS/EUR/SAS).

    Source panel format (tab-separated):
        sample  pop  super_pop  gender
        HG00096 GBR  EUR        male
    """
    logger.info(f"Downloading 1KGP population panel from {PANEL_URL} ...")
    try:
        with urllib.request.urlopen(PANEL_URL, timeout=60) as resp:
            raw = resp.read().decode('utf-8')
    except Exception as e:
        logger.error(f"Download failed: {e}")
        logger.info("Please download manually from:\n  " + PANEL_URL)
        sys.exit(1)

    lines = raw.strip().split('\n')
    rows = ['sample\tsuperpopulation\tpopulation']
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        sample, pop, super_pop = parts[0], parts[1], parts[2]
        rows.append(f"{sample}\t{super_pop}\t{pop}")

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write('\n'.join(rows) + '\n')
    logger.info(f"Population file written to {output_path}  ({len(rows)-1} samples)")


def load_population_file(filepath: str) -> Dict[str, str]:
    """
    Load a sample → superpopulation mapping from a TSV file.

    Accepts both the output of download_population_file() and the raw 1KGP panel:
        sample  superpopulation  [population]

    Returns {sample_id: superpopulation_code}
    """
    pop_map: Dict[str, str] = {}
    with open(filepath, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if header is None:
                header = [p.lower() for p in parts]
                continue
            if len(parts) < 2:
                continue
            sample = parts[0].strip()
            # second column is superpopulation in our format, or population in raw panel
            # If raw panel: parts[2] is super_pop; our format: parts[1] is superpopulation
            if len(parts) >= 3 and parts[2].strip() in SUPERPOPULATIONS:
                # Raw panel format: sample, pop, super_pop, gender
                superpop = parts[2].strip()
            else:
                superpop = parts[1].strip()
            if superpop in SUPERPOPULATIONS:
                pop_map[sample] = superpop
    logger.info(f"Loaded population labels for {len(pop_map)} samples "
                f"(superpopulations: {sorted(set(pop_map.values()))})")
    return pop_map


def write_per_population_tables(
    per_pop_accuracy: Dict[str, Dict[str, Dict[str, float]]],
    per_pop_counts: Dict[str, Dict[str, int]],
    data_type: str,
    output_dir: str,
    genes: Optional[List[str]] = None,
) -> None:
    """
    Write one accuracy TSV per superpopulation for supplementary tables.
    Also writes a summary table with all superpopulations side-by-side.
    """
    if not per_pop_accuracy:
        return
    report_genes = genes if genes else CLASSICAL_GENES
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    all_tools = sorted({t for pop_acc in per_pop_accuracy.values() for t in pop_acc})

    # Per-population files
    for superpop, accuracy in per_pop_accuracy.items():
        out_path = str(Path(output_dir) / f"tool_accuracy_{data_type}_{superpop}.tsv")
        n_samples = max(per_pop_counts.get(superpop, {}).values(), default=0)
        with open(out_path, 'w') as f:
            gene_cols = '\t'.join(report_genes)
            f.write(f"Tool\tData_type\tSuperpopulation\tN_samples\t{gene_cols}\tMean_concordance\n")
            for tool in sorted(accuracy):
                tool_acc = accuracy[tool]
                n = per_pop_counts.get(superpop, {}).get(tool, 0)
                gene_vals = [f"{tool_acc.get(g) or 0.0:.3f}" if tool_acc.get(g) is not None
                             else 'NA' for g in report_genes]
                typed = [tool_acc[g] for g in report_genes if tool_acc.get(g) is not None]
                mean_c = f"{mean(typed):.3f}" if typed else 'NA'  # type: ignore[misc]
                f.write(f"{tool}\t{data_type}\t{superpop}\t{n}\t"
                        + '\t'.join(gene_vals) + f"\t{mean_c}\n")
        logger.info(f"Per-population accuracy written to {out_path}")

    # Combined summary: rows = tool×superpop, cols = genes
    summary_path = str(Path(output_dir) / f"tool_accuracy_{data_type}_by_population.tsv")
    pop_list = sorted(per_pop_accuracy.keys())
    with open(summary_path, 'w') as f:
        gene_cols = '\t'.join(
            f"{pop}_{gene}" for pop in pop_list for gene in report_genes
        )
        f.write(f"Tool\t{gene_cols}\n")
        for tool in all_tools:
            vals = []
            for pop in pop_list:
                pop_acc = per_pop_accuracy.get(pop, {}).get(tool, {})
                for gene in report_genes:
                    v = pop_acc.get(gene)
                    vals.append(f"{v:.3f}" if v is not None else 'NA')
            f.write(f"{tool}\t" + '\t'.join(vals) + "\n")
    logger.info(f"Population summary written to {summary_path}")


# ---------------------------------------------------------------------------
# Population-stratified calibration (extends calibrate())
# ---------------------------------------------------------------------------

def calibrate_stratified(
    ground_truth: Dict[str, Dict[str, List[str]]],
    results_dir: str,
    data_type: str,
    resolution: str = '2-field',
    genes: Optional[List[str]] = None,
    population_map: Optional[Dict[str, str]] = None,
) -> Tuple[
    Dict[str, Dict[str, float]],                   # overall accuracy
    Dict[str, Dict[str, Dict[str, float]]],         # per_pop_accuracy[pop][tool][gene]
    Dict[str, int],                                 # total_counts[tool]
    Dict[str, Dict[str, int]],                      # per_pop_counts[pop][tool]
]:
    """
    Extends calibrate() to simultaneously compute overall and per-superpopulation accuracy.
    """
    eval_genes: Set[str] = set(genes) if genes else set(CLASSICAL_GENES)
    results_dir_path = Path(results_dir)

    available_tools = [t for t in TOOL_PARSERS if (results_dir_path / t).is_dir()]
    if not available_tools:
        logger.error(f"No tool subdirectories found in {results_dir}")
        sys.exit(1)
    logger.info(f"Tools found: {available_tools}")

    # Overall tallies
    correct: Dict[str, Dict[str, int]] = {t: defaultdict(int) for t in available_tools}
    total:   Dict[str, Dict[str, int]] = {t: defaultdict(int) for t in available_tools}

    # Per-population tallies
    pop_correct: Dict[str, Dict[str, Dict[str, int]]] = {}
    pop_total:   Dict[str, Dict[str, Dict[str, int]]] = {}

    for sample, truth_genes in ground_truth.items():
        superpop = population_map.get(sample) if population_map else None

        for tool in available_tools:
            result_file = results_dir_path / tool / f"{sample}_{tool}.txt"
            if not result_file.exists():
                continue

            tool_results = TOOL_PARSERS[tool](str(result_file), resolution)

            for gene_short, truth_alleles in truth_genes.items():
                if not truth_alleles or gene_short not in eval_genes:
                    continue
                if gene_short not in TOOL_GENE_COVERAGE.get(tool, set()):
                    continue

                hla_gene = f"HLA-{gene_short}"
                tool_calls = tool_results.get(hla_gene, [])

                total[tool][gene_short] += 1
                is_correct = locus_correct(tool_calls, truth_alleles)
                if is_correct:
                    correct[tool][gene_short] += 1

                # Population-stratified
                if superpop:
                    if superpop not in pop_correct:
                        pop_correct[superpop] = {t: defaultdict(int) for t in available_tools}
                        pop_total[superpop] = {t: defaultdict(int) for t in available_tools}
                    pop_total[superpop][tool][gene_short] += 1
                    if is_correct:
                        pop_correct[superpop][tool][gene_short] += 1

    # Compute overall accuracy rates
    accuracy: Dict[str, Dict[str, float]] = {}
    for tool in available_tools:
        accuracy[tool] = {}
        for gene in CLASSICAL_GENES:
            n = total[tool].get(gene, 0)
            c = correct[tool].get(gene, 0)
            accuracy[tool][gene] = (c / n) if n > 0 else None  # type: ignore

        typed_genes = [g for g in CLASSICAL_GENES if accuracy[tool][g] is not None]
        if typed_genes:
            m = mean(accuracy[tool][g] for g in typed_genes)  # type: ignore[misc]
            logger.info(
                f"  {tool:12s}  mean={m:.3f}  n_genes={len(typed_genes)} "
                f"n_samples={max(total[tool].values(), default=0)}"
            )

    # Per-population accuracy rates
    per_pop_accuracy: Dict[str, Dict[str, Dict[str, float]]] = {}
    per_pop_counts: Dict[str, Dict[str, int]] = {}
    for superpop in pop_correct:
        per_pop_accuracy[superpop] = {}
        per_pop_counts[superpop] = {}
        for tool in available_tools:
            per_pop_accuracy[superpop][tool] = {}
            per_pop_counts[superpop][tool] = max(pop_total[superpop][tool].values(), default=0)
            for gene in CLASSICAL_GENES:
                n = pop_total[superpop][tool].get(gene, 0)
                c = pop_correct[superpop][tool].get(gene, 0)
                per_pop_accuracy[superpop][tool][gene] = (c / n) if n > 0 else None  # type: ignore

    total_counts = {t: max(total[t].values(), default=0) for t in available_tools}
    return accuracy, per_pop_accuracy, total_counts, per_pop_counts


# ---------------------------------------------------------------------------
# Voting strategy comparison
# ---------------------------------------------------------------------------

def compare_voting_strategies(
    ground_truth: Dict[str, Dict[str, List[str]]],
    results_dir: str,
    resolution: str = '2-field',
    genes: Optional[List[str]] = None,
    weights_file: Optional[str] = None,
    expected_reads: int = 1000,
) -> Dict[str, Dict[str, List[float]]]:
    """
    For each sample in results_dir, compute per-gene concordance under each
    voting strategy (equal / read_confidence / calibrated) by running the
    weighted_consensus_vote() function directly.

    Returns:
        {strategy: {gene: [concordance_per_sample]}}
    """
    eval_genes: Set[str] = set(genes) if genes else set(CLASSICAL_GENES)
    results_dir_path = Path(results_dir)
    available_tools = [t for t in TOOL_PARSERS if (results_dir_path / t).is_dir()]

    calibrated_weights: Optional[Dict] = None
    modes = ['equal', 'read_confidence']
    if weights_file and Path(weights_file).exists():
        calibrated_weights = load_calibrated_weights(weights_file)
        modes.append('calibrated')
        logger.info(f"Loaded calibrated weights from {weights_file}")

    # {mode: {gene: [concordances]}}
    concordances: Dict[str, Dict[str, List[float]]] = {
        m: {g: [] for g in eval_genes} for m in modes
    }

    samples_processed = 0
    for sample, truth_genes in ground_truth.items():
        # Parse all available tool results for this sample
        all_results: Dict[str, Dict] = {}
        for tool in available_tools:
            rf = results_dir_path / tool / f"{sample}_{tool}.txt"
            if rf.exists():
                all_results[tool] = TOOL_PARSERS[tool](str(rf), resolution)

        if not all_results:
            continue
        samples_processed += 1

        for gene_short in eval_genes:
            if gene_short not in truth_genes or not truth_genes[gene_short]:
                continue
            hla_gene = f"HLA-{gene_short}"
            truth = truth_genes[gene_short]

            gene_data = {
                tool: res.get(hla_gene, [])
                for tool, res in all_results.items()
                if res.get(hla_gene)
            }
            if not gene_data:
                continue

            for mode in modes:
                cw = calibrated_weights if mode == 'calibrated' else None
                a1, a2, _conf, _ = weighted_consensus_vote(
                    gene_data, min_tools=1, expected_reads=expected_reads,
                    weighting=mode, calibrated_weights=cw, gene=hla_gene,
                )
                called = {a for a in [a1, a2] if a}
                truth_set = set(truth)
                n_correct = len(called & truth_set)
                concordance = n_correct / len(truth_set) if truth_set else 0.0
                concordances[mode][gene_short].append(concordance)

    logger.info(f"Strategy comparison: {samples_processed} samples, modes={modes}")
    return concordances


def compute_wilcoxon_stats(
    concordances: Dict[str, Dict[str, List[float]]],
    ref_mode: str = 'equal',
    test_mode: str = 'calibrated',
) -> Dict[str, Dict]:
    """
    Run one-sided Wilcoxon signed-rank tests (test_mode > ref_mode) per gene.

    Requires scipy; if not installed, p-values are omitted with a warning.
    Returns {gene: {p_value, mean_ref, mean_test, delta, n, note}}.
    """
    try:
        from scipy.stats import wilcoxon as _scipy_wilcoxon
        _has_scipy = True
    except ImportError:
        _has_scipy = False
        logger.warning(
            "scipy not installed — Wilcoxon p-values skipped. "
            "Install with: pip install scipy"
        )

    out: Dict[str, Dict] = {}
    for gene in concordances.get(ref_mode, {}):
        ref_v = concordances[ref_mode].get(gene, [])
        test_v = concordances.get(test_mode, {}).get(gene, [])
        n = min(len(ref_v), len(test_v))
        mean_ref  = mean(ref_v[:n])  if n > 0 else 0.0  # type: ignore[misc]
        mean_test = mean(test_v[:n]) if n > 0 else 0.0  # type: ignore[misc]
        delta = mean_test - mean_ref

        p_value = None
        note = ''
        if n < 5:
            note = f'too few samples (n={n}) for reliable Wilcoxon test'
        elif _has_scipy:
            diffs = [t - r for t, r in zip(test_v[:n], ref_v[:n])]
            if all(d == 0 for d in diffs):
                p_value = 1.0
                note = 'no differences between strategies'
            else:
                try:
                    _, p_value = _scipy_wilcoxon(
                        test_v[:n], ref_v[:n],
                        alternative='greater', zero_method='zsplit'
                    )
                except Exception as exc:
                    note = f'Wilcoxon failed: {exc}'

        out[gene] = {
            'p_value':   p_value,
            'mean_ref':  mean_ref,
            'mean_test': mean_test,
            'delta':     delta,
            'n':         n,
            'note':      note,
        }
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_download_gt(args: argparse.Namespace) -> None:
    download_1kgp_ground_truth(args.output)


def cmd_calibrate(args: argparse.Namespace) -> None:
    # Parse gene list
    genes: Optional[List[str]] = None
    if args.genes and args.genes.upper() != 'ALL':
        genes = [g.strip().upper() for g in args.genes.split(',') if g.strip()]
        logger.info(f"Restricting calibration to genes: {genes}")
    else:
        logger.info("Evaluating all classical genes")

    logger.info(f"Loading ground truth from {args.ground_truth} ...")
    ground_truth = load_ground_truth(args.ground_truth, args.resolution)

    # Optionally load population labels for stratified analysis
    population_map: Optional[Dict[str, str]] = None
    if getattr(args, 'population_file', None) and Path(args.population_file).exists():
        population_map = load_population_file(args.population_file)
        logger.info("Population-stratified analysis enabled")

    logger.info(f"Running calibration (data_type={args.data_type}) ...")
    accuracy, per_pop_accuracy, total_counts, per_pop_counts = calibrate_stratified(
        ground_truth, args.results_dir, args.data_type, args.resolution,
        genes, population_map
    )

    n_samples_total = max(total_counts.values(), default=0)
    weights = derive_weights(accuracy, total_counts, args.data_type, args.resolution, genes)

    write_weights_summary(weights, accuracy, genes)
    write_progress_summary(accuracy, total_counts, genes, target_n=args.target_n)

    if per_pop_accuracy:
        logger.info(f"Per-population accuracy computed for: {sorted(per_pop_accuracy.keys())}")

    if args.output_weights:
        write_weights_json(
            weights, accuracy, n_samples_total,
            args.data_type, args.resolution, args.output_weights
        )

    if args.output_table:
        write_accuracy_table(
            accuracy, weights, total_counts, args.data_type, args.output_table, genes
        )
        # Write per-population tables in the same directory
        if per_pop_accuracy:
            out_dir = str(Path(args.output_table).parent)
            write_per_population_tables(
                per_pop_accuracy, per_pop_counts, args.data_type, out_dir, genes
            )


def cmd_compare_strategies(args: argparse.Namespace) -> None:
    """Compare voting strategy concordances with optional Wilcoxon significance tests."""
    genes: Optional[List[str]] = None
    if args.genes and args.genes.upper() != 'ALL':
        genes = [g.strip().upper() for g in args.genes.split(',') if g.strip()]

    logger.info(f"Loading ground truth from {args.ground_truth} ...")
    ground_truth = load_ground_truth(args.ground_truth, args.resolution)

    logger.info("Running per-sample voting strategy comparison ...")
    concordances = compare_voting_strategies(
        ground_truth, args.results_dir, args.resolution, genes,
        getattr(args, 'weights_file', None), args.expected_reads,
    )

    modes = list(concordances.keys())
    report_genes = genes if genes else CLASSICAL_GENES
    ref_mode  = args.ref_mode  if args.ref_mode  in modes else modes[0]
    test_mode = args.test_mode if args.test_mode in modes else modes[-1]

    stats = compute_wilcoxon_stats(concordances, ref_mode, test_mode)

    # Print summary table
    print(f"\n=== Voting Strategy Comparison: {test_mode} vs {ref_mode} ===")
    col_w = max(len(ref_mode), len(test_mode), 10)
    header = (f"{'Gene':<8}  {'N':>4}  {ref_mode:>{col_w}}  "
              f"{test_mode:>{col_w}}  {'Delta':>8}  {'p-value':>10}  {'Sig':>4}")
    print(header)
    print('-' * len(header))

    for gene in report_genes:
        if gene not in stats:
            continue
        s = stats[gene]
        delta   = s['delta']
        p_val   = s['p_value']
        p_str   = f"{p_val:.4f}" if p_val is not None else 'N/A'
        sig_str = '*' if (p_val is not None and p_val < 0.05) else ''
        note    = f"  [{s['note']}]" if s.get('note') else ''
        print(f"{gene:<8}  {s['n']:>4}  {s['mean_ref']:>{col_w}.3f}  "
              f"{s['mean_test']:>{col_w}.3f}  {delta:>+8.3f}  {p_str:>10}  {sig_str:>4}{note}")

    if len(modes) > 2:
        # Also print all modes for context
        print(f"\n=== Mean concordance across all strategies ===")
        gene_header = ''.join(f"{g:>8}" for g in report_genes)
        print(f"{'Strategy':<20}  {gene_header}")
        for mode in modes:
            gene_vals = ''.join(
                f"{mean(concordances[mode].get(g, [0])):>8.3f}"  # type: ignore[misc]
                if concordances[mode].get(g) else f"{'NA':>8}"
                for g in report_genes
            )
            print(f"{mode:<20}  {gene_vals}")

    # Write output TSV
    if args.output:
        with open(args.output, 'w') as f:
            f.write(f"Gene\tN\t{ref_mode}_mean\t{test_mode}_mean\t"
                    f"Delta\tWilcoxon_p\tSignificant_p05\tNote\n")
            for gene in report_genes:
                if gene not in stats:
                    continue
                s = stats[gene]
                p_val = s['p_value']
                p_str = f"{p_val:.6f}" if p_val is not None else 'NA'
                sig   = 'yes' if (p_val is not None and p_val < 0.05) else 'no'
                f.write(f"{gene}\t{s['n']}\t{s['mean_ref']:.4f}\t{s['mean_test']:.4f}\t"
                        f"{s['delta']:+.4f}\t{p_str}\t{sig}\t{s.get('note','')}\n")
        logger.info(f"Strategy comparison written to {args.output}")


def cmd_generate_pop_file(args: argparse.Namespace) -> None:
    download_population_file(args.output)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description='Calibrate HLA tool weights from 1000 Genomes ground truth',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = p.add_subparsers(dest='command', required=True)

    # --- download-gt ---
    dl = sub.add_parser('download-gt', help='Download and reformat 1KGP HLA ground truth')
    dl.add_argument('--output', required=True, help='Output TSV path')
    dl.set_defaults(func=cmd_download_gt)

    # --- generate-population-file ---
    gp = sub.add_parser(
        'generate-population-file',
        help='Download 1KGP panel and write sample → superpopulation TSV',
    )
    gp.add_argument('--output', required=True,
                    help='Output TSV path (sample, superpopulation, population)')
    gp.set_defaults(func=cmd_generate_pop_file)

    # --- calibrate ---
    cal = sub.add_parser('calibrate', help='Run empirical calibration and derive weights')
    cal.add_argument('--ground-truth', required=True,
                     help='Ground-truth TSV (from download-gt or custom)')
    cal.add_argument('--results-dir', required=True,
                     help='Directory with per-tool subdirs of pipeline result files')
    cal.add_argument('--data-type', required=True, choices=['wgs', 'rna'],
                     help='Sequencing type (affects tool selection)')
    cal.add_argument('--resolution', default='2-field', choices=['2-field', '4-field'],
                     help='Allele resolution for comparison')
    cal.add_argument('--genes', default='A,B,C,DRB1,DQB1',
                     help='Comma-separated genes to evaluate. '
                          'Other genes receive equal placeholder weights. '
                          'Use ALL to evaluate all classical genes.')
    cal.add_argument('--population-file', default=None,
                     help='Sample → superpopulation TSV (from generate-population-file). '
                          'Enables population-stratified accuracy tables.')
    cal.add_argument('--output-weights',
                     help='Output JSON weights file (used by consensus_voting.py)')
    cal.add_argument('--output-table',
                     help='Output TSV accuracy table (for manuscript supplementary). '
                          'Per-population tables are written to the same directory.')
    cal.add_argument('--target-n', type=int, default=50,
                     help='Target sample count per tool for progress summary')
    cal.set_defaults(func=cmd_calibrate)

    # --- compare-strategies ---
    cmp = sub.add_parser(
        'compare-strategies',
        help='Compare voting strategy concordances (equal/read_confidence/calibrated) '
             'with Wilcoxon signed-rank tests',
    )
    cmp.add_argument('--ground-truth', required=True,
                     help='Ground-truth TSV (from download-gt or custom)')
    cmp.add_argument('--results-dir', required=True,
                     help='Directory with per-tool subdirs of pipeline result files')
    cmp.add_argument('--resolution', default='2-field', choices=['2-field', '4-field'],
                     help='Allele resolution for comparison')
    cmp.add_argument('--genes', default='A,B,C,DRB1,DQB1',
                     help='Comma-separated genes to compare. Use ALL for all classical.')
    cmp.add_argument('--weights-file', default=None,
                     help='Calibrated weights JSON (enables the "calibrated" mode). '
                          'Produced by the calibrate subcommand.')
    cmp.add_argument('--expected-reads', type=int, default=1000,
                     help='Expected reads for read_confidence weighting')
    cmp.add_argument('--ref-mode', default='equal',
                     choices=['equal', 'read_confidence', 'calibrated'],
                     help='Reference (baseline) voting strategy')
    cmp.add_argument('--test-mode', default='calibrated',
                     choices=['equal', 'read_confidence', 'calibrated'],
                     help='Test voting strategy to compare against reference')
    cmp.add_argument('--output', default=None,
                     help='Output TSV with per-gene statistics (optional)')
    cmp.set_defaults(func=cmd_compare_strategies)

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
