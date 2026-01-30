#!/bin/bash
#SBATCH --job-name=hla_array
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/hla_array_%A_%a.out
#SBATCH --error=logs/hla_array_%A_%a.err
# NOTE: Set account and array size with:
#   sbatch --account=YOUR_PROJECT_ID --array=1-N submit_hla_array.sh

#
# HLA Typing Pipeline - SLURM Array Job
# Processes multiple samples in parallel (one job per sample)
#
# Usage:
#   1. Create sample list: ls input/*.bam | xargs -n1 basename | sed 's/.bam$//' > sample_list.txt
#   2. Submit: sbatch --array=1-$(wc -l < sample_list.txt) submit_hla_array.sh
#
# Required files:
#   - sample_list.txt: One sample ID per line
#   - input/<sample_id>.bam: BAM files
#   OR
#   - input/<sample_id>_R1.fastq.gz and input/<sample_id>_R2.fastq.gz: FASTQ files
#

set -euo pipefail

# ============================================
# Configuration
# ============================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

PROJECT_ID="${SLURM_JOB_ACCOUNT:-${CSC_PROJECT:-}}"
if [[ -z "$PROJECT_ID" ]]; then
    echo "ERROR: No project ID found. Submit with: sbatch --account=YOUR_PROJECT_ID ..."
    exit 1
fi
BASE_DIR="/scratch/${PROJECT_ID}/${USER}/hla_analysis"

# Sample list file
SAMPLE_LIST="${SAMPLE_LIST:-${BASE_DIR}/sample_list.txt}"

# Input directories
BAM_DIR="${BAM_DIR:-${BASE_DIR}/input/bam}"
FASTQ_DIR="${FASTQ_DIR:-${BASE_DIR}/input/fastq}"

# Default parameters
TOOLS="${TOOLS:-spechla,hlahd}"
REFERENCE="${REFERENCE:-hg38}"
RESOLUTION="${RESOLUTION:-2-field}"
SEQ_TYPE="${SEQ_TYPE:-dna}"
OUTDIR="${OUTDIR:-${BASE_DIR}/results}"

# ============================================
# Get Sample for This Array Task
# ============================================

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: This script must be run as a SLURM array job"
    echo "Usage: sbatch --array=1-N submit_hla_array.sh"
    exit 1
fi

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    echo ""
    echo "Create it with:"
    echo "  ls ${BAM_DIR}/*.bam | xargs -n1 basename | sed 's/.bam\$//' > $SAMPLE_LIST"
    exit 1
fi

# Get the sample ID for this task
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# ============================================
# Setup Environment
# ============================================

echo "=============================================="
echo "HLA Typing Pipeline - Array Job"
echo "=============================================="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample ID: ${SAMPLE_ID}"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

# Load modules
module purge
module load nextflow/23.10.0 2>/dev/null || module load nextflow
module load singularity 2>/dev/null || true

# Set Singularity cache
export SINGULARITY_CACHEDIR="${BASE_DIR}/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
mkdir -p "$SINGULARITY_CACHEDIR"

# Create directories
mkdir -p "$OUTDIR"
mkdir -p "${BASE_DIR}/logs"

# ============================================
# Determine Input Type
# ============================================

BAM_FILE="${BAM_DIR}/${SAMPLE_ID}.bam"
FASTQ_R1="${FASTQ_DIR}/${SAMPLE_ID}_R1.fastq.gz"
FASTQ_R2="${FASTQ_DIR}/${SAMPLE_ID}_R2.fastq.gz"

# Alternative FASTQ naming
if [[ ! -f "$FASTQ_R1" ]]; then
    FASTQ_R1="${FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz"
    FASTQ_R2="${FASTQ_DIR}/${SAMPLE_ID}_2.fastq.gz"
fi

if [[ -f "$BAM_FILE" ]]; then
    INPUT_TYPE="bam"
    INPUT_FILE="$BAM_FILE"
    echo "Input type: BAM"
    echo "Input file: $BAM_FILE"

    # Check/create BAM index
    if [[ ! -f "${BAM_FILE}.bai" ]] && [[ ! -f "${BAM_FILE%.bam}.bai" ]]; then
        echo "Creating BAM index..."
        module load samtools 2>/dev/null || module load biokit
        samtools index "$BAM_FILE"
    fi

    INPUT_ARGS="--input_bam $BAM_FILE"

elif [[ -f "$FASTQ_R1" ]] && [[ -f "$FASTQ_R2" ]]; then
    INPUT_TYPE="fastq"
    echo "Input type: FASTQ"
    echo "Input R1: $FASTQ_R1"
    echo "Input R2: $FASTQ_R2"

    INPUT_ARGS="--input_fastq_1 $FASTQ_R1 --input_fastq_2 $FASTQ_R2"
else
    echo "ERROR: No input files found for sample: $SAMPLE_ID"
    echo "  Checked: $BAM_FILE"
    echo "  Checked: $FASTQ_R1"
    exit 1
fi

echo ""

# ============================================
# Run Nextflow Pipeline
# ============================================

echo "Starting HLA typing for sample: $SAMPLE_ID"
echo "=============================================="

# Use sample-specific work directory to avoid conflicts
WORK_DIR="${BASE_DIR}/work/${SAMPLE_ID}"
mkdir -p "$WORK_DIR"

cd "$BASE_DIR"

nextflow run "${PIPELINE_DIR}/main.nf" \
    $INPUT_ARGS \
    --outdir "$OUTDIR" \
    --tools "$TOOLS" \
    --reference "$REFERENCE" \
    --resolution "$RESOLUTION" \
    --seq_type "$SEQ_TYPE" \
    --max_cpus ${SLURM_CPUS_PER_TASK:-20} \
    --max_memory "${SLURM_MEM_PER_NODE:-64000}M" \
    -profile singularity \
    -work-dir "$WORK_DIR" \
    -c "${PIPELINE_DIR}/conf/puhti.config" \
    -resume

# ============================================
# Completion
# ============================================

echo ""
echo "=============================================="
echo "Sample completed: $SAMPLE_ID"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Results: ${OUTDIR}/${SAMPLE_ID}/"
echo ""

# Show consensus results
CONSENSUS_FILE="${OUTDIR}/${SAMPLE_ID}/${SAMPLE_ID}_consensus.txt"
if [[ -f "$CONSENSUS_FILE" ]]; then
    echo "Consensus HLA types:"
    cat "$CONSENSUS_FILE"
fi

# Clean up work directory to save space
if [[ "${CLEANUP_WORK:-false}" == "true" ]]; then
    echo ""
    echo "Cleaning up work directory..."
    rm -rf "$WORK_DIR"
fi
