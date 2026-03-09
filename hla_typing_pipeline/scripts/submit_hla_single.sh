#!/bin/bash
#SBATCH --job-name=hla_typing
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/hla_%j.out
#SBATCH --error=logs/hla_%j.err
# NOTE: Set account with: sbatch --account=YOUR_PROJECT_ID submit_hla_single.sh

#
# HLA Typing Pipeline - Single Sample SLURM Submission Script
# Usage: sbatch submit_hla_single.sh <input_file> [options]
#
# Examples:
#   sbatch submit_hla_single.sh sample.bam
#   sbatch submit_hla_single.sh sample.bam --tools spechla,hlahd,arcashla
#   sbatch submit_hla_single.sh sample_R1.fastq.gz --fastq2 sample_R2.fastq.gz
#

set -euo pipefail

# ============================================
# Configuration
# ============================================

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Project settings - auto-detect from SLURM or environment
PROJECT_ID="${SLURM_JOB_ACCOUNT:-${CSC_PROJECT:-}}"
if [[ -z "$PROJECT_ID" ]]; then
    echo "ERROR: No project ID found. Submit with: sbatch --account=YOUR_PROJECT_ID ..."
    exit 1
fi
BASE_DIR="/scratch/${PROJECT_ID}/${USER}/hla_analysis"

# Default parameters
TOOLS="hlahd,spechla,arcashla"
REFERENCE="hg38"
RESOLUTION="2-field"
SEQ_TYPE="dna"
OUTDIR="${BASE_DIR}/results"

# ============================================
# Parse Arguments
# ============================================

INPUT_FILE=""
FASTQ2=""
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --fastq2)
            FASTQ2="$2"
            shift 2
            ;;
        --tools)
            TOOLS="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        --seq_type)
            SEQ_TYPE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -*)
            EXTRA_ARGS="$EXTRA_ARGS $1 $2"
            shift 2
            ;;
        *)
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
            fi
            shift
            ;;
    esac
done

# Validate input
if [[ -z "$INPUT_FILE" ]]; then
    echo "ERROR: No input file specified"
    echo ""
    echo "Usage: sbatch submit_hla_single.sh <input_file> [options]"
    echo ""
    echo "Options:"
    echo "  --fastq2 <file>      R2 FASTQ file (required for FASTQ input)"
    echo "  --tools <list>       HLA tools to use (default: spechla,hlahd)"
    echo "  --reference <ref>    Reference genome: hg38 or hg19 (default: hg38)"
    echo "  --resolution <res>   Output resolution: 2-field or 4-field (default: 2-field)"
    echo "  --seq_type <type>    Sequence type for OptiType: dna or rna (default: dna)"
    echo "  --outdir <dir>       Output directory (default: ${BASE_DIR}/results)"
    exit 1
fi

# Make input path absolute
if [[ ! "$INPUT_FILE" = /* ]]; then
    INPUT_FILE="${SLURM_SUBMIT_DIR:-$PWD}/$INPUT_FILE"
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# ============================================
# Setup Environment
# ============================================

echo "=============================================="
echo "HLA Typing Pipeline - Puhti"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
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

# Create output and log directories
mkdir -p "$OUTDIR"
mkdir -p "${BASE_DIR}/logs"

echo "Configuration:"
echo "  Input: $INPUT_FILE"
echo "  Output: $OUTDIR"
echo "  Tools: $TOOLS"
echo "  Reference: $REFERENCE"
echo "  Resolution: $RESOLUTION"
echo ""

# ============================================
# Determine Input Type and Run Pipeline
# ============================================

cd "$BASE_DIR"

# Check if input is BAM or FASTQ
if [[ "$INPUT_FILE" == *.bam ]]; then
    echo "Input type: BAM"

    # Check for BAM index
    if [[ ! -f "${INPUT_FILE}.bai" ]] && [[ ! -f "${INPUT_FILE%.bam}.bai" ]]; then
        echo "Creating BAM index..."
        module load samtools 2>/dev/null || module load biokit
        samtools index "$INPUT_FILE"
    fi

    INPUT_ARGS="--input_bam $INPUT_FILE"

elif [[ "$INPUT_FILE" == *.fastq.gz ]] || [[ "$INPUT_FILE" == *.fq.gz ]] || [[ "$INPUT_FILE" == *.fastq ]] || [[ "$INPUT_FILE" == *.fq ]]; then
    echo "Input type: FASTQ"

    if [[ -z "$FASTQ2" ]]; then
        echo "ERROR: FASTQ input requires --fastq2 option"
        exit 1
    fi

    # Make FASTQ2 path absolute
    if [[ ! "$FASTQ2" = /* ]]; then
        FASTQ2="${SLURM_SUBMIT_DIR:-$PWD}/$FASTQ2"
    fi

    if [[ ! -f "$FASTQ2" ]]; then
        echo "ERROR: R2 FASTQ file not found: $FASTQ2"
        exit 1
    fi

    INPUT_ARGS="--input_fastq_1 $INPUT_FILE --input_fastq_2 $FASTQ2"
else
    echo "ERROR: Unrecognized input file type. Must be .bam or .fastq.gz"
    exit 1
fi

# ============================================
# Run Nextflow Pipeline
# ============================================

echo ""
echo "Starting HLA typing pipeline..."
echo "=============================================="

nextflow run "${PIPELINE_DIR}/main.nf" \
    $INPUT_ARGS \
    --outdir "$OUTDIR" \
    --tools "$TOOLS" \
    --reference "$REFERENCE" \
    --resolution "$RESOLUTION" \
    --seq_type "$SEQ_TYPE" \
    --max_cpus ${SLURM_CPUS_PER_TASK:-40} \
    --max_memory "${SLURM_MEM_PER_NODE:-180000}M" \
    -profile singularity \
    -c "${PIPELINE_DIR}/conf/puhti.config" \
    -resume \
    $EXTRA_ARGS

# ============================================
# Completion
# ============================================

echo ""
echo "=============================================="
echo "Pipeline completed!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Results available at: $OUTDIR"
echo ""
echo "Key output files:"
echo "  - Consensus HLA types: ${OUTDIR}/<sample>/<sample>_consensus.txt"
echo "  - Visualizations: ${OUTDIR}/<sample>/visualizations/"
echo "  - Summary report: ${OUTDIR}/summary/hla_summary_report.html"
echo "  - MultiQC report: ${OUTDIR}/multiqc/multiqc_report.html"
echo ""
