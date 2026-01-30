#!/bin/bash
#SBATCH --job-name=hla_batch
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/hla_batch_%j.out
#SBATCH --error=logs/hla_batch_%j.err
# NOTE: Set account with: sbatch --account=YOUR_PROJECT_ID submit_hla_batch.sh

#
# HLA Typing Pipeline - Batch SLURM Submission Script
# Processes multiple samples from a samplesheet
#
# Usage: sbatch submit_hla_batch.sh <samplesheet.csv> [options]
#
# Samplesheet format (BAM):
#   sample_id,bam_path
#   Sample1,/path/to/Sample1.bam
#   Sample2,/path/to/Sample2.bam
#
# Samplesheet format (FASTQ):
#   sample_id,fastq_1,fastq_2
#   Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
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

# Default parameters
TOOLS="spechla,hlahd"
REFERENCE="hg38"
RESOLUTION="2-field"
SEQ_TYPE="dna"
OUTDIR="${BASE_DIR}/results"

# ============================================
# Parse Arguments
# ============================================

SAMPLESHEET=""
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
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
            if [[ -z "$SAMPLESHEET" ]]; then
                SAMPLESHEET="$1"
            fi
            shift
            ;;
    esac
done

# Validate samplesheet
if [[ -z "$SAMPLESHEET" ]]; then
    echo "ERROR: No samplesheet specified"
    echo ""
    echo "Usage: sbatch submit_hla_batch.sh <samplesheet.csv> [options]"
    echo ""
    echo "Options:"
    echo "  --tools <list>       HLA tools to use (default: spechla,hlahd)"
    echo "  --reference <ref>    Reference genome: hg38 or hg19 (default: hg38)"
    echo "  --resolution <res>   Output resolution (default: 2-field)"
    echo "  --seq_type <type>    Sequence type for OptiType (default: dna)"
    echo "  --outdir <dir>       Output directory"
    exit 1
fi

# Make samplesheet path absolute
if [[ ! "$SAMPLESHEET" = /* ]]; then
    SAMPLESHEET="${SLURM_SUBMIT_DIR:-$PWD}/$SAMPLESHEET"
fi

if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "ERROR: Samplesheet not found: $SAMPLESHEET"
    exit 1
fi

# ============================================
# Setup Environment
# ============================================

echo "=============================================="
echo "HLA Typing Pipeline - Batch Mode"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

# Load modules
module purge
module load nextflow/23.10.0 2>/dev/null || module load nextflow
module load singularity 2>/dev/null || true
module load samtools 2>/dev/null || module load biokit 2>/dev/null || true

# Set Singularity cache
export SINGULARITY_CACHEDIR="${BASE_DIR}/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
mkdir -p "$SINGULARITY_CACHEDIR"

# Create directories
mkdir -p "$OUTDIR"
mkdir -p "${BASE_DIR}/logs"

# Count samples
SAMPLE_COUNT=$(tail -n +2 "$SAMPLESHEET" | wc -l)

echo "Configuration:"
echo "  Samplesheet: $SAMPLESHEET"
echo "  Sample count: $SAMPLE_COUNT"
echo "  Output: $OUTDIR"
echo "  Tools: $TOOLS"
echo "  Reference: $REFERENCE"
echo ""

# ============================================
# Validate and Index BAM Files
# ============================================

echo "Validating input files..."

# Check header to determine input type
HEADER=$(head -1 "$SAMPLESHEET")

if [[ "$HEADER" == *"bam"* ]]; then
    echo "Input type: BAM files"

    # Check for BAM indices
    MISSING_INDEX=0
    while IFS=',' read -r sample_id bam_path; do
        [[ "$sample_id" == "sample_id" ]] && continue  # Skip header

        bam_path=$(echo "$bam_path" | tr -d '\r')  # Remove carriage return if present

        if [[ ! -f "$bam_path" ]]; then
            echo "  WARNING: BAM not found: $bam_path"
        elif [[ ! -f "${bam_path}.bai" ]] && [[ ! -f "${bam_path%.bam}.bai" ]]; then
            echo "  Creating index for: $bam_path"
            samtools index "$bam_path" &
            MISSING_INDEX=$((MISSING_INDEX + 1))
        fi
    done < "$SAMPLESHEET"

    # Wait for indexing to complete
    if [[ $MISSING_INDEX -gt 0 ]]; then
        echo "  Waiting for ${MISSING_INDEX} BAM indices to be created..."
        wait
        echo "  Indexing complete"
    fi
else
    echo "Input type: FASTQ files"
fi

echo ""

# ============================================
# Run Nextflow Pipeline
# ============================================

echo "Starting HLA typing pipeline..."
echo "=============================================="

cd "$BASE_DIR"

nextflow run "${PIPELINE_DIR}/main.nf" \
    --input_samplesheet "$SAMPLESHEET" \
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
echo "Samples processed: $SAMPLE_COUNT"
echo ""
echo "Results available at: $OUTDIR"
echo ""
echo "Key output files:"
echo "  - Per-sample results: ${OUTDIR}/<sample>/"
echo "  - Summary report: ${OUTDIR}/summary/hla_summary_report.html"
echo "  - Summary statistics: ${OUTDIR}/summary/hla_summary_statistics.tsv"
echo "  - MultiQC report: ${OUTDIR}/multiqc/multiqc_report.html"
echo ""
echo "Quick summary:"
if [[ -f "${OUTDIR}/summary/hla_summary_statistics.tsv" ]]; then
    echo ""
    head -20 "${OUTDIR}/summary/hla_summary_statistics.tsv"
fi
