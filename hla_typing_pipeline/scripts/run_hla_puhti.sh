#!/bin/bash
#
# HLA Typing Pipeline - Puhti Run Script
# Convenient wrapper for running HLA typing on CSC Puhti
#
# Usage:
#   ./run_hla_puhti.sh --bam sample.bam
#   ./run_hla_puhti.sh --fastq sample_R1.fq.gz sample_R2.fq.gz
#   ./run_hla_puhti.sh --samplesheet samples.csv
#   ./run_hla_puhti.sh --bam-dir /path/to/bams
#

set -euo pipefail

# ============================================
# Configuration
# ============================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Auto-detect project ID from path or environment
if [[ "$PWD" =~ /scratch/(project_[0-9]+)/ ]]; then
    PROJECT_ID="${BASH_REMATCH[1]}"
elif [[ -n "${CSC_PROJECT:-}" ]]; then
    PROJECT_ID="$CSC_PROJECT"
else
    echo "ERROR: Could not detect CSC project ID."
    echo "Either:"
    echo "  1. Run from /scratch/project_XXXXXXX/ directory"
    echo "  2. Set CSC_PROJECT environment variable: export CSC_PROJECT=project_XXXXXXX"
    exit 1
fi

BASE_DIR="/scratch/${PROJECT_ID}/${USER}/hla_analysis"

# Default SLURM settings
PARTITION="small"
TIME="12:00:00"
MEM="180G"
CPUS="40"

# Default pipeline settings
TOOLS="spechla,hlahd"
REFERENCE="hg38"
OUTDIR="${BASE_DIR}/results"

# ============================================
# Help Message
# ============================================

show_help() {
    cat << EOF
HLA Typing Pipeline - Puhti Runner
==================================

Usage: $(basename "$0") [INPUT_OPTIONS] [PIPELINE_OPTIONS] [SLURM_OPTIONS]

INPUT OPTIONS (choose one):
  --bam <file>                Single BAM file
  --fastq <R1> <R2>          Paired FASTQ files
  --samplesheet <csv>         Samplesheet with multiple samples
  --bam-dir <dir>             Process all BAM files in directory

PIPELINE OPTIONS:
  --tools <list>              HLA tools (default: spechla,hlahd)
                              Options: spechla,hlahd,hlala,arcashla,optitype
  --reference <ref>           Reference: hg38 or hg19 (default: hg38)
  --resolution <res>          Resolution: 2-field or 4-field (default: 2-field)
  --seq-type <type>           Sequence type for OptiType: dna or rna (default: dna)
  --outdir <dir>              Output directory (default: ${BASE_DIR}/results)

SLURM OPTIONS:
  --partition <name>          SLURM partition (default: small)
  --time <HH:MM:SS>           Time limit (default: 12:00:00)
  --mem <size>                Memory (default: 180G)
  --cpus <n>                  CPUs (default: 40)
  --interactive               Run interactively instead of submitting

OTHER OPTIONS:
  --dry-run                   Show command without executing
  --help                      Show this help message

EXAMPLES:
  # Single BAM file
  $(basename "$0") --bam sample.bam

  # Paired FASTQ files
  $(basename "$0") --fastq sample_R1.fastq.gz sample_R2.fastq.gz

  # Multiple samples from samplesheet
  $(basename "$0") --samplesheet samples.csv

  # All BAM files in a directory
  $(basename "$0") --bam-dir /path/to/bam_files/

  # Interactive mode for testing
  $(basename "$0") --bam sample.bam --interactive --time 04:00:00

  # Use multiple tools
  $(basename "$0") --bam sample.bam --tools spechla,hlahd,arcashla,optitype

EOF
    exit 0
}

# ============================================
# Parse Arguments
# ============================================

INPUT_TYPE=""
INPUT_FILE=""
FASTQ_R1=""
FASTQ_R2=""
SAMPLESHEET=""
BAM_DIR=""
INTERACTIVE=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam)
            INPUT_TYPE="bam"
            INPUT_FILE="$2"
            shift 2
            ;;
        --fastq)
            INPUT_TYPE="fastq"
            FASTQ_R1="$2"
            FASTQ_R2="$3"
            shift 3
            ;;
        --samplesheet)
            INPUT_TYPE="samplesheet"
            SAMPLESHEET="$2"
            shift 2
            ;;
        --bam-dir)
            INPUT_TYPE="bam-dir"
            BAM_DIR="$2"
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
        --seq-type)
            SEQ_TYPE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --interactive)
            INTERACTIVE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ============================================
# Validate Input
# ============================================

if [[ -z "$INPUT_TYPE" ]]; then
    echo "ERROR: No input specified"
    echo "Use --help for usage information"
    exit 1
fi

# Make paths absolute
make_absolute() {
    local path="$1"
    if [[ ! "$path" = /* ]]; then
        echo "$PWD/$path"
    else
        echo "$path"
    fi
}

case $INPUT_TYPE in
    bam)
        INPUT_FILE=$(make_absolute "$INPUT_FILE")
        if [[ ! -f "$INPUT_FILE" ]]; then
            echo "ERROR: BAM file not found: $INPUT_FILE"
            exit 1
        fi
        ;;
    fastq)
        FASTQ_R1=$(make_absolute "$FASTQ_R1")
        FASTQ_R2=$(make_absolute "$FASTQ_R2")
        if [[ ! -f "$FASTQ_R1" ]]; then
            echo "ERROR: R1 FASTQ not found: $FASTQ_R1"
            exit 1
        fi
        if [[ ! -f "$FASTQ_R2" ]]; then
            echo "ERROR: R2 FASTQ not found: $FASTQ_R2"
            exit 1
        fi
        ;;
    samplesheet)
        SAMPLESHEET=$(make_absolute "$SAMPLESHEET")
        if [[ ! -f "$SAMPLESHEET" ]]; then
            echo "ERROR: Samplesheet not found: $SAMPLESHEET"
            exit 1
        fi
        ;;
    bam-dir)
        BAM_DIR=$(make_absolute "$BAM_DIR")
        if [[ ! -d "$BAM_DIR" ]]; then
            echo "ERROR: BAM directory not found: $BAM_DIR"
            exit 1
        fi
        ;;
esac

# ============================================
# Create Directories
# ============================================

mkdir -p "${BASE_DIR}/logs"
mkdir -p "$OUTDIR"

# ============================================
# Build Command
# ============================================

echo "=============================================="
echo "HLA Typing Pipeline - Puhti"
echo "=============================================="
echo "Project: $PROJECT_ID"
echo "Input type: $INPUT_TYPE"
echo "Tools: $TOOLS"
echo "Output: $OUTDIR"
echo ""

case $INPUT_TYPE in
    bam)
        echo "Input: $INPUT_FILE"
        SUBMIT_SCRIPT="${PIPELINE_DIR}/scripts/submit_hla_single.sh"
        SUBMIT_ARGS="$INPUT_FILE --tools $TOOLS --reference $REFERENCE --outdir $OUTDIR"
        ;;
    fastq)
        echo "Input R1: $FASTQ_R1"
        echo "Input R2: $FASTQ_R2"
        SUBMIT_SCRIPT="${PIPELINE_DIR}/scripts/submit_hla_single.sh"
        SUBMIT_ARGS="$FASTQ_R1 --fastq2 $FASTQ_R2 --tools $TOOLS --reference $REFERENCE --outdir $OUTDIR"
        ;;
    samplesheet)
        echo "Samplesheet: $SAMPLESHEET"
        SAMPLE_COUNT=$(tail -n +2 "$SAMPLESHEET" | wc -l)
        echo "Samples: $SAMPLE_COUNT"
        SUBMIT_SCRIPT="${PIPELINE_DIR}/scripts/submit_hla_batch.sh"
        SUBMIT_ARGS="$SAMPLESHEET --tools $TOOLS --reference $REFERENCE --outdir $OUTDIR"
        ;;
    bam-dir)
        echo "BAM directory: $BAM_DIR"

        # Generate samplesheet from BAM directory
        GENERATED_SAMPLESHEET="${BASE_DIR}/generated_samplesheet.csv"
        echo "sample_id,bam_path" > "$GENERATED_SAMPLESHEET"

        for bam in "$BAM_DIR"/*.bam; do
            if [[ -f "$bam" ]]; then
                sample_id=$(basename "$bam" .bam)
                echo "${sample_id},${bam}" >> "$GENERATED_SAMPLESHEET"
            fi
        done

        SAMPLE_COUNT=$(tail -n +2 "$GENERATED_SAMPLESHEET" | wc -l)
        echo "Found $SAMPLE_COUNT BAM files"
        echo "Generated samplesheet: $GENERATED_SAMPLESHEET"

        SUBMIT_SCRIPT="${PIPELINE_DIR}/scripts/submit_hla_batch.sh"
        SUBMIT_ARGS="$GENERATED_SAMPLESHEET --tools $TOOLS --reference $REFERENCE --outdir $OUTDIR"
        ;;
esac

echo ""

# ============================================
# Submit or Run
# ============================================

if [[ "$DRY_RUN" == true ]]; then
    echo "DRY RUN - Would execute:"
    echo ""
    if [[ "$INTERACTIVE" == true ]]; then
        echo "sinteractive --account $PROJECT_ID --partition $PARTITION --time $TIME --mem $MEM --cpus-per-task $CPUS"
        echo "bash $SUBMIT_SCRIPT $SUBMIT_ARGS"
    else
        echo "sbatch --account $PROJECT_ID --partition $PARTITION --time $TIME --mem $MEM --cpus-per-task $CPUS $SUBMIT_SCRIPT $SUBMIT_ARGS"
    fi
    exit 0
fi

if [[ "$INTERACTIVE" == true ]]; then
    echo "Starting interactive session..."
    echo "After session starts, run:"
    echo "  bash $SUBMIT_SCRIPT $SUBMIT_ARGS"
    echo ""

    sinteractive \
        --account "$PROJECT_ID" \
        --partition "$PARTITION" \
        --time "$TIME" \
        --mem "$MEM" \
        --cpus-per-task "$CPUS"
else
    echo "Submitting SLURM job..."

    JOB_ID=$(sbatch \
        --account "$PROJECT_ID" \
        --partition "$PARTITION" \
        --time "$TIME" \
        --mem "$MEM" \
        --cpus-per-task "$CPUS" \
        --parsable \
        "$SUBMIT_SCRIPT" $SUBMIT_ARGS)

    echo ""
    echo "=============================================="
    echo "Job submitted successfully!"
    echo "=============================================="
    echo "Job ID: $JOB_ID"
    echo ""
    echo "Monitor with:"
    echo "  squeue -j $JOB_ID"
    echo "  tail -f ${BASE_DIR}/logs/hla_${JOB_ID}.out"
    echo ""
    echo "Results will be in: $OUTDIR"
fi
