#!/bin/bash
#
# HLA*LA wrapper script for BAM files
# HLA*LA performs HLA typing from WGS/WES data using graph-based alignment
#
# Usage: ./run_hlala.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-j THREADS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER="${PROJECT_DIR}/hla_references/containers/hlala.sif"
HLALA_GRAPHS="${PROJECT_DIR}/hla_references/databases/hlala_graphs"

# Default values
THREADS=8
GRAPH="PRG_MHC_GRCh38_withIMGT"

usage() {
    echo "HLA*LA wrapper for BAM files (graph-based WGS/WES typing)"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name (alphanumeric characters only)"
    echo "  -b    Input BAM file (sorted and indexed)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -j    Number of threads (default: 8)"
    echo "  -g    Graph to use (default: PRG_MHC_GRCh38_withIMGT)"
    echo "  -h    Show this help message"
    echo ""
    exit 1
}

while getopts "n:b:o:j:g:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        g) GRAPH="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$BAM_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate sample name (alphanumeric only)
if [[ ! "$SAMPLE_NAME" =~ ^[a-zA-Z0-9_]+$ ]]; then
    echo "Error: Sample name must contain only alphanumeric characters and underscores"
    exit 1
fi

# Check files exist
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    echo "Error: Container not found: $CONTAINER"
    exit 1
fi

# Check graph database exists
if [ ! -d "${HLALA_GRAPHS}/${GRAPH}" ]; then
    echo "Error: HLA*LA graph not found: ${HLALA_GRAPHS}/${GRAPH}"
    echo "Please run: ./scripts/setup_hla_databases.sh"
    exit 1
fi

# Get absolute paths
BAM_FILE=$(realpath "$BAM_FILE")
BAM_DIR=$(dirname "$BAM_FILE")
BAM_NAME=$(basename "$BAM_FILE")

# Check for BAM index
if [ ! -f "${BAM_FILE}.bai" ] && [ ! -f "${BAM_FILE%.*}.bai" ]; then
    echo "Error: BAM index not found. Please index your BAM file."
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "HLA*LA HLA Typing Pipeline"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Threads:    $THREADS"
echo "Graph:      $GRAPH"
echo "=============================================="
echo ""

# Run HLA*LA
echo "[Running HLA*LA...]"

singularity exec \
    --bind "${BAM_DIR}:/bam_dir" \
    --bind "${WORK_DIR}:/output" \
    --bind "${HLALA_GRAPHS}:/usr/local/bin/HLA-LA/graphs" \
    "$CONTAINER" \
    perl /usr/local/bin/HLA-LA/src/HLA-LA.pl \
    --BAM "/bam_dir/${BAM_NAME}" \
    --graph "$GRAPH" \
    --sampleID "$SAMPLE_NAME" \
    --maxThreads "$THREADS" \
    --workingDir /output

echo ""
echo "=============================================="
echo "HLA*LA Complete!"
echo "=============================================="
echo "Results directory: ${WORK_DIR}/"

# Find and display results
RESULT_FILE="${WORK_DIR}/${SAMPLE_NAME}/hla/R1_bestguess_G.txt"
if [ -f "$RESULT_FILE" ]; then
    echo ""
    echo "Best guess results (G group):"
    cat "$RESULT_FILE"
fi

RESULT_FILE2="${WORK_DIR}/${SAMPLE_NAME}/hla/R1_bestguess.txt"
if [ -f "$RESULT_FILE2" ]; then
    echo ""
    echo "Best guess results:"
    cat "$RESULT_FILE2"
fi
echo "=============================================="
