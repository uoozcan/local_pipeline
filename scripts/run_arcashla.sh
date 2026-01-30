#!/bin/bash
#
# arcasHLA wrapper script for BAM files
# arcasHLA performs HLA typing from RNA-seq data
#
# Usage: ./run_arcashla.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-j THREADS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER="${PROJECT_DIR}/hla_references/containers/arcashla.sif"

# Default values
THREADS=8
GENES="A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1"

usage() {
    echo "arcasHLA wrapper for BAM files (optimized for RNA-seq)"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name"
    echo "  -b    Input BAM file (RNA-seq aligned to reference)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -j    Number of threads (default: 8)"
    echo "  -g    HLA genes to type (comma-separated, default: A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1)"
    echo "  -h    Show this help message"
    echo ""
    echo "Note: arcasHLA is optimized for RNA-seq data"
    echo ""
    exit 1
}

while getopts "n:b:o:j:g:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        g) GENES="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$BAM_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
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

# Get absolute paths
BAM_FILE=$(realpath "$BAM_FILE")
BAM_DIR=$(dirname "$BAM_FILE")
BAM_NAME=$(basename "$BAM_FILE")

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "arcasHLA HLA Typing Pipeline (RNA-seq)"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Threads:    $THREADS"
echo "Genes:      $GENES"
echo "=============================================="
echo ""

# Step 1: Extract chromosome 6 reads
echo "[Step 1/2] Extracting chromosome 6 reads..."

singularity exec \
    --bind "${BAM_DIR}:/bam_dir" \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    /home/arcasHLA-master/arcasHLA extract \
    "/bam_dir/${BAM_NAME}" \
    -o /output \
    -t "$THREADS" \
    -v

# Step 2: Genotype HLA
echo ""
echo "[Step 2/2] Running HLA genotyping..."

# Find extracted FASTQ files
EXTRACTED_R1=$(find "${WORK_DIR}" -name "*.extracted.1.fq.gz" 2>/dev/null | head -1)
EXTRACTED_R2=$(find "${WORK_DIR}" -name "*.extracted.2.fq.gz" 2>/dev/null | head -1)

if [ -z "$EXTRACTED_R1" ]; then
    echo "Error: No extracted reads found. Check if BAM contains chromosome 6 reads."
    exit 1
fi

R1_NAME=$(basename "$EXTRACTED_R1")
R2_NAME=$(basename "$EXTRACTED_R2")

singularity exec \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    /home/arcasHLA-master/arcasHLA genotype \
    "/output/${R1_NAME}" \
    "/output/${R2_NAME}" \
    -g "$GENES" \
    -o /output \
    -t "$THREADS" \
    -v

echo ""
echo "=============================================="
echo "arcasHLA Complete!"
echo "=============================================="
echo "Results directory: ${WORK_DIR}/"

# Find and display results
RESULT_FILE=$(find "${WORK_DIR}" -name "*.genotype.json" 2>/dev/null | head -1)
if [ -n "$RESULT_FILE" ] && [ -f "$RESULT_FILE" ]; then
    echo ""
    echo "Genotype results:"
    cat "$RESULT_FILE" | python3 -m json.tool 2>/dev/null || cat "$RESULT_FILE"
fi

# Also check for TSV results
TSV_FILE=$(find "${WORK_DIR}" -name "*genes.tsv" 2>/dev/null | head -1)
if [ -n "$TSV_FILE" ] && [ -f "$TSV_FILE" ]; then
    echo ""
    echo "Gene-level results:"
    cat "$TSV_FILE"
fi
echo "=============================================="
