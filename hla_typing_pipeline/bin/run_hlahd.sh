#!/bin/bash
#
# HLA-HD Wrapper Script
# Performs HLA typing from BAM files using HLA-HD
#
# Usage: ./run_hlahd.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-d DATABASE] [-j THREADS]
#

set -euo pipefail

# Default values
THREADS=8
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect paths (can be overridden)
CONTAINER="${CONTAINER:-${SCRIPT_DIR}/../../hla_references/containers/hlahd.sif}"
HLAHD_DB="${HLAHD_DB:-${SCRIPT_DIR}/../../hla_references/databases/hlahd_db}"

usage() {
    cat << EOF
HLA-HD HLA Typing Wrapper

Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]

Required arguments:
  -n    Sample name (will be used for output files)
  -b    Input BAM file (sorted and indexed)
  -o    Output directory

Optional arguments:
  -d    HLA-HD database directory (default: auto-detect)
  -c    Singularity container path (default: auto-detect)
  -j    Number of threads (default: 8)
  -h    Show this help message

Environment variables:
  CONTAINER    Path to HLA-HD Singularity container
  HLAHD_DB     Path to HLA-HD database directory

Example:
  $0 -n Sample1 -b Sample1.bam -o ./results -j 4

EOF
    exit 1
}

while getopts "n:b:o:d:c:j:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) HLAHD_DB="$OPTARG" ;;
        c) CONTAINER="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Validate required arguments
if [ -z "${SAMPLE_NAME:-}" ] || [ -z "${BAM_FILE:-}" ] || [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate files
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    echo "Error: Container not found: $CONTAINER"
    exit 1
fi

if [ ! -d "$HLAHD_DB" ]; then
    echo "Error: HLA-HD database not found: $HLAHD_DB"
    exit 1
fi

# Get absolute paths
BAM_FILE=$(realpath "$BAM_FILE")
BAM_DIR=$(dirname "$BAM_FILE")
BAM_NAME=$(basename "$BAM_FILE")
HLAHD_DB=$(realpath "$HLAHD_DB")

# Check for BAM index
if [ ! -f "${BAM_FILE}.bai" ] && [ ! -f "${BAM_FILE%.bam}.bai" ]; then
    echo "Warning: BAM index not found. Creating index..."
    singularity exec "$CONTAINER" samtools index "$BAM_FILE" 2>/dev/null || \
    samtools index "$BAM_FILE"
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "HLA-HD HLA Typing Pipeline"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Database:   $HLAHD_DB"
echo "Threads:    $THREADS"
echo "=============================================="
echo ""

# Create temporary directory for intermediate files
TEMP_DIR="${WORK_DIR}/temp"
mkdir -p "$TEMP_DIR"

# Step 1: Extract HLA reads
echo "[Step 1/3] Extracting HLA reads from BAM..."
singularity exec \
    --bind "${BAM_DIR}:/bam_dir" \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    bash -c "
        cd /output/temp
        # Extract MHC region (chr6:28-34Mb)
        samtools view -b -h /bam_dir/${BAM_NAME} chr6:28000000-34000000 > hla_region.bam 2>/dev/null || \
        samtools view -b -h /bam_dir/${BAM_NAME} 6:28000000-34000000 > hla_region.bam

        # Extract unmapped reads
        samtools view -b -f 4 /bam_dir/${BAM_NAME} > unmapped.bam

        # Merge and sort by name
        samtools merge -f merged.bam hla_region.bam unmapped.bam
        samtools sort -n -@ ${THREADS} merged.bam -o sorted.bam

        # Convert to FASTQ
        samtools fastq -@ ${THREADS} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam
    "
echo "  Extracted HLA reads to FASTQ files"

# Step 2: Run HLA-HD
echo ""
echo "[Step 2/3] Running HLA-HD..."
singularity exec \
    --bind "${WORK_DIR}:/output" \
    --bind "${HLAHD_DB}:/database" \
    "$CONTAINER" \
    bash -c "
        cd /output
        hlahd.sh -t ${THREADS} -m 100 -c 0.95 -f /database/freq_data \
            /output/temp/R1.fastq /output/temp/R2.fastq \
            /database/HLA_gene.split /database/dictionary \
            ${SAMPLE_NAME} .
    "

# Step 3: Parse and format results
echo ""
echo "[Step 3/3] Formatting results..."

RESULT_FILE="${WORK_DIR}/${SAMPLE_NAME}/result/${SAMPLE_NAME}_final.result.txt"
OUTPUT_FILE="${WORK_DIR}/${SAMPLE_NAME}_hlahd.txt"

if [ -f "$RESULT_FILE" ]; then
    cp "$RESULT_FILE" "$OUTPUT_FILE"
    echo "  Results saved to: $OUTPUT_FILE"
else
    echo "# HLA-HD results for ${SAMPLE_NAME}" > "$OUTPUT_FILE"
    echo "# Error: No results generated" >> "$OUTPUT_FILE"
    echo "  Warning: No results file found"
fi

# Clean up temp files
rm -rf "$TEMP_DIR"

echo ""
echo "=============================================="
echo "HLA-HD Complete!"
echo "=============================================="
echo "Results: $OUTPUT_FILE"
echo ""

if [ -f "$OUTPUT_FILE" ]; then
    echo "HLA Typing Results:"
    echo "-------------------"
    cat "$OUTPUT_FILE"
fi
echo "=============================================="
