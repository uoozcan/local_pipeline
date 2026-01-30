#!/bin/bash
#
# HLA-HD wrapper script for BAM files
# HLA-HD performs HLA typing from WGS/WES/RNA-seq data
#
# Usage: ./run_hlahd.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-r REFERENCE] [-j THREADS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER="${PROJECT_DIR}/hla_references/containers/hlahd.sif"
LOCAL_BIN="${PROJECT_DIR}/spechla_local/local_bin"

# Default values
THREADS=8
REFERENCE="hg38"

usage() {
    echo "HLA-HD wrapper for BAM files"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name"
    echo "  -b    Input BAM file (sorted and indexed)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -r    Reference genome: hg38 or hg19 (default: hg38)"
    echo "  -j    Number of threads (default: 8)"
    echo "  -h    Show this help message"
    echo ""
    exit 1
}

while getopts "n:b:o:r:j:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
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

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"

echo "=============================================="
echo "HLA-HD HLA Typing Pipeline"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $OUTPUT_DIR"
echo "Threads:    $THREADS"
echo "=============================================="
echo ""

# Set PATH for local tools
export PATH="${LOCAL_BIN}:/home/umut/miniconda3/bin:$PATH"

# Step 1: Extract HLA reads and convert to FASTQ
echo "[Step 1/2] Extracting HLA reads from BAM..."

# Define HLA region
if [ "$REFERENCE" == "hg38" ]; then
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28510120-33480577"
    else
        HLA_REGION="6:28510120-33480577"
    fi
else
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28477797-33448354"
    else
        HLA_REGION="6:28477797-33448354"
    fi
fi

# Extract and convert to FASTQ
samtools view -@ "$THREADS" -b "$BAM_FILE" "$HLA_REGION" | \
    samtools sort -n -@ "$THREADS" -o "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam"

samtools fastq -@ "$THREADS" \
    -1 "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" \
    -2 "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz" \
    -0 /dev/null -s /dev/null \
    "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam" 2>&1 | grep -E "processed|discarded" || true

# Decompress for HLA-HD (requires uncompressed FASTQ)
gunzip -c "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" > "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq"
gunzip -c "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz" > "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq"

echo "  Extracted HLA reads to FASTQ files"

# Step 2: Run HLA-HD
echo ""
echo "[Step 2/2] Running HLA-HD..."

singularity exec \
    --bind "${WORK_DIR}:/data" \
    "$CONTAINER" \
    /app/hlahd.1.4.0/bin/hlahd.sh \
    -t "$THREADS" \
    -f /app/hlahd.1.4.0/freq_data \
    "/data/${SAMPLE_NAME}_R1.fastq" \
    "/data/${SAMPLE_NAME}_R2.fastq" \
    /app/hlahd.1.4.0/HLA_gene.split.txt \
    /app/hlahd.1.4.0/dictionary \
    "$SAMPLE_NAME" \
    /data

# Cleanup
rm -f "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam"
rm -f "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq" "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq"

echo ""
echo "=============================================="
echo "HLA-HD Complete!"
echo "=============================================="
echo "Results: ${WORK_DIR}/result/"
if [ -f "${WORK_DIR}/result/${SAMPLE_NAME}_final.result.txt" ]; then
    echo ""
    cat "${WORK_DIR}/result/${SAMPLE_NAME}_final.result.txt"
fi
echo "=============================================="
