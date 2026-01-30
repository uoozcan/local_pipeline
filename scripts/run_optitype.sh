#!/bin/bash
#
# OptiType wrapper script for BAM files
# OptiType performs HLA Class I typing (HLA-A, -B, -C) from WGS/WES/RNA-seq data
#
# Usage: ./run_optitype.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-t TYPE] [-r REFERENCE] [-j THREADS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER="${PROJECT_DIR}/hla_references/containers/optitype.sif"
LOCAL_BIN="${PROJECT_DIR}/spechla_local/local_bin"

# Default values
THREADS=8
REFERENCE="hg38"
SEQ_TYPE="dna"  # dna or rna

usage() {
    echo "OptiType wrapper for BAM files (HLA Class I only: A, B, C)"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name"
    echo "  -b    Input BAM file (sorted and indexed)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -t    Sequencing type: dna or rna (default: dna)"
    echo "  -r    Reference genome: hg38 or hg19 (default: hg38)"
    echo "  -j    Number of threads (default: 8)"
    echo "  -h    Show this help message"
    echo ""
    echo "Note: OptiType only types HLA Class I genes (A, B, C)"
    echo ""
    exit 1
}

while getopts "n:b:o:t:r:j:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) SEQ_TYPE="$OPTARG" ;;
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

# Validate seq type
if [ "$SEQ_TYPE" != "dna" ] && [ "$SEQ_TYPE" != "rna" ]; then
    echo "Error: Sequencing type must be 'dna' or 'rna'"
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

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "OptiType HLA Typing Pipeline (Class I only)"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Seq Type:   $SEQ_TYPE"
echo "=============================================="
echo ""

# Set PATH for local tools
export PATH="${LOCAL_BIN}:/home/umut/miniconda3/bin:$PATH"

# Step 1: Extract HLA reads and convert to FASTQ
echo "[Step 1/2] Extracting HLA reads from BAM..."

# Define HLA region (Class I is in 6p21.3)
if [ "$REFERENCE" == "hg38" ]; then
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:29941260-31357179"  # Class I region
    else
        HLA_REGION="6:29941260-31357179"
    fi
else
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:29909037-31324989"
    else
        HLA_REGION="6:29909037-31324989"
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

echo "  Extracted HLA reads to FASTQ files"

# Step 2: Run OptiType
echo ""
echo "[Step 2/2] Running OptiType..."

# Determine the flag for seq type
if [ "$SEQ_TYPE" == "dna" ]; then
    TYPE_FLAG="--dna"
else
    TYPE_FLAG="--rna"
fi

singularity exec \
    --bind "${WORK_DIR}:/data" \
    "$CONTAINER" \
    python /usr/local/bin/OptiType/OptiTypePipeline.py \
    --input "/data/${SAMPLE_NAME}_R1.fastq.gz" "/data/${SAMPLE_NAME}_R2.fastq.gz" \
    $TYPE_FLAG \
    --outdir /data \
    --prefix "$SAMPLE_NAME"

# Cleanup
rm -f "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam"

echo ""
echo "=============================================="
echo "OptiType Complete!"
echo "=============================================="
echo "Results directory: ${WORK_DIR}/"

# Find and display results
RESULT_FILE=$(find "${WORK_DIR}" -name "${SAMPLE_NAME}*_result.tsv" 2>/dev/null | head -1)
if [ -n "$RESULT_FILE" ] && [ -f "$RESULT_FILE" ]; then
    echo ""
    cat "$RESULT_FILE"
fi
echo "=============================================="
