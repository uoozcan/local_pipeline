#!/bin/bash
#
# SpecHLA wrapper script for whole-genome/exome BAM files
# Automates: HLA read extraction -> FASTQ conversion -> SpecHLA typing
#
# Usage: ./run_spechla_bam.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-r REFERENCE] [-j THREADS]
#
# Example:
#   ./run_spechla_bam.sh -n FH_7087 -b /path/to/sample.bam -o ./results -r hg38 -j 8
#

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set PATH to include local tools and miniconda
export PATH="${SCRIPT_DIR}/local_bin:/home/umut/miniconda3/bin:$PATH"

# Default values
THREADS=8
REFERENCE="hg38"
KEEP_TEMP=0

# Usage function
usage() {
    echo "SpecHLA wrapper for BAM files"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name (e.g., FH_7087)"
    echo "  -b    Input BAM file (sorted and indexed, aligned to hg38 or hg19)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -r    Reference genome: hg38 or hg19 (default: hg38)"
    echo "  -j    Number of threads (default: 8)"
    echo "  -k    Keep temporary files (default: remove)"
    echo "  -h    Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 -n FH_7087 -b /path/to/sample.bam -o ./results -r hg38 -j 8"
    echo ""
    exit 1
}

# Parse arguments
while getopts "n:b:o:r:j:kh" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        k) KEEP_TEMP=1 ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$BAM_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

# Check BAM index exists
if [ ! -f "${BAM_FILE}.bai" ] && [ ! -f "${BAM_FILE%.*}.bai" ]; then
    echo "Error: BAM index not found. Please index your BAM file with: samtools index $BAM_FILE"
    exit 1
fi

# Validate reference
if [ "$REFERENCE" != "hg38" ] && [ "$REFERENCE" != "hg19" ]; then
    echo "Error: Reference must be 'hg38' or 'hg19'"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"

echo "=============================================="
echo "SpecHLA HLA Typing Pipeline"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $OUTPUT_DIR"
echo "Reference:  $REFERENCE"
echo "Threads:    $THREADS"
echo "=============================================="
echo ""

# Step 1: Extract HLA reads from BAM
echo "[Step 1/3] Extracting HLA reads from BAM file..."

# Define HLA region based on reference genome
if [ "$REFERENCE" == "hg38" ]; then
    # Check chromosome naming convention
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28510120-33480577"
    else
        HLA_REGION="6:28510120-33480577"
    fi
else
    # hg19
    CHR_PREFIX=$(samtools view -H "$BAM_FILE" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28477797-33448354"
    else
        HLA_REGION="6:28477797-33448354"
    fi
fi

echo "  HLA region: $HLA_REGION"

# Extract reads from HLA region
samtools view -@ "$THREADS" -b "$BAM_FILE" "$HLA_REGION" > "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam" 2>/dev/null

# Index extracted BAM
samtools index "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam"

# Count extracted reads
READ_COUNT=$(samtools view -c "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam")
echo "  Extracted $READ_COUNT reads from HLA region"

if [ "$READ_COUNT" -lt 1000 ]; then
    echo "Warning: Low number of reads extracted. HLA typing results may be unreliable."
fi

# Step 2: Convert BAM to FASTQ
echo ""
echo "[Step 2/3] Converting BAM to FASTQ..."

# Sort by read name for proper paired-end extraction
samtools sort -n -@ "$THREADS" "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam" -o "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam"

# Convert to FASTQ
samtools fastq -@ "$THREADS" \
    -1 "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" \
    -2 "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz" \
    -0 /dev/null \
    -s /dev/null \
    "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam" 2>&1 | grep -E "processed|discarded" || true

# Check FASTQ files
R1_SIZE=$(stat -c%s "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" 2>/dev/null || echo "0")
R2_SIZE=$(stat -c%s "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz" 2>/dev/null || echo "0")

if [ "$R1_SIZE" -lt 1000 ] || [ "$R2_SIZE" -lt 1000 ]; then
    echo "Error: FASTQ files are too small. Check your BAM file."
    exit 1
fi

echo "  Created FASTQ files: ${SAMPLE_NAME}_R1.fastq.gz, ${SAMPLE_NAME}_R2.fastq.gz"

# Step 3: Run SpecHLA
echo ""
echo "[Step 3/3] Running SpecHLA HLA typing..."

bash "${SCRIPT_DIR}/script/whole/SpecHLA.sh" \
    -n "$SAMPLE_NAME" \
    -1 "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" \
    -2 "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz" \
    -o "$OUTPUT_DIR" \
    -j "$THREADS" \
    -u 1

# Cleanup temporary files
if [ "$KEEP_TEMP" -eq 0 ]; then
    echo ""
    echo "Cleaning up temporary files..."
    rm -f "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam"
    rm -f "${WORK_DIR}/${SAMPLE_NAME}.hla_extract.bam.bai"
    rm -f "${WORK_DIR}/${SAMPLE_NAME}.namesort.bam"
fi

# Display results
echo ""
echo "=============================================="
echo "HLA Typing Complete!"
echo "=============================================="
echo ""
echo "Results file: ${WORK_DIR}/hla.result.txt"
echo ""
if [ -f "${WORK_DIR}/hla.result.txt" ]; then
    cat "${WORK_DIR}/hla.result.txt"
fi
echo ""
echo "=============================================="
