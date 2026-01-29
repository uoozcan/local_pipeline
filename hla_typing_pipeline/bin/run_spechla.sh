#!/bin/bash
#
# SpecHLA wrapper script for BAM files
# Uses Singularity container for execution
#
# Usage: ./run_spechla.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-r REFERENCE] [-j THREADS]
#

set -euo pipefail

# Default values
THREADS=8
REFERENCE="hg38"
KEEP_TEMP=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Detect paths (can be overridden)
CONTAINER="${CONTAINER:-${SCRIPT_DIR}/../../hla_references/containers/spechla_1.0.7-3.sif}"

usage() {
    cat << EOF
SpecHLA HLA Typing Wrapper (Container-based)

Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]

Required arguments:
  -n    Sample name (e.g., FH_7087)
  -b    Input BAM file (sorted and indexed)
  -o    Output directory

Optional arguments:
  -r    Reference genome: hg38 or hg19 (default: hg38)
  -c    Singularity container path (default: auto-detect)
  -j    Number of threads (default: 8)
  -k    Keep temporary files (default: remove)
  -h    Show this help message

Environment variables:
  CONTAINER    Path to SpecHLA Singularity container

Example:
  $0 -n Sample1 -b Sample1.bam -o ./results -r hg38 -j 8

EOF
    exit 1
}

while getopts "n:b:o:r:c:j:kh" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        c) CONTAINER="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        k) KEEP_TEMP=1 ;;
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
    echo "Please ensure the SpecHLA container is available."
    exit 1
fi

# Validate reference
if [ "$REFERENCE" != "hg38" ] && [ "$REFERENCE" != "hg19" ]; then
    echo "Error: Reference must be 'hg38' or 'hg19'"
    exit 1
fi

# Get absolute paths
BAM_FILE=$(realpath "$BAM_FILE")
BAM_DIR=$(dirname "$BAM_FILE")
BAM_NAME=$(basename "$BAM_FILE")

# Check for BAM index
if [ ! -f "${BAM_FILE}.bai" ] && [ ! -f "${BAM_FILE%.bam}.bai" ]; then
    echo "Warning: BAM index not found. Creating index..."
    singularity exec "$CONTAINER" samtools index "$BAM_FILE"
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "SpecHLA HLA Typing Pipeline (Container)"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Reference:  $REFERENCE"
echo "Threads:    $THREADS"
echo "Container:  $CONTAINER"
echo "=============================================="
echo ""

# Define HLA region based on reference genome
if [ "$REFERENCE" == "hg38" ]; then
    # Check chromosome naming convention
    CHR_PREFIX=$(singularity exec --bind "${BAM_DIR}:/bam_dir" "$CONTAINER" \
        samtools view -H "/bam_dir/${BAM_NAME}" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28510120-33480577"
    else
        HLA_REGION="6:28510120-33480577"
    fi
else
    CHR_PREFIX=$(singularity exec --bind "${BAM_DIR}:/bam_dir" "$CONTAINER" \
        samtools view -H "/bam_dir/${BAM_NAME}" | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    if [ -n "$CHR_PREFIX" ]; then
        HLA_REGION="chr6:28477797-33448354"
    else
        HLA_REGION="6:28477797-33448354"
    fi
fi

# Step 1: Extract HLA reads
echo "[Step 1/4] Extracting HLA reads from $HLA_REGION..."
singularity exec \
    --bind "${BAM_DIR}:/bam_dir" \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    bash -c "
        samtools view -@ ${THREADS} -b /bam_dir/${BAM_NAME} ${HLA_REGION} > /output/hla_extract.bam
        samtools index /output/hla_extract.bam
    "

READ_COUNT=$(singularity exec --bind "${WORK_DIR}:/output" "$CONTAINER" \
    samtools view -c /output/hla_extract.bam)
echo "  Extracted $READ_COUNT reads from HLA region"

if [ "$READ_COUNT" -lt 1000 ]; then
    echo "Warning: Low number of reads extracted. HLA typing results may be unreliable."
fi

# Step 2: Convert BAM to FASTQ
echo ""
echo "[Step 2/4] Converting BAM to FASTQ..."
singularity exec \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    bash -c "
        samtools sort -n -@ ${THREADS} /output/hla_extract.bam -o /output/namesort.bam
        samtools fastq -@ ${THREADS} \
            -1 /output/${SAMPLE_NAME}_R1.fastq.gz \
            -2 /output/${SAMPLE_NAME}_R2.fastq.gz \
            -0 /dev/null -s /dev/null \
            /output/namesort.bam
    "
echo "  Created FASTQ files"

# Step 3: Run SpecHLA
echo ""
echo "[Step 3/4] Running SpecHLA..."
singularity exec \
    --bind "${WORK_DIR}:/output" \
    "$CONTAINER" \
    bash -c "
        cd /output
        bash /opt/SpecHLA/script/whole/SpecHLA.sh \
            -n ${SAMPLE_NAME} \
            -1 ${SAMPLE_NAME}_R1.fastq.gz \
            -2 ${SAMPLE_NAME}_R2.fastq.gz \
            -o . \
            -j ${THREADS} \
            -u 1
    "

# Step 4: Parse results
echo ""
echo "[Step 4/4] Formatting results..."

# Find result file
RESULT_FILE=""
if [ -f "${WORK_DIR}/hla.result.txt" ]; then
    RESULT_FILE="${WORK_DIR}/hla.result.txt"
elif [ -f "${WORK_DIR}/${SAMPLE_NAME}/hla.result.txt" ]; then
    RESULT_FILE="${WORK_DIR}/${SAMPLE_NAME}/hla.result.txt"
fi

OUTPUT_FILE="${WORK_DIR}/${SAMPLE_NAME}_spechla.txt"
if [ -n "$RESULT_FILE" ]; then
    cp "$RESULT_FILE" "$OUTPUT_FILE"
    echo "  Results saved to: $OUTPUT_FILE"
else
    echo "# SpecHLA results for ${SAMPLE_NAME}" > "$OUTPUT_FILE"
    echo "# Error: No results generated" >> "$OUTPUT_FILE"
    echo "  Warning: No results file found"
fi

# Cleanup temporary files
if [ "$KEEP_TEMP" -eq 0 ]; then
    echo ""
    echo "Cleaning up temporary files..."
    rm -f "${WORK_DIR}/hla_extract.bam" "${WORK_DIR}/hla_extract.bam.bai"
    rm -f "${WORK_DIR}/namesort.bam"
    rm -f "${WORK_DIR}/${SAMPLE_NAME}_R1.fastq.gz" "${WORK_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
fi

echo ""
echo "=============================================="
echo "SpecHLA Complete!"
echo "=============================================="
echo "Results: $OUTPUT_FILE"
echo ""

if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    echo "HLA Typing Results:"
    echo "-------------------"
    cat "$OUTPUT_FILE"
fi
echo "=============================================="
