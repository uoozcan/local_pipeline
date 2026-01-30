#!/bin/bash
#
# HLAscan wrapper script for BAM files
# HLAscan performs HLA typing from WGS/WES data
#
# Usage: ./run_hlascan.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [-r REFERENCE] [-j THREADS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONTAINER="${PROJECT_DIR}/hla_references/containers/hlascan.sif"
LOCAL_BIN="${PROJECT_DIR}/spechla_local/local_bin"
HLASCAN_DB="${PROJECT_DIR}/hla_references/databases/hlascan_db"

# Default values
THREADS=8
REFERENCE="hg38"
GENES="HLA-A,HLA-B,HLA-C,HLA-DPA1,HLA-DPB1,HLA-DQA1,HLA-DQB1,HLA-DRB1"

usage() {
    echo "HLAscan wrapper for BAM files"
    echo ""
    echo "Usage: $0 -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]"
    echo ""
    echo "Required arguments:"
    echo "  -n    Sample name"
    echo "  -b    Input BAM file (sorted and indexed)"
    echo "  -o    Output directory"
    echo ""
    echo "Optional arguments:"
    echo "  -r    Reference genome: hg38 (38) or hg19 (37) (default: hg38)"
    echo "  -j    Number of threads (default: 8)"
    echo "  -g    HLA genes to type (comma-separated, default: HLA-A,HLA-B,HLA-C,...)"
    echo "  -h    Show this help message"
    echo ""
    exit 1
}

while getopts "n:b:o:r:j:g:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
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

# Convert reference to version number
if [ "$REFERENCE" == "hg38" ]; then
    REF_VERSION="38"
elif [ "$REFERENCE" == "hg19" ]; then
    REF_VERSION="37"
else
    REF_VERSION="$REFERENCE"
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

# Check database exists
if [ ! -d "$HLASCAN_DB/IMGT/HLA" ]; then
    echo "Error: HLAscan database not found: $HLASCAN_DB/IMGT/HLA"
    echo "Please run: ./scripts/setup_hla_databases.sh"
    exit 1
fi

# Get absolute path of BAM file
BAM_FILE=$(realpath "$BAM_FILE")
BAM_DIR=$(dirname "$BAM_FILE")
BAM_NAME=$(basename "$BAM_FILE")

# Create output directory
mkdir -p "${OUTPUT_DIR}/${SAMPLE_NAME}"
WORK_DIR=$(realpath "${OUTPUT_DIR}/${SAMPLE_NAME}")

echo "=============================================="
echo "HLAscan HLA Typing Pipeline"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $WORK_DIR"
echo "Reference:  $REF_VERSION"
echo "Threads:    $THREADS"
echo "Genes:      $GENES"
echo "=============================================="
echo ""

# Run HLAscan for each gene
echo "[Running HLAscan...]"

# Convert comma-separated genes to array
IFS=',' read -ra GENE_ARRAY <<< "$GENES"

for gene in "${GENE_ARRAY[@]}"; do
    echo "  Typing $gene..."

    singularity exec \
        --bind "${BAM_DIR}:/bam_dir" \
        --bind "${WORK_DIR}:/output" \
        --bind "${HLASCAN_DB}:/database" \
        "$CONTAINER" \
        /usr/local/bin/hlascan \
        -b "/bam_dir/${BAM_NAME}" \
        -d /database/IMGT/HLA \
        -v "$REF_VERSION" \
        -g "$gene" \
        -t "$THREADS" \
        > "${WORK_DIR}/${SAMPLE_NAME}.${gene}.txt" 2>&1 || true
done

# Combine results
echo ""
echo "[Combining results...]"

{
    echo "# HLAscan Results for ${SAMPLE_NAME}"
    echo "# Reference: GRCh${REF_VERSION}"
    echo ""
    echo "Gene	Allele1	Allele2	Score1	Score2"

    for gene in "${GENE_ARRAY[@]}"; do
        if [ -f "${WORK_DIR}/${SAMPLE_NAME}.${gene}.txt" ]; then
            # Parse HLAscan output
            result=$(grep -E "^${gene}" "${WORK_DIR}/${SAMPLE_NAME}.${gene}.txt" 2>/dev/null | head -1 || echo "${gene}	-	-	-	-")
            if [ -n "$result" ]; then
                echo "$result"
            else
                echo "${gene}	-	-	-	-"
            fi
        else
            echo "${gene}	-	-	-	-"
        fi
    done
} > "${WORK_DIR}/${SAMPLE_NAME}.hlascan.result.txt"

echo ""
echo "=============================================="
echo "HLAscan Complete!"
echo "=============================================="
echo "Results: ${WORK_DIR}/${SAMPLE_NAME}.hlascan.result.txt"
echo ""
cat "${WORK_DIR}/${SAMPLE_NAME}.hlascan.result.txt"
echo "=============================================="
