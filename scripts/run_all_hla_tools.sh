#!/bin/bash
#
# Master script to run multiple HLA typing tools on a BAM file
# Runs: SpecHLA, HLA-HD, HLAscan, OptiType, arcasHLA, HLA*LA
#
# Usage: ./run_all_hla_tools.sh -n SAMPLE_NAME -b BAM_FILE -o OUTPUT_DIR [OPTIONS]
#

set -e

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Default values
THREADS=8
REFERENCE="hg38"
SEQ_TYPE="dna"
TOOLS="spechla,hlahd,optitype"  # Default tools to run

usage() {
    echo "Master HLA Typing Pipeline - Run multiple tools"
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
    echo "  -t    Sequencing type: dna or rna (default: dna)"
    echo "  -T    Tools to run (comma-separated, default: spechla,hlahd,optitype)"
    echo "        Available: spechla, hlahd, hlascan, optitype, arcashla, hlala, all"
    echo "  -h    Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Run default tools (SpecHLA, HLA-HD, OptiType)"
    echo "  $0 -n SAMPLE -b sample.bam -o ./results"
    echo ""
    echo "  # Run all tools"
    echo "  $0 -n SAMPLE -b sample.bam -o ./results -T all"
    echo ""
    echo "  # Run specific tools"
    echo "  $0 -n SAMPLE -b sample.bam -o ./results -T spechla,optitype,hlahd"
    echo ""
    exit 1
}

while getopts "n:b:o:r:j:t:T:h" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        b) BAM_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        t) SEQ_TYPE="$OPTARG" ;;
        T) TOOLS="$OPTARG" ;;
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

# Expand "all" to all tools
if [ "$TOOLS" == "all" ]; then
    TOOLS="spechla,hlahd,hlascan,optitype,arcashla,hlala"
fi

# Get absolute path
BAM_FILE=$(realpath "$BAM_FILE")

# Create output directory
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

echo "=============================================="
echo "HLA Typing Pipeline - Multiple Tools"
echo "=============================================="
echo "Sample:     $SAMPLE_NAME"
echo "BAM file:   $BAM_FILE"
echo "Output:     $OUTPUT_DIR"
echo "Reference:  $REFERENCE"
echo "Seq Type:   $SEQ_TYPE"
echo "Threads:    $THREADS"
echo "Tools:      $TOOLS"
echo "=============================================="
echo ""

# Parse tools
IFS=',' read -ra TOOL_ARRAY <<< "$TOOLS"

# Track results
declare -A RESULTS

# Run each tool
for tool in "${TOOL_ARRAY[@]}"; do
    tool=$(echo "$tool" | tr '[:upper:]' '[:lower:]' | xargs)  # lowercase and trim

    echo ""
    echo "======================================================"
    echo "Running: $tool"
    echo "======================================================"

    case "$tool" in
        spechla)
            if bash "${PROJECT_DIR}/spechla_local/run_spechla_bam.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/spechla" \
                -r "$REFERENCE" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        hlahd)
            if bash "${SCRIPT_DIR}/run_hlahd.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/hlahd" \
                -r "$REFERENCE" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        hlascan)
            if bash "${SCRIPT_DIR}/run_hlascan.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/hlascan" \
                -r "$REFERENCE" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        optitype)
            if bash "${SCRIPT_DIR}/run_optitype.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/optitype" \
                -t "$SEQ_TYPE" \
                -r "$REFERENCE" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        arcashla)
            if bash "${SCRIPT_DIR}/run_arcashla.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/arcashla" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        hlala)
            if bash "${SCRIPT_DIR}/run_hlala.sh" \
                -n "$SAMPLE_NAME" \
                -b "$BAM_FILE" \
                -o "${OUTPUT_DIR}/hlala" \
                -j "$THREADS" 2>&1; then
                RESULTS[$tool]="SUCCESS"
            else
                RESULTS[$tool]="FAILED"
            fi
            ;;

        *)
            echo "Warning: Unknown tool '$tool', skipping..."
            RESULTS[$tool]="SKIPPED"
            ;;
    esac
done

# Summary
echo ""
echo "======================================================"
echo "HLA Typing Complete - Summary"
echo "======================================================"
echo ""
echo "Tool         Status"
echo "------------ --------"
for tool in "${TOOL_ARRAY[@]}"; do
    tool=$(echo "$tool" | tr '[:upper:]' '[:lower:]' | xargs)
    printf "%-12s %s\n" "$tool" "${RESULTS[$tool]:-UNKNOWN}"
done
echo ""
echo "Results saved to: $OUTPUT_DIR"
echo "======================================================"

# Create combined results summary
SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_hla_summary.txt"
{
    echo "# HLA Typing Summary for ${SAMPLE_NAME}"
    echo "# Date: $(date)"
    echo "# BAM: $BAM_FILE"
    echo ""

    # SpecHLA results
    if [ -f "${OUTPUT_DIR}/spechla/${SAMPLE_NAME}/hla.result.txt" ]; then
        echo "## SpecHLA Results"
        cat "${OUTPUT_DIR}/spechla/${SAMPLE_NAME}/hla.result.txt"
        echo ""
    fi

    # HLA-HD results
    if [ -f "${OUTPUT_DIR}/hlahd/${SAMPLE_NAME}/result/${SAMPLE_NAME}_final.result.txt" ]; then
        echo "## HLA-HD Results"
        cat "${OUTPUT_DIR}/hlahd/${SAMPLE_NAME}/result/${SAMPLE_NAME}_final.result.txt"
        echo ""
    fi

    # OptiType results
    OPTI_RESULT=$(find "${OUTPUT_DIR}/optitype/${SAMPLE_NAME}" -name "*_result.tsv" 2>/dev/null | head -1)
    if [ -n "$OPTI_RESULT" ] && [ -f "$OPTI_RESULT" ]; then
        echo "## OptiType Results (Class I only)"
        cat "$OPTI_RESULT"
        echo ""
    fi

    # HLAscan results
    if [ -f "${OUTPUT_DIR}/hlascan/${SAMPLE_NAME}/${SAMPLE_NAME}.hlascan.result.txt" ]; then
        echo "## HLAscan Results"
        cat "${OUTPUT_DIR}/hlascan/${SAMPLE_NAME}/${SAMPLE_NAME}.hlascan.result.txt"
        echo ""
    fi

} > "$SUMMARY_FILE"

echo ""
echo "Combined summary saved to: $SUMMARY_FILE"
