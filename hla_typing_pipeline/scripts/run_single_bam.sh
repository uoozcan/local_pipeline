#!/bin/bash
# Run HLA typing pipeline on a single BAM file
# Usage: ./run_single_bam.sh <input.bam> [output_dir] [tools]
#
# Examples:
#   ./run_single_bam.sh sample.bam
#   ./run_single_bam.sh sample.bam ./results spechla,arcashla,optitype

set -euo pipefail

# Parse arguments
INPUT_BAM="${1:-}"
OUTPUT_DIR="${2:-./results}"
TOOLS="${3:-spechla,arcashla,optitype}"

if [ -z "$INPUT_BAM" ]; then
    echo "Usage: $0 <input.bam> [output_dir] [tools]"
    echo ""
    echo "Arguments:"
    echo "  input.bam     Path to input BAM file (required)"
    echo "  output_dir    Output directory (default: ./results)"
    echo "  tools         Comma-separated list of tools (default: spechla,arcashla,optitype)"
    echo ""
    echo "Available tools: spechla, hlahd, arcashla, optitype, hlala, xhla"
    echo ""
    echo "Examples:"
    echo "  $0 sample.bam"
    echo "  $0 sample.bam ./my_results spechla,optitype"
    exit 1
fi

# Get absolute path
INPUT_BAM="$(realpath "$INPUT_BAM")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Check if BAM file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

# Check for user config
USER_CONFIG="${PIPELINE_DIR}/conf/user.config"
if [ ! -f "$USER_CONFIG" ]; then
    echo "Warning: user.config not found at $USER_CONFIG"
    echo "Copying template..."
    cp "${PIPELINE_DIR}/conf/user.config.template" "$USER_CONFIG"
    echo "Please edit $USER_CONFIG with your settings and run again."
    exit 1
fi

# Extract sample name from BAM file
SAMPLE_NAME=$(basename "$INPUT_BAM" .bam)
SAMPLE_NAME="${SAMPLE_NAME%.sorted}"
SAMPLE_NAME="${SAMPLE_NAME%.dedup}"

echo "==========================================="
echo "HLA Typing Pipeline - Single BAM"
echo "==========================================="
echo "Input BAM:    $INPUT_BAM"
echo "Sample name:  $SAMPLE_NAME"
echo "Output dir:   $OUTPUT_DIR"
echo "Tools:        $TOOLS"
echo "==========================================="

# Check if nextflow is available
if ! command -v nextflow &> /dev/null; then
    echo "Error: nextflow is not installed or not in PATH"
    echo "Install with: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# Determine profile to use
PROFILE="singularity"
if [ -n "${NXF_SINGULARITY_CACHEDIR:-}" ] || command -v singularity &> /dev/null; then
    PROFILE="singularity"
elif command -v docker &> /dev/null; then
    PROFILE="docker"
else
    echo "Warning: Neither Singularity nor Docker found. Trying local execution..."
    PROFILE="local"
fi

echo "Using profile: $PROFILE"
echo ""

# Run the pipeline
cd "$PIPELINE_DIR"
nextflow run main.nf \
    -c conf/user.config \
    -profile "$PROFILE" \
    --input_bam "$INPUT_BAM" \
    --outdir "$OUTPUT_DIR" \
    --tools "$TOOLS" \
    -resume

echo ""
echo "==========================================="
echo "Pipeline completed!"
echo "Results are in: $OUTPUT_DIR/$SAMPLE_NAME"
echo "==========================================="
