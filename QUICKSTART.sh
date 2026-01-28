#!/bin/bash
# Quick Start Script for HLA Typing Pipeline
# This script helps you get started quickly

echo "================================"
echo "HLA Typing Pipeline - Quick Start"
echo "================================"
echo ""

# Check Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed."
    echo "Install with: curl -s https://get.nextflow.io | bash"
    exit 1
fi

echo "✓ Nextflow found: $(nextflow -version | head -1)"
echo ""

# Check for input data
if [ ! -d "test_data" ] && [ ! -d "samples" ]; then
    echo "No input data found. Please:"
    echo "  1. Create a 'samples/' directory"
    echo "  2. Add your FASTQ or BAM files"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Ask for input type
echo "Select input type:"
echo "  1) FASTQ files"
echo "  2) BAM files"
read -p "Choice (1 or 2): " input_choice

if [ "$input_choice" == "1" ]; then
    INPUT_TYPE="fastq"
elif [ "$input_choice" == "2" ]; then
    INPUT_TYPE="bam"
else
    echo "Invalid choice"
    exit 1
fi

# Ask for tools
echo ""
echo "Select HLA typing tools (comma-separated):"
echo "  Available: optitype, arcashla, spechla, xhla, hlahd, hlala"
echo "  Example: optitype,arcashla"
read -p "Tools: " TOOLS

if [ -z "$TOOLS" ]; then
    TOOLS="optitype,arcashla"
    echo "Using default: $TOOLS"
fi

# Ask for consensus
echo ""
read -p "Enable majority voting/consensus? (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    CONSENSUS="--enable_majority_voting"
else
    CONSENSUS=""
fi

# Build command
CMD="nextflow run main.nf \
    --input samples/ \
    --input_type $INPUT_TYPE \
    --tools $TOOLS \
    $CONSENSUS \
    --outdir results/ \
    -profile docker"

echo ""
echo "================================"
echo "Running command:"
echo "$CMD"
echo "================================"
echo ""

# Ask for confirmation
read -p "Execute? (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    eval $CMD
else
    echo "Command not executed. You can run it manually:"
    echo "$CMD"
fi
