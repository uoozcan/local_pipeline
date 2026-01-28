#!/bin/bash
# SpecHLA Setup Script for HPC Systems
# This script sets up SpecHLA without conda, using system-installed tools
#
# Requirements:
# - bowtie2 (or use local_bin/bowtie2)
# - bwa (or use local_bin/bwa)
# - samtools (or use local_bin/samtools)
# - freebayes (or use local_bin/freebayes)
# - Python 3 with packages: pysam, biopython, numpy, pandas, scipy, pulp
# - For SpecHap compilation: cmake, gfortran, htslib-dev, arpack (libarpack-dev)
#
# Usage:
#   bash setup_spechla.sh [--with-spechap]
#
# On CSC Puhti, load modules:
#   module load gcc cmake samtools bowtie2 bwa python-data

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "SpecHLA Setup Script"
echo "=========================================="

# Check if we should compile SpecHap
COMPILE_SPECHAP=false
if [[ "${1:-}" == "--with-spechap" ]]; then
    COMPILE_SPECHAP=true
fi

# Create config files for each HLA gene
echo "Creating HLA config files..."
for hla in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
    config_file="$SCRIPT_DIR/db/HLA/HLA_${hla}.config.txt"
    echo "bwa=$SCRIPT_DIR/db/HLA/HLA_${hla}/HLA_${hla}.fa" > "$config_file"
    echo "freebayes=$SCRIPT_DIR/db/HLA/HLA_${hla}/HLA_${hla}.fa" >> "$config_file"
    echo "blat=$SCRIPT_DIR/db/HLA/HLA_${hla}/" >> "$config_file"
done
echo "Config files created."

# Check for bowtie2 and build indexes
BOWTIE2_BUILD=""
if command -v bowtie2-build &>/dev/null; then
    BOWTIE2_BUILD="bowtie2-build"
elif [[ -x "$SCRIPT_DIR/local_bin/bowtie2-build" ]]; then
    BOWTIE2_BUILD="$SCRIPT_DIR/local_bin/bowtie2-build"
else
    echo "ERROR: bowtie2-build not found. Please install bowtie2 or load the module."
    exit 1
fi

# Build bowtie2 indexes if not present
REF1="$SCRIPT_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
REF2="$SCRIPT_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"

if [[ ! -f "${REF1}.1.bt2" ]]; then
    echo "Building bowtie2 index for HLA database (this may take a few minutes)..."
    "$BOWTIE2_BUILD" -q "$REF1" "$REF1"
    echo "Index 1 built."
else
    echo "Bowtie2 index for REF1 already exists."
fi

if [[ ! -f "${REF2}.1.bt2" ]]; then
    echo "Building bowtie2 index for HLA database v2..."
    "$BOWTIE2_BUILD" -q "$REF2" "$REF2"
    echo "Index 2 built."
else
    echo "Bowtie2 index for REF2 already exists."
fi

# Make bin files executable
chmod +x -R bin/* 2>/dev/null || true

# Compile SpecHap if requested and gfortran is available
if $COMPILE_SPECHAP; then
    echo ""
    echo "Compiling SpecHap..."

    # Check for required tools
    if ! command -v cmake &>/dev/null && [[ ! -x "$SCRIPT_DIR/local_bin/cmake" ]]; then
        echo "ERROR: cmake not found. Please install cmake or load the module."
        exit 1
    fi
    CMAKE_CMD="${CMAKE:-$(command -v cmake 2>/dev/null || echo "$SCRIPT_DIR/local_bin/cmake")}"

    if ! command -v gfortran &>/dev/null; then
        echo "ERROR: gfortran not found. SpecHap requires ARPACK which needs Fortran."
        echo "On HPC, try: module load gcc"
        exit 1
    fi

    # Find htslib
    HTSLIB_DIR=""
    if [[ -d "$SCRIPT_DIR/local_tools/samtools-1.17/htslib-1.17" ]]; then
        HTSLIB_DIR="$SCRIPT_DIR/local_tools/samtools-1.17/htslib-1.17"
    elif pkg-config --exists htslib 2>/dev/null; then
        HTSLIB_DIR="$(pkg-config --variable=includedir htslib)"
    fi

    # Compile SpecHap
    mkdir -p "$SCRIPT_DIR/bin/SpecHap/build"
    cd "$SCRIPT_DIR/bin/SpecHap/build"

    if [[ -n "$HTSLIB_DIR" ]]; then
        "$CMAKE_CMD" .. -DHTSlib_INCLUDE_DIR="$HTSLIB_DIR" 2>&1 | tail -10
    else
        "$CMAKE_CMD" .. 2>&1 | tail -10
    fi

    make -j4 2>&1 | tail -10

    if [[ -f "$SCRIPT_DIR/bin/SpecHap/build/SpecHap" ]]; then
        echo "SpecHap compiled successfully!"
    else
        echo "WARNING: SpecHap compilation may have failed."
    fi

    # Compile extractHairs
    mkdir -p "$SCRIPT_DIR/bin/extractHairs/build"
    cd "$SCRIPT_DIR/bin/extractHairs/build"
    "$CMAKE_CMD" .. 2>&1 | tail -5
    make -j4 2>&1 | tail -5

    cd "$SCRIPT_DIR"
fi

echo ""
echo "=========================================="
echo "SpecHLA setup complete!"
echo "=========================================="
echo ""
echo "To use SpecHLA, add local_bin to your PATH:"
echo "  export PATH=\"$SCRIPT_DIR/local_bin:\$PATH\""
echo ""
echo "Then run:"
echo "  bash script/whole/SpecHLA.sh -h"
echo ""
if [[ ! -f "$SCRIPT_DIR/bin/SpecHap/build/SpecHap" ]]; then
    echo "NOTE: SpecHap is not compiled. For full-resolution typing, run:"
    echo "  bash setup_spechla.sh --with-spechap"
    echo "This requires gfortran and ARPACK."
fi
