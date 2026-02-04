#!/bin/bash
# Test that all required containers and tools are available
# Usage: ./test_containers.sh [container_dir]

set -euo pipefail

CONTAINER_DIR="${1:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# If no container dir specified, try to read from user.config
if [ -z "$CONTAINER_DIR" ]; then
    USER_CONFIG="${PIPELINE_DIR}/conf/user.config"
    if [ -f "$USER_CONFIG" ]; then
        CONTAINER_DIR=$(grep -oP "container_dir\s*=\s*'\K[^']*" "$USER_CONFIG" 2>/dev/null || echo "")
    fi
fi

if [ -z "$CONTAINER_DIR" ] || [ ! -d "$CONTAINER_DIR" ]; then
    echo "Usage: $0 <container_directory>"
    echo ""
    echo "Example: $0 /scratch/project_xxx/containers"
    exit 1
fi

echo "==========================================="
echo "HLA Pipeline Container Test"
echo "==========================================="
echo "Container directory: $CONTAINER_DIR"
echo ""

# Check if singularity is available
if ! command -v singularity &> /dev/null; then
    echo "Error: Singularity is not installed"
    exit 1
fi

# Function to test a container
test_container() {
    local container=$1
    local commands=$2
    local name=$(basename "$container" .sif)

    echo "Testing $name..."

    if [ ! -f "$container" ]; then
        echo "  [MISSING] Container not found: $container"
        return 1
    fi

    local failed=0
    for cmd in $commands; do
        if singularity exec "$container" which "$cmd" &>/dev/null; then
            echo "  [OK] $cmd"
        else
            echo "  [FAIL] $cmd not found in container"
            failed=1
        fi
    done

    return $failed
}

echo ""
echo "=== Required Containers ==="
echo ""

# Track overall status
ALL_OK=1

# Test basetools container
echo "--- basetools.sif ---"
if test_container "${CONTAINER_DIR}/basetools.sif" "samtools fastqc python3 bc multiqc"; then
    echo "  Status: OK"
else
    echo "  Status: FAILED - Build with: bash scripts/build_basetools.sh"
    ALL_OK=0
fi
echo ""

# Test SpecHLA container
echo "--- spechla_1.0.7-3.sif ---"
if [ -f "${CONTAINER_DIR}/spechla_1.0.7-3.sif" ]; then
    if singularity exec "${CONTAINER_DIR}/spechla_1.0.7-3.sif" test -f /opt/SpecHLA/script/whole/SpecHLA.sh &>/dev/null; then
        echo "  [OK] SpecHLA installation found"
        if singularity exec "${CONTAINER_DIR}/spechla_1.0.7-3.sif" which samtools &>/dev/null; then
            echo "  [OK] samtools"
        else
            echo "  [WARN] samtools not in container (pipeline will use basetools for preprocessing)"
        fi
        echo "  Status: OK"
    else
        echo "  [FAIL] SpecHLA not found at /opt/SpecHLA"
        ALL_OK=0
    fi
else
    echo "  [MISSING] Container not found"
    ALL_OK=0
fi
echo ""

# Test arcasHLA container
echo "--- arcashla.sif ---"
if [ -f "${CONTAINER_DIR}/arcashla.sif" ]; then
    if test_container "${CONTAINER_DIR}/arcashla.sif" "arcasHLA samtools"; then
        echo "  Status: OK"
    else
        ALL_OK=0
    fi
else
    echo "  [MISSING] Container not found"
    ALL_OK=0
fi
echo ""

# Test OptiType container
echo "--- optitype.sif ---"
if [ -f "${CONTAINER_DIR}/optitype.sif" ]; then
    if singularity exec "${CONTAINER_DIR}/optitype.sif" which OptiTypePipeline.py &>/dev/null || \
       singularity exec "${CONTAINER_DIR}/optitype.sif" which python &>/dev/null; then
        echo "  [OK] OptiType installation found"
        if singularity exec "${CONTAINER_DIR}/optitype.sif" which samtools &>/dev/null; then
            echo "  [OK] samtools"
        else
            echo "  [WARN] samtools not in container (pipeline will use basetools for preprocessing)"
        fi
        echo "  Status: OK"
    else
        echo "  [FAIL] OptiType not properly installed"
        ALL_OK=0
    fi
else
    echo "  [MISSING] Container not found"
    ALL_OK=0
fi
echo ""

# Test HLA-HD container (optional)
echo "--- hlahd.sif (optional) ---"
if [ -f "${CONTAINER_DIR}/hlahd.sif" ]; then
    if singularity exec "${CONTAINER_DIR}/hlahd.sif" which hlahd.sh &>/dev/null; then
        echo "  [OK] HLA-HD found"
        echo "  Status: OK"
    else
        echo "  [FAIL] hlahd.sh not found"
    fi
else
    echo "  [SKIP] Container not present (optional tool)"
fi
echo ""

# Test HLA*LA container (optional)
echo "--- hlala.sif (optional) ---"
if [ -f "${CONTAINER_DIR}/hlala.sif" ]; then
    echo "  [OK] Container present"
else
    echo "  [SKIP] Container not present (optional tool)"
fi
echo ""

# Summary
echo "==========================================="
if [ $ALL_OK -eq 1 ]; then
    echo "All required containers are properly configured!"
    echo "You can run the pipeline with:"
    echo "  nextflow run main.nf -c conf/user.config -profile singularity --input_bam your_file.bam"
else
    echo "Some containers are missing or misconfigured."
    echo "Please check the messages above and fix the issues."
fi
echo "==========================================="
