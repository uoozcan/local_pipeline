#!/bin/bash
#
# HLA Pipeline - Quick Fix Script
# Run this script to fix common issues
#

echo "=========================================="
echo "HLA Pipeline Quick Fix"
echo "=========================================="
echo ""

# 1. Clean up previous run artifacts
echo "Step 1: Cleaning up previous run artifacts..."
rm -f ./results/pipeline_info/execution_trace.txt
rm -f ./results/pipeline_info/execution_timeline.html
rm -f ./results/pipeline_info/execution_report.html
rm -f ./results/pipeline_info/pipeline_dag.html
echo "✓ Cleaned up old trace files"
echo ""

# 2. Clean work directory (optional - comment out if you want to keep cache)
# echo "Step 2: Cleaning work directory..."
# rm -rf work/
# echo "✓ Cleaned work directory"
# echo ""

# 3. Verify module files exist
echo "Step 2: Verifying module files..."
MISSING=0

if [ ! -f "modules/optitype.nf" ]; then
    echo "✗ Missing: modules/optitype.nf"
    MISSING=1
fi

if [ ! -f "modules/arcashla.nf" ]; then
    echo "✗ Missing: modules/arcashla.nf"
    MISSING=1
fi

if [ ! -f "modules/spechla.nf" ]; then
    echo "✗ Missing: modules/spechla.nf"
    MISSING=1
fi

if [ ! -f "modules/bam_to_fastq.nf" ]; then
    echo "✗ Missing: modules/bam_to_fastq.nf"
    MISSING=1
fi

if [ ! -f "modules/aggregation.nf" ]; then
    echo "✗ Missing: modules/aggregation.nf"
    MISSING=1
fi

if [ ! -f "modules/majority_voting.nf" ]; then
    echo "✗ Missing: modules/majority_voting.nf"
    MISSING=1
fi

if [ $MISSING -eq 0 ]; then
    echo "✓ All module files present"
else
    echo ""
    echo "ERROR: Some module files are missing!"
    echo "Please ensure you extracted all files from the tar.gz archive"
    exit 1
fi
echo ""

# 4. Verify container files exist
echo "Step 3: Verifying Singularity containers..."
CONTAINER_DIR="/scratch/project_2008084/hla_references/singularity_cache/containers"

if [ -f "${CONTAINER_DIR}/optitype.sif" ]; then
    echo "✓ optitype.sif found"
else
    echo "✗ optitype.sif NOT found at ${CONTAINER_DIR}"
fi

if [ -f "${CONTAINER_DIR}/arcashla.sif" ]; then
    echo "✓ arcashla.sif found"
else
    echo "✗ arcashla.sif NOT found at ${CONTAINER_DIR}"
fi

if [ -f "${CONTAINER_DIR}/spechla.sif" ]; then
    echo "✓ spechla.sif found"
else
    echo "✗ spechla.sif NOT found at ${CONTAINER_DIR}"
fi
echo ""

# 5. Check main files
echo "Step 4: Verifying main pipeline files..."
if [ -f "main.nf" ]; then
    echo "✓ main.nf found"
else
    echo "✗ main.nf NOT found"
    exit 1
fi

if [ -f "nextflow.config" ]; then
    echo "✓ nextflow.config found"
else
    echo "✗ nextflow.config NOT found"
    exit 1
fi
echo ""

# 6. Make scripts executable
echo "Step 5: Making scripts executable..."
chmod +x submit_slurm.sh 2>/dev/null && echo "✓ submit_slurm.sh is executable" || echo "✗ Could not make submit_slurm.sh executable"
echo ""

# 7. Show Nextflow version
echo "Step 6: Checking Nextflow..."
if command -v nextflow &> /dev/null; then
    NEXTFLOW_VERSION=$(nextflow -version 2>&1 | head -1)
    echo "✓ Nextflow found: $NEXTFLOW_VERSION"
else
    echo "✗ Nextflow not found!"
    echo "  Load module: module load nextflow/25.10.0"
fi
echo ""

echo "=========================================="
echo "Fix script completed!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Edit submit_slurm.sh with your data paths"
echo "2. Run: sbatch submit_slurm.sh"
echo ""
echo "To test directly:"
echo "nextflow run main.nf --help"
echo ""
