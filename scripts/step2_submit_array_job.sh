#!/bin/bash

# STEP 2: Master Script - BAM Preparation (Array Job)
# This script submits both the array job and aggregation with proper dependencies

set -e

echo "=========================================="
echo "BAM Preparation Pipeline Submission"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Create logs directory
mkdir -p "${LOGS_DIR}"

# Count BAM files to determine array size
echo "Counting BAM files..."
BAM_COUNT=$(find ${RAW_BAM_DIR} -name "*_tumor.bam" -type f | wc -l)

if [ ${BAM_COUNT} -eq 0 ]; then
    echo "ERROR: No BAM files found in ${RAW_BAM_DIR}"
    echo ""
    echo "Please run step1_transfer_bam.sh first"
    exit 1
fi

echo "Found ${BAM_COUNT} BAM files"
echo ""

# Update array size in the script
ARRAY_SCRIPT="${SCRIPTS_DIR}/step2_prepare_bam_array.sh"
if [ -f "${ARRAY_SCRIPT}" ]; then
    # Temporarily update array size (create a copy)
    sed "s/#SBATCH --array=1-100%20/#SBATCH --array=1-${BAM_COUNT}%20/" "${ARRAY_SCRIPT}" > "${ARRAY_SCRIPT}.tmp"
    ARRAY_SCRIPT_SUBMIT="${ARRAY_SCRIPT}.tmp"
else
    echo "ERROR: Array script not found: ${ARRAY_SCRIPT}"
    exit 1
fi

# Submit array job
echo "Submitting array job for ${BAM_COUNT} samples..."
echo "  Resources per task: 8 CPUs, 32GB RAM"
echo "  Max concurrent tasks: 20"
echo ""

ARRAY_JOB_ID=$(sbatch --parsable "${ARRAY_SCRIPT_SUBMIT}")

if [ -z "${ARRAY_JOB_ID}" ]; then
    echo "ERROR: Failed to submit array job"
    rm -f "${ARRAY_SCRIPT_SUBMIT}"
    exit 1
fi

echo "✓ Array job submitted: Job ID ${ARRAY_JOB_ID}"
echo ""

# Submit aggregation job with dependency
echo "Submitting aggregation job (depends on array job)..."
AGG_SCRIPT="${SCRIPTS_DIR}/step2b_aggregate_results.sh"
AGG_JOB_ID=$(sbatch --parsable --dependency=afterok:${ARRAY_JOB_ID} "${AGG_SCRIPT}")

if [ -z "${AGG_JOB_ID}" ]; then
    echo "WARNING: Failed to submit aggregation job"
    echo "You can manually submit it after array job completes:"
    echo "  sbatch ${AGG_SCRIPT}"
else
    echo "✓ Aggregation job submitted: Job ID ${AGG_JOB_ID}"
fi

echo ""
echo "=========================================="
echo "Submission Complete"
echo "=========================================="
echo ""
echo "Job Chain:"
echo "  1. Array Job: ${ARRAY_JOB_ID} (${BAM_COUNT} tasks)"
echo "  2. Aggregation: ${AGG_JOB_ID} (runs after array completes)"
echo ""
echo "Monitor jobs:"
echo "  squeue -u \$USER"
echo "  squeue -j ${ARRAY_JOB_ID}"
echo ""
echo "Check progress:"
echo "  ls -lh ${LOGS_DIR}/prepare_bam_array_${ARRAY_JOB_ID}_*.log"
echo "  tail -f ${LOGS_DIR}/prepare_bam_array_${ARRAY_JOB_ID}_1.log"
echo ""
echo "After completion:"
echo "  ls ${INPUT_DIR}/samplesheet*.csv"
echo "  cat ${INPUT_DIR}/preparation_summary.txt"
echo ""
echo "Cancel jobs if needed:"
echo "  scancel ${ARRAY_JOB_ID}  # Cancel array job"
echo "  scancel ${AGG_JOB_ID}    # Cancel aggregation"
echo ""

# Cleanup temporary script
rm -f "${ARRAY_SCRIPT_SUBMIT}"
