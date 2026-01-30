#!/bin/bash
# Submit Direct HLA Typing (Bypass Buggy Pipeline)

set -e

echo "=========================================="
echo "Direct HLA Typing Submission"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Count samples
SAMPLE_COUNT=$(tail -n +2 "${SAMPLE_SHEET}" | wc -l)

echo "Configuration:"
echo "  Batch: ${BATCH_NAME}"
echo "  Samples: ${SAMPLE_COUNT}"
echo "  Output: ${RESULTS_DIR}"
echo ""

# Create necessary directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOGS_DIR}"

# Update array size in script
DIRECT_SCRIPT="${SCRIPTS_DIR}/run_hla_direct_array.sh"

if [ ! -f "${DIRECT_SCRIPT}" ]; then
    echo "✗ ERROR: Script not found: ${DIRECT_SCRIPT}"
    echo "Please copy run_hla_direct_array.sh to scripts directory"
    exit 1
fi

# Create temporary script with correct array size
TEMP_SCRIPT="${DIRECT_SCRIPT}.tmp"
sed "s/#SBATCH --array=1-24/#SBATCH --array=1-${SAMPLE_COUNT}/" "${DIRECT_SCRIPT}" > "${TEMP_SCRIPT}"

echo "=========================================="
echo "Submitting Array Job"
echo "=========================================="
echo ""
echo "  Samples: ${SAMPLE_COUNT}"
echo "  CPUs per task: 8"
echo "  Memory per task: 32GB"
echo "  Time limit: 8 hours"
echo ""

# Submit the job
JOB_ID=$(sbatch --parsable "${TEMP_SCRIPT}")

if [ -z "${JOB_ID}" ]; then
    echo "✗ ERROR: Job submission failed"
    rm -f "${TEMP_SCRIPT}"
    exit 1
fi

echo "✓ Job submitted: ${JOB_ID}"
echo ""

# Cleanup
rm -f "${TEMP_SCRIPT}"

echo "=========================================="
echo "Monitoring"
echo "=========================================="
echo ""
echo "Check job status:"
echo "  squeue -j ${JOB_ID}"
echo "  squeue -u \$USER"
echo ""
echo "Monitor logs:"
echo "  ls -lh ${LOGS_DIR}/hla_direct_${JOB_ID}_*.log"
echo "  tail -f ${LOGS_DIR}/hla_direct_${JOB_ID}_1.log"
echo ""
echo "After completion, create consensus:"
echo "  bash ${SCRIPTS_DIR}/create_consensus.sh"
echo ""
