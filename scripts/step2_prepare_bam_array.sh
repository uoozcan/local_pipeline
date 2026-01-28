#!/bin/bash
#SBATCH --job-name=prep_bam_array
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --array=1-100%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/prepare_bam_array_%A_%a.log
#SBATCH --error=logs/prepare_bam_array_%A_%a.err

# STEP 2: Prepare Batch BAM Files for Pipeline (Array Job Version)
# Purpose: Validate BAM files, create indices in parallel
# Each array task processes one BAM file

set -e

echo "=========================================="
echo "Array Task: Preparing BAM File for HLA Typing"
echo "=========================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Load required modules
echo "Loading modules..."
module load samtools
echo "✓ Modules loaded"
echo ""

# Create necessary directories
mkdir -p "${INPUT_DIR}"
mkdir -p "${INPUT_DIR}/array_results"
mkdir -p "${LOGS_DIR}"

# Get list of BAM files
BAM_LIST_FILE="${INPUT_DIR}/bam_file_list.txt"

# Create BAM list if this is task 1 or if it doesn't exist
if [ ! -f "${BAM_LIST_FILE}" ] || [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]; then
    find ${RAW_BAM_DIR} -name "*_tumor.bam" -type f | sort > "${BAM_LIST_FILE}"
    total_bams=$(cat "${BAM_LIST_FILE}" | wc -l)
    echo "Found ${total_bams} BAM files total"
fi

# Wait a moment for file system sync if this is not task 1
if [ ${SLURM_ARRAY_TASK_ID} -ne 1 ]; then
    sleep 2
fi

# Check if BAM list exists
if [ ! -f "${BAM_LIST_FILE}" ]; then
    echo "ERROR: BAM list file not found: ${BAM_LIST_FILE}"
    exit 1
fi

# Get the BAM file for this array task
BAM_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST_FILE}")

if [ -z "${BAM_PATH}" ]; then
    echo "No BAM file assigned to task ${SLURM_ARRAY_TASK_ID}"
    echo "This task has no work to do (probably array size > number of files)"
    exit 0
fi

# Extract sample ID
SAMPLE_ID=$(basename $(dirname ${BAM_PATH}))

echo "=========================================="
echo "Processing Sample: ${SAMPLE_ID}"
echo "=========================================="
echo "BAM file: ${BAM_PATH}"
echo ""

# Output file for this task's results
TASK_RESULT="${INPUT_DIR}/array_results/task_${SLURM_ARRAY_TASK_ID}.txt"
TASK_STATUS="${INPUT_DIR}/array_results/status_${SLURM_ARRAY_TASK_ID}.txt"

# Function to log status
log_status() {
    local status=$1
    local message=$2
    echo "${status}|${SAMPLE_ID}|${BAM_PATH}|${message}" > "${TASK_STATUS}"
}

# Check if BAM exists
if [ ! -f "${BAM_PATH}" ]; then
    echo "✗ ERROR: BAM file not found: ${BAM_PATH}"
    log_status "FAILED" "BAM file not found"
    exit 1
fi

echo "✓ BAM file found"
echo "  Size: $(du -h ${BAM_PATH} | cut -f1)"
echo ""

# Check/create BAM index with threading
BAI_PATH="${BAM_PATH}.bai"
if [ ! -f "${BAI_PATH}" ]; then
    echo "Creating BAM index (using ${SLURM_CPUS_PER_TASK} threads)..."
    if samtools index -@ ${SLURM_CPUS_PER_TASK} "${BAM_PATH}"; then
        echo "  ✓ Index created successfully"
    else
        echo "  ✗ ERROR: Index creation failed"
        log_status "FAILED" "Index creation failed"
        exit 1
    fi
else
    echo "✓ BAM index already exists"
fi
echo ""

# Validate BAM structure
echo "Validating BAM structure..."
if samtools quickcheck "${BAM_PATH}"; then
    echo "  ✓ BAM file is valid"
else
    echo "  ✗ ERROR: BAM file is corrupted"
    log_status "FAILED" "BAM file corrupted"
    exit 1
fi
echo ""

# Detect chromosome naming convention
echo "Detecting chromosome naming convention..."
CHR_WITH_PREFIX=$(samtools view -H "${BAM_PATH}" | grep -c "^@SQ.*SN:chr" || echo "0")

if [ ${CHR_WITH_PREFIX} -gt 0 ]; then
    CHR6_NAME="chr6"
    echo "  Chromosomes: WITH 'chr' prefix (${CHR6_NAME})"
else
    CHR6_NAME="6"
    echo "  Chromosomes: WITHOUT 'chr' prefix (${CHR6_NAME})"
fi
echo ""

# Check HLA region coverage (using threads for faster processing)
echo "Analyzing HLA region coverage (${CHR6_NAME}:28510120-33480577)..."
HLA_READS=$(samtools view -@ ${SLURM_CPUS_PER_TASK} -c "${BAM_PATH}" ${CHR6_NAME}:28510120-33480577 2>/dev/null || echo "0")

if [ "${HLA_READS}" -gt 1000 ]; then
    echo "  ✓ Found ${HLA_READS} reads in HLA region (GOOD)"
    HLA_STATUS="GOOD"
elif [ "${HLA_READS}" -gt 100 ]; then
    echo "  ⚠ Found ${HLA_READS} reads in HLA region (LOW)"
    echo "    Pipeline may produce low-confidence results"
    HLA_STATUS="LOW"
else
    echo "  ✗ Found only ${HLA_READS} reads in HLA region (VERY LOW)"
    echo "    Results will likely be unreliable"
    HLA_STATUS="VERY_LOW"
fi
echo ""

# Get comprehensive BAM statistics (parallel processing)
echo "Computing BAM statistics..."
TOTAL_READS=$(samtools view -@ ${SLURM_CPUS_PER_TASK} -c "${BAM_PATH}" 2>/dev/null || echo "0")
MAPPED_READS=$(samtools view -@ ${SLURM_CPUS_PER_TASK} -c -F 4 "${BAM_PATH}" 2>/dev/null || echo "0")

echo "  Total reads: ${TOTAL_READS}"
echo "  Mapped reads: ${MAPPED_READS}"

if [ ${TOTAL_READS} -gt 0 ]; then
    HLA_PERCENT=$(echo "scale=4; ${HLA_READS} * 100 / ${TOTAL_READS}" | bc)
    MAPPED_PERCENT=$(echo "scale=2; ${MAPPED_READS} * 100 / ${TOTAL_READS}" | bc)
    echo "  Mapping rate: ${MAPPED_PERCENT}%"
    echo "  HLA reads: ${HLA_PERCENT}% of total"
else
    echo "  ✗ ERROR: No reads found in BAM file"
    log_status "FAILED" "Empty BAM file"
    exit 1
fi
echo ""

# Write results for this sample
cat > "${TASK_RESULT}" << EOF
sample,bam,bai,total_reads,mapped_reads,hla_reads,hla_status
${SAMPLE_ID},${BAM_PATH},${BAI_PATH},${TOTAL_READS},${MAPPED_READS},${HLA_READS},${HLA_STATUS}
EOF

# Log success status
log_status "SUCCESS" "HLA_reads=${HLA_READS}|Total_reads=${TOTAL_READS}"

echo "=========================================="
echo "✓ Sample Prepared Successfully"
echo "=========================================="
echo "Sample ID: ${SAMPLE_ID}"
echo "HLA reads: ${HLA_READS} (${HLA_STATUS})"
echo "Results saved: ${TASK_RESULT}"
echo ""
echo "Completed: $(date)"
echo "=========================================="
