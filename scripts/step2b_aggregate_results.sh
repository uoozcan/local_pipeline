#!/bin/bash
#SBATCH --job-name=aggregate_bam_prep
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/aggregate_bam_prep_%j.log
#SBATCH --error=logs/aggregate_bam_prep_%j.err

# STEP 2b: Aggregate BAM Preparation Results
# Purpose: Combine results from array job into final samplesheet
# Run this after step2_prepare_bam_array.sh completes

set -e

echo "=========================================="
echo "Aggregating BAM Preparation Results"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

ARRAY_RESULTS_DIR="${INPUT_DIR}/array_results"

# Check if array results exist
if [ ! -d "${ARRAY_RESULTS_DIR}" ]; then
    echo "ERROR: Array results directory not found: ${ARRAY_RESULTS_DIR}"
    echo "Please run step2_prepare_bam_array.sh first"
    exit 1
fi

# Count results
TOTAL_TASKS=$(ls ${ARRAY_RESULTS_DIR}/task_*.txt 2>/dev/null | wc -l)
SUCCESS_TASKS=$(grep -l "SUCCESS" ${ARRAY_RESULTS_DIR}/status_*.txt 2>/dev/null | wc -l || echo "0")
FAILED_TASKS=$((TOTAL_TASKS - SUCCESS_TASKS))

echo "Array Job Results:"
echo "  Total tasks completed: ${TOTAL_TASKS}"
echo "  Successful: ${SUCCESS_TASKS}"
echo "  Failed: ${FAILED_TASKS}"
echo ""

if [ ${TOTAL_TASKS} -eq 0 ]; then
    echo "ERROR: No array results found"
    exit 1
fi

# Create detailed samplesheet header (with statistics)
echo "sample,bam,bai,total_reads,mapped_reads,hla_reads,hla_status" > "${DETAILED_SHEET}"

# Aggregate successful results
echo "Aggregating results..."
for result_file in ${ARRAY_RESULTS_DIR}/task_*.txt; do
    if [ -f "${result_file}" ]; then
        # Skip header and append data
        tail -n +2 "${result_file}" >> "${DETAILED_SHEET}"
    fi
done

# Sort by sample name
TEMP_SHEET="${DETAILED_SHEET}.tmp"
head -n 1 "${DETAILED_SHEET}" > "${TEMP_SHEET}"
tail -n +2 "${DETAILED_SHEET}" | sort >> "${TEMP_SHEET}"
mv "${TEMP_SHEET}" "${DETAILED_SHEET}"

VALID_SAMPLES=$(tail -n +2 "${DETAILED_SHEET}" | wc -l)

echo "✓ Detailed samplesheet created: ${DETAILED_SHEET}"
echo "  Valid samples: ${VALID_SAMPLES}"
echo ""

# Create simplified samplesheet for pipeline (without statistics)
echo "sample,bam,bai" > "${SAMPLE_SHEET}"
tail -n +2 "${DETAILED_SHEET}" | cut -d',' -f1,2,3 >> "${SAMPLE_SHEET}"

echo "✓ Pipeline samplesheet created: ${SAMPLE_SHEET}"
echo ""

# Create detailed summary
SUMMARY_FILE="${INPUT_DIR}/preparation_summary.txt"
cat > ${SUMMARY_FILE} << EOF
BAM Preparation Summary (Array Job)
====================================
Date: $(date)
Batch: ${BATCH_NAME}

Array Job Statistics:
--------------------
Total tasks: ${TOTAL_TASKS}
Successful tasks: ${SUCCESS_TASKS}
Failed tasks: ${FAILED_TASKS}
Valid samples: ${VALID_SAMPLES}

Sample Coverage Analysis:
------------------------
EOF

# Analyze HLA coverage distribution
echo "" >> ${SUMMARY_FILE}
echo "HLA Coverage Summary:" >> ${SUMMARY_FILE}
GOOD_COVERAGE=$(tail -n +2 "${DETAILED_SHEET}" | grep -c ",GOOD" || echo "0")
LOW_COVERAGE=$(tail -n +2 "${DETAILED_SHEET}" | grep -c ",LOW" || echo "0")
VERY_LOW_COVERAGE=$(tail -n +2 "${DETAILED_SHEET}" | grep -c ",VERY_LOW" || echo "0")

echo "  GOOD coverage (>1000 reads): ${GOOD_COVERAGE}" >> ${SUMMARY_FILE}
echo "  LOW coverage (100-1000 reads): ${LOW_COVERAGE}" >> ${SUMMARY_FILE}
echo "  VERY LOW coverage (<100 reads): ${VERY_LOW_COVERAGE}" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# List failed samples if any
if [ ${FAILED_TASKS} -gt 0 ]; then
    echo "Failed Samples:" >> ${SUMMARY_FILE}
    echo "--------------" >> ${SUMMARY_FILE}
    for status_file in ${ARRAY_RESULTS_DIR}/status_*.txt; do
        if [ -f "${status_file}" ]; then
            if grep -q "FAILED" "${status_file}"; then
                SAMPLE=$(cut -d'|' -f2 "${status_file}")
                REASON=$(cut -d'|' -f4 "${status_file}")
                echo "  ${SAMPLE}: ${REASON}" >> ${SUMMARY_FILE}
            fi
        fi
    done
    echo "" >> ${SUMMARY_FILE}
fi

# List all valid samples
echo "Valid Sample List:" >> ${SUMMARY_FILE}
echo "-----------------" >> ${SUMMARY_FILE}
tail -n +2 "${SAMPLE_SHEET}" | cut -d',' -f1 | sort >> ${SUMMARY_FILE}

echo "✓ Summary saved: ${SUMMARY_FILE}"
echo ""

# Display summary
echo "=========================================="
echo "Preparation Summary"
echo "=========================================="
echo ""
echo "Valid samples: ${VALID_SAMPLES}"
echo ""
echo "HLA Coverage Distribution:"
echo "  GOOD (>1000 reads): ${GOOD_COVERAGE}"
echo "  LOW (100-1000 reads): ${LOW_COVERAGE}"
echo "  VERY LOW (<100 reads): ${VERY_LOW_COVERAGE}"
echo ""

if [ ${FAILED_TASKS} -gt 0 ]; then
    echo "⚠ WARNING: ${FAILED_TASKS} samples failed preparation"
    echo "  See ${SUMMARY_FILE} for details"
    echo ""
fi

# Check if we have enough valid samples
if [ ${VALID_SAMPLES} -eq 0 ]; then
    echo "ERROR: No valid samples prepared!"
    echo "Cannot proceed with pipeline"
    exit 1
fi

# Preview samplesheet
echo "Sample Sheet Preview (first 10 samples):"
echo "----------------------------------------"
head -11 "${SAMPLE_SHEET}" | column -t -s','
echo ""

if [ ${VALID_SAMPLES} -gt 10 ]; then
    echo "... and $((VALID_SAMPLES - 10)) more samples"
    echo ""
fi

echo "=========================================="
echo "Aggregation Complete!"
echo "=========================================="
echo ""
echo "Files Created:"
echo "  ${DETAILED_SHEET} (detailed with statistics)"
echo "  ${SAMPLE_SHEET} (simplified for pipeline)"
echo "  ${SUMMARY_FILE} (summary report)"
echo ""
echo "Completed: $(date)"
echo ""
echo "NEXT STEP:"
echo "  sbatch ${SCRIPTS_DIR}/step3_run_pipeline_bam.sh"
echo ""
