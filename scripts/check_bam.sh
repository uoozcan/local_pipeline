#!/bin/bash
# Quick validation of transferred BAM files

WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

echo "=========================================="
echo "BAM Files Quick Check"
echo "=========================================="
echo ""

# Find all BAM files
bam_files=$(find ${RAW_BAM_DIR} -name "*_tumor.bam" -type f | sort)

if [ -z "$bam_files" ]; then
    echo "No BAM files found in ${RAW_BAM_DIR}"
    exit 1
fi

total=$(echo "$bam_files" | wc -l)
with_index=0
without_index=0

echo "Found ${total} BAM files"
echo ""

for bam in $bam_files; do
    sample=$(basename $(dirname $bam))
    size=$(du -h $bam | cut -f1)
    
    if [ -f "${bam}.bai" ]; then
        status="✓ indexed"
        ((with_index++))
    else
        status="✗ no index"
        ((without_index++))
    fi
    
    echo "${sample}: ${size} ${status}"
done

echo ""
echo "Summary:"
echo "  With index: ${with_index}"
echo "  Without index: ${without_index}"
echo ""

if [ ${without_index} -gt 0 ]; then
    echo "To create missing indices, run:"
    echo "  bash ${SCRIPTS_DIR}/step2_submit_array_job.sh"
fi
