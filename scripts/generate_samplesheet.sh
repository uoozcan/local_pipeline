#!/bin/bash

# Generate Nextflow Samplesheet for Batch 1
BATCH_NAME="batch1_VenEx_DNA"
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
FASTQ_DIR="${WORK_DIR}/raw_fastq/${BATCH_NAME}"
SAMPLESHEET="${WORK_DIR}/samplesheet_${BATCH_NAME}.csv"

echo "Generating Nextflow samplesheet..."
echo "Output: ${SAMPLESHEET}"

# Create header
echo "sample,fastq_1,fastq_2,strandedness" > "${SAMPLESHEET}"

# Find all R1 files and create entries
cd "${FASTQ_DIR}"
find . -name "*_R1_*.fastq.gz" | sort | while read r1; do
    # Get R2 file
    r2="${r1/_R1_/_R2_}"
    
    # Check R2 exists
    if [ ! -f "$r2" ]; then
        echo "WARNING: Missing R2 for $r1" >&2
        continue
    fi
    
    # Extract sample name
    # For files like: VenEx_DNA_NONHUS_Batch1/VX_100_2_D1/VX_100_2_D1_S38_R1_001.fastq.gz
    # Extract: VX_100_2_D1
    sample=$(basename "$r1" | sed -E 's/(_S[0-9]+)?(_L[0-9]+)?_R[12]_001\.fastq\.gz//')
    
    # Full paths
    r1_full="${FASTQ_DIR}/${r1#./}"
    r2_full="${FASTQ_DIR}/${r2#./}"
    
    # Add to samplesheet (DNA-seq is unstranded)
    echo "${sample},${r1_full},${r2_full},unstranded" >> "${SAMPLESHEET}"
done

# Show summary
sample_count=$(tail -n +2 "${SAMPLESHEET}" | wc -l)
echo ""
echo "Samplesheet created: ${SAMPLESHEET}"
echo "Total samples: ${sample_count}"
echo ""
echo "First 5 entries:"
head -6 "${SAMPLESHEET}"