#!/bin/bash
#SBATCH --job-name=prepare_batch_bam
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/prepare_batch_bam_%j.log
#SBATCH --error=logs/prepare_batch_bam_%j.err

# STEP 2: Prepare Batch BAM Files for Pipeline
# Purpose: Validate BAM files, create indices, generate samplesheet

set -e

echo "=========================================="
echo "Preparing Batch BAM Files for HLA Typing"
echo "=========================================="
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

# Create input directory
mkdir -p "${INPUT_DIR}"

# Initialize samplesheet
echo "sample,bam,bai" > "${SAMPLE_SHEET}"

# Counters
total_samples=0
valid_samples=0
invalid_samples=0

# Function to validate and prepare single BAM file
prepare_bam_sample() {
    local bam_path=$1
    local sample_id=$(basename $(dirname ${bam_path}))
    
    echo "----------------------------------------"
    echo "Processing: ${sample_id}"
    echo "----------------------------------------"
    
    # Check if BAM exists
    if [ ! -f "${bam_path}" ]; then
        echo "  ✗ BAM file not found: ${bam_path}"
        ((invalid_samples++))
        return 1
    fi
    
    echo "  ✓ BAM file found"
    echo "    Size: $(du -h ${bam_path} | cut -f1)"
    
    # Check/create BAM index
    bai_path="${bam_path}.bai"
    if [ ! -f "${bai_path}" ]; then
        echo "  Creating BAM index..."
        if samtools index "${bam_path}"; then
            echo "    ✓ Index created"
        else
            echo "    ✗ Index creation failed"
            ((invalid_samples++))
            return 1
        fi
    else
        echo "  ✓ Index found"
    fi
    
    # Quick validation
    echo "  Validating BAM structure..."
    if samtools quickcheck "${bam_path}"; then
        echo "    ✓ BAM file is valid"
    else
        echo "    ✗ BAM file is corrupted"
        ((invalid_samples++))
        return 1
    fi
    
    # Check chromosome naming
    echo "  Detecting chromosome naming..."
    chr_with_prefix=$(samtools view -H "${bam_path}" | grep -c "^@SQ.*SN:chr" || echo "0")
    
    if [ ${chr_with_prefix} -gt 0 ]; then
        chr6_name="chr6"
        echo "    Chromosomes: WITH 'chr' prefix"
    else
        chr6_name="6"
        echo "    Chromosomes: WITHOUT 'chr' prefix"
    fi
    
    # Check HLA region coverage
    echo "  Checking HLA region (${chr6_name}:28510120-33480577)..."
    hla_reads=$(samtools view -c "${bam_path}" ${chr6_name}:28510120-33480577 2>/dev/null || echo "0")
    
    if [ "${hla_reads}" -gt 1000 ]; then
        echo "    ✓ Found ${hla_reads} reads in HLA region"
    else
        echo "    ⚠ WARNING: Only ${hla_reads} reads in HLA region"
        echo "      Pipeline may produce low-confidence results"
    fi
    
    # Get basic statistics
    echo "  BAM statistics:"
    total_reads=$(samtools view -c "${bam_path}" 2>/dev/null || echo "0")
    mapped_reads=$(samtools view -c -F 4 "${bam_path}" 2>/dev/null || echo "0")
    echo "    Total reads: ${total_reads}"
    echo "    Mapped reads: ${mapped_reads}"
    
    if [ ${total_reads} -gt 0 ]; then
        hla_percent=$(echo "scale=2; ${hla_reads} * 100 / ${total_reads}" | bc)
        echo "    HLA reads: ${hla_percent}% of total"
    fi
    
    # Add to samplesheet
    echo "${sample_id},${bam_path},${bai_path}" >> "${SAMPLE_SHEET}"
    
    echo "  ✓ Sample prepared successfully"
    ((valid_samples++))
    echo ""
    
    return 0
}

# Find all BAM files
echo "Searching for BAM files in: ${RAW_BAM_DIR}"
echo ""

bam_files=$(find ${RAW_BAM_DIR} -name "*_tumor.bam" -type f | sort)

if [ -z "$bam_files" ]; then
    echo "ERROR: No BAM files found in ${RAW_BAM_DIR}"
    echo ""
    echo "Please run step1_transfer_bam.sh first"
    exit 1
fi

total_samples=$(echo "$bam_files" | wc -l)
echo "Found ${total_samples} BAM files"
echo ""

# Process each BAM file
for bam_file in $bam_files; do
    prepare_bam_sample "$bam_file" || true
done

# Summary
echo "=========================================="
echo "Preparation Summary"
echo "=========================================="
echo ""
echo "Total samples found: ${total_samples}"
echo "Valid samples: ${valid_samples}"
echo "Invalid samples: ${invalid_samples}"
echo ""
echo "Sample sheet created: ${SAMPLE_SHEET}"
echo "Sample sheet entries: $(tail -n +2 ${SAMPLE_SHEET} | wc -l)"
echo ""

if [ $valid_samples -eq 0 ]; then
    echo "ERROR: No valid samples prepared!"
    echo "Cannot proceed with pipeline"
    exit 1
fi

# Show sample sheet preview
echo "Sample sheet preview (first 10 samples):"
echo "----------------------------------------"
head -11 "${SAMPLE_SHEET}"
echo ""

# Create summary file
summary_file="${INPUT_DIR}/preparation_summary.txt"
cat > ${summary_file} << EOF
BAM Preparation Summary
======================
Date: $(date)
Batch: ${BATCH_NAME}

Statistics:
-----------
Total samples: ${total_samples}
Valid samples: ${valid_samples}
Invalid samples: ${invalid_samples}

Input Directory: ${INPUT_DIR}
Sample Sheet: ${SAMPLE_SHEET}

Sample List:
-----------
EOF

tail -n +2 "${SAMPLE_SHEET}" | cut -d',' -f1 >> ${summary_file}

echo "✓ Summary saved: ${summary_file}"
echo ""

echo "=========================================="
echo "Preparation Complete!"
echo "=========================================="
echo ""
echo "Completed: $(date)"
echo ""
echo "NEXT STEP:"
echo "  sbatch step3_run_batch_pipeline.sh"
echo ""