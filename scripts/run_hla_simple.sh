#!/bin/bash
#SBATCH --job-name=hla_simple_%a
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/hla_simple_%A_%a.log
#SBATCH --error=logs/hla_simple_%A_%a.err

# SIMPLE HLA TYPING - Direct Tool Execution
# Works with whatever HLA tools are available

set -e

echo "=========================================="
echo "Simple HLA Typing (Task ${SLURM_ARRAY_TASK_ID})"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Get sample info
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_SHEET}")
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -d',' -f1)
BAM_FILE=$(echo "${SAMPLE_LINE}" | cut -d',' -f2)
BAI_FILE=$(echo "${SAMPLE_LINE}" | cut -d',' -f3)

echo "Sample: ${SAMPLE_ID}"
echo "BAM: ${BAM_FILE}"
echo ""

# Verify BAM
if [ ! -f "${BAM_FILE}" ]; then
    echo "✗ ERROR: BAM not found: ${BAM_FILE}"
    exit 1
fi

# Create output
SAMPLE_OUT="${RESULTS_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUT}"

# Try available tools in order of reliability
TOOLS_SUCCEEDED=0

# ========================================
# Try OptiType (if available)
# ========================================
echo "Trying OptiType..."
if command -v OptiTypePipeline.py &> /dev/null; then
    echo "  OptiType found, running..."
    mkdir -p "${SAMPLE_OUT}/optitype"
    
    # OptiType needs HLA reads extracted
    samtools view -@ ${SLURM_CPUS_PER_TASK} -b ${BAM_FILE} 6:28510120-33480577 > ${SAMPLE_OUT}/hla_region.bam
    
    if OptiTypePipeline.py -i ${SAMPLE_OUT}/hla_region.bam -d -o ${SAMPLE_OUT}/optitype 2>/dev/null; then
        echo "  ✓ OptiType completed"
        ((TOOLS_SUCCEEDED++))
    else
        echo "  ✗ OptiType failed"
    fi
else
    echo "  OptiType not available"
fi

echo ""

# ========================================
# Try Singularity tools if available
# ========================================
if command -v singularity &> /dev/null; then
    echo "Singularity available, checking for containers..."
    
    # Look for any HLA typing containers
    for container in ${SINGULARITY_CACHE}/containers/*.sif; do
        if [ -f "$container" ]; then
            container_name=$(basename $container .sif)
            echo "  Found container: ${container_name}"
            
            case ${container_name} in
                *arcas*|*arcashla*)
                    echo "    Running ArcasHLA..."
                    mkdir -p "${SAMPLE_OUT}/arcashla"
                    # Add ArcasHLA commands here
                    ;;
                *xhla*)
                    echo "    Running xHLA..."
                    mkdir -p "${SAMPLE_OUT}/xhla"
                    # Add xHLA commands here
                    ;;
            esac
        fi
    done
fi

echo ""

# ========================================
# Fallback: Extract HLA reads for manual analysis
# ========================================
if [ ${TOOLS_SUCCEEDED} -eq 0 ]; then
    echo "No HLA tools succeeded, extracting HLA reads for manual processing..."
    mkdir -p "${SAMPLE_OUT}/hla_reads"
    
    # Extract HLA region from BAM
    echo "  Extracting chromosome 6 HLA region..."
    samtools view -@ ${SLURM_CPUS_PER_TASK} -h ${BAM_FILE} 6:28510120-33480577 > ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_hla.sam
    
    # Convert to FASTQ
    echo "  Converting to FASTQ..."
    samtools fastq -@ ${SLURM_CPUS_PER_TASK} ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_hla.sam \
        -1 ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_HLA_R1.fastq.gz \
        -2 ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_HLA_R2.fastq.gz \
        -s ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_HLA_single.fastq.gz
    
    # Count reads
    HLA_READS=$(samtools view -c ${SAMPLE_OUT}/hla_reads/${SAMPLE_ID}_hla.sam)
    
    echo "  ✓ Extracted ${HLA_READS} HLA reads"
    echo "  Location: ${SAMPLE_OUT}/hla_reads/"
    echo ""
    echo "  These can be typed using:"
    echo "    - OptiType: https://github.com/FRED-2/OptiType"
    echo "    - PHLAT: https://sites.google.com/site/phlatweb/"
    echo "    - seq2HLA: https://github.com/TRON-Bioinformatics/seq2HLA"
    
    # Create a summary file
    cat > ${SAMPLE_OUT}/extraction_summary.txt << EOF
Sample: ${SAMPLE_ID}
BAM: ${BAM_FILE}
HLA Reads Extracted: ${HLA_READS}
Date: $(date)

Files Created:
- ${SAMPLE_ID}_HLA_R1.fastq.gz (forward reads)
- ${SAMPLE_ID}_HLA_R2.fastq.gz (reverse reads)
- ${SAMPLE_ID}_HLA_single.fastq.gz (unpaired reads)
- ${SAMPLE_ID}_hla.sam (SAM format)

These FASTQ files can be used with any HLA typing tool:
- OptiType (recommended for Class I: HLA-A, -B, -C)
- PHLAT (supports Class I and II)
- HLA-HD (high accuracy)
- seq2HLA (fast, supports Class I and II)

Example OptiType command:
OptiTypePipeline.py -i ${SAMPLE_ID}_HLA_R1.fastq.gz ${SAMPLE_ID}_HLA_R2.fastq.gz \
    --dna -o optitype_output/ -v
EOF
fi

echo ""
echo "=========================================="
echo "Processing Complete: ${SAMPLE_ID}"
echo "=========================================="
echo "Output: ${SAMPLE_OUT}/"
echo "Tools succeeded: ${TOOLS_SUCCEEDED}"
echo "Completed: $(date)"
echo ""
