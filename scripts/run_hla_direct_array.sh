#!/bin/bash
#SBATCH --job-name=hla_direct_%a
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/hla_direct_%A_%a.log
#SBATCH --error=logs/hla_direct_%A_%a.err

# DIRECT HLA TYPING - Bypass Buggy Pipeline
# Runs HLA typing tools directly on BAM files using array job

set -e

echo "=========================================="
echo "Direct HLA Typing (Array Task)"
echo "=========================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Load required modules
echo "Loading modules..."
module load biokit
module load samtools
module load python-data

echo "✓ Modules loaded"
echo ""

# Set up Singularity
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"
export SINGULARITY_BIND="/scratch/${PROJECT_ID}"

# Get sample info from samplesheet
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_SHEET}")
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -d',' -f1)
BAM_FILE=$(echo "${SAMPLE_LINE}" | cut -d',' -f2)
BAI_FILE=$(echo "${SAMPLE_LINE}" | cut -d',' -f3)

echo "Processing sample: ${SAMPLE_ID}"
echo "  BAM: ${BAM_FILE}"
echo "  BAI: ${BAI_FILE}"
echo ""

# Create output directories
SAMPLE_OUTPUT="${RESULTS_DIR}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUTPUT}"/{arcashla,xhla,hlahd}

# Verify BAM exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "✗ ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "✓ BAM file found ($(du -h ${BAM_FILE} | cut -f1))"
echo ""

# ========================================
# Run ArcasHLA
# ========================================
echo "=========================================="
echo "Running ArcasHLA"
echo "=========================================="
echo ""

ARCAS_OUTPUT="${SAMPLE_OUTPUT}/arcashla"
cd "${ARCAS_OUTPUT}"

# ArcasHLA container
ARCAS_SIF="${SINGULARITY_CACHE}/containers/arcashla.sif"

if [ -f "${ARCAS_SIF}" ]; then
    echo "Using ArcasHLA container..."
    
    # Extract HLA reads from BAM
    singularity exec ${ARCAS_SIF} arcasHLA extract \
        ${BAM_FILE} \
        -o . \
        -t ${SLURM_CPUS_PER_TASK} \
        --log arcashla_extract.log || echo "ArcasHLA extract failed"
    
    # Run genotyping
    if [ -f "${SAMPLE_ID}.extracted.fq.gz" ]; then
        singularity exec ${ARCAS_SIF} arcasHLA genotype \
            ${SAMPLE_ID}.extracted.fq.gz \
            -o . \
            -t ${SLURM_CPUS_PER_TASK} \
            --log arcashla_genotype.log || echo "ArcasHLA genotype failed"
        
        echo "✓ ArcasHLA completed"
    else
        echo "✗ ArcasHLA extraction failed"
    fi
else
    echo "⚠ ArcasHLA container not found, skipping"
fi

echo ""

# ========================================
# Run xHLA
# ========================================
echo "=========================================="
echo "Running xHLA"
echo "=========================================="
echo ""

XHLA_OUTPUT="${SAMPLE_OUTPUT}/xhla"
cd "${XHLA_OUTPUT}"

XHLA_SIF="${SINGULARITY_CACHE}/containers/xhla.sif"

if [ -f "${XHLA_SIF}" ]; then
    echo "Using xHLA container..."
    
    # xHLA requires JSON config
    cat > xhla_config.json << EOF
{
    "sample_id": "${SAMPLE_ID}",
    "bam": "${BAM_FILE}"
}
EOF
    
    singularity exec ${XHLA_SIF} \
        run.py \
        --sample_id ${SAMPLE_ID} \
        --input_bam_path ${BAM_FILE} \
        --output_path . \
        --num_threads ${SLURM_CPUS_PER_TASK} || echo "xHLA failed"
    
    echo "✓ xHLA completed"
else
    echo "⚠ xHLA container not found, skipping"
fi

echo ""

# ========================================
# Run HLA-HD
# ========================================
echo "=========================================="
echo "Running HLA-HD"
echo "=========================================="
echo ""

HLAHD_OUTPUT="${SAMPLE_OUTPUT}/hlahd"
cd "${HLAHD_OUTPUT}"

HLAHD_SIF="${SINGULARITY_CACHE}/containers/hlahd.sif"

if [ -f "${HLAHD_SIF}" ]; then
    echo "Using HLA-HD container..."
    
    # Convert BAM to FASTQ for HLA-HD
    echo "Converting BAM to FASTQ..."
    samtools fastq -@ ${SLURM_CPUS_PER_TASK} ${BAM_FILE} \
        -1 ${SAMPLE_ID}_R1.fastq.gz \
        -2 ${SAMPLE_ID}_R2.fastq.gz \
        -s ${SAMPLE_ID}_single.fastq.gz || echo "BAM to FASTQ conversion failed"
    
    if [ -f "${SAMPLE_ID}_R1.fastq.gz" ]; then
        singularity exec ${HLAHD_SIF} \
            hlahd.sh \
            -t ${SLURM_CPUS_PER_TASK} \
            -f ${HLA_REFERENCES}/hlahd_freq_data \
            ${SAMPLE_ID}_R1.fastq.gz \
            ${SAMPLE_ID}_R2.fastq.gz \
            ${HLA_REFERENCES}/hlahd_dictionary \
            ${SAMPLE_ID} \
            . || echo "HLA-HD failed"
        
        echo "✓ HLA-HD completed"
    else
        echo "✗ FASTQ conversion failed, skipping HLA-HD"
    fi
else
    echo "⚠ HLA-HD container not found, skipping"
fi

echo ""

# ========================================
# Summary
# ========================================
echo "=========================================="
echo "Sample Processing Complete"
echo "=========================================="
echo ""
echo "Sample: ${SAMPLE_ID}"
echo "Results: ${SAMPLE_OUTPUT}/"
echo ""

# Check what succeeded
SUCCESS_TOOLS=""
if [ -f "${SAMPLE_OUTPUT}/arcashla/${SAMPLE_ID}.genotype.json" ]; then
    echo "✓ ArcasHLA results found"
    SUCCESS_TOOLS="${SUCCESS_TOOLS},arcashla"
fi

if [ -f "${SAMPLE_OUTPUT}/xhla/report-${SAMPLE_ID}.json" ]; then
    echo "✓ xHLA results found"
    SUCCESS_TOOLS="${SUCCESS_TOOLS},xhla"
fi

if [ -f "${SAMPLE_OUTPUT}/hlahd/${SAMPLE_ID}/result/${SAMPLE_ID}_final.result.txt" ]; then
    echo "✓ HLA-HD results found"
    SUCCESS_TOOLS="${SUCCESS_TOOLS},hlahd"
fi

echo ""
echo "Successful tools: ${SUCCESS_TOOLS}"
echo "Completed: $(date)"
echo ""
