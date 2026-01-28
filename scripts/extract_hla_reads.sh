#!/bin/bash
#SBATCH --job-name=extract_hla_%a
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/extract_hla_%A_%a.log
#SBATCH --error=logs/extract_hla_%A_%a.err

# SIMPLEST SOLUTION: Extract HLA Reads from BAM Files
# This ALWAYS works - no dependencies on buggy pipelines

set -e

echo "=========================================="
echo "HLA Read Extraction"
echo "=========================================="
echo "Task: ${SLURM_ARRAY_TASK_ID}/24"
echo "Started: $(date)"
echo ""

# Configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
SAMPLESHEET="${WORK_DIR}/pipeline_input_bam/batch1_VenEx_DNA_BAM/samplesheet_pipeline.csv"
RESULTS="${WORK_DIR}/results_bam/batch1_VenEx_DNA_BAM"

# Load samtools
module load biokit

# Get sample info
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLESHEET}")
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -d',' -f1)
BAM_FILE=$(echo "${SAMPLE_LINE}" | cut -d',' -f2)

echo "Sample: ${SAMPLE_ID}"
echo "BAM: ${BAM_FILE}"
echo ""

# Check BAM exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "✗ ERROR: BAM file not found"
    exit 1
fi

echo "✓ BAM found ($(du -h ${BAM_FILE} | cut -f1))"
echo ""

# Create output directory
OUTPUT="${RESULTS}/${SAMPLE_ID}/hla_reads"
mkdir -p "${OUTPUT}"

# Extract HLA region from chromosome 6
echo "Extracting HLA region (chr6:28510120-33480577)..."

# Try with "6" first (most common)
if samtools view -H "${BAM_FILE}" | grep -q "SN:6"; then
    CHR="6"
    echo "  Using chromosome: 6"
elif samtools view -H "${BAM_FILE}" | grep -q "SN:chr6"; then
    CHR="chr6"
    echo "  Using chromosome: chr6"
else
    echo "✗ ERROR: Could not find chromosome 6 in BAM header"
    exit 1
fi

# Extract HLA reads
samtools view -@ ${SLURM_CPUS_PER_TASK} -h "${BAM_FILE}" ${CHR}:28510120-33480577 \
    > "${OUTPUT}/${SAMPLE_ID}_hla.sam"

# Check if we got reads
HLA_READS=$(samtools view -c "${OUTPUT}/${SAMPLE_ID}_hla.sam")

if [ ${HLA_READS} -eq 0 ]; then
    echo "✗ WARNING: No HLA reads extracted!"
    echo "  BAM may not have coverage in HLA region"
else
    echo "✓ Extracted ${HLA_READS} reads"
fi

# Convert to FASTQ for use with HLA typing tools
echo ""
echo "Converting to FASTQ format..."

samtools fastq -@ ${SLURM_CPUS_PER_TASK} "${OUTPUT}/${SAMPLE_ID}_hla.sam" \
    -1 "${OUTPUT}/${SAMPLE_ID}_R1.fastq.gz" \
    -2 "${OUTPUT}/${SAMPLE_ID}_R2.fastq.gz" \
    -s "${OUTPUT}/${SAMPLE_ID}_single.fastq.gz" \
    -0 /dev/null -n

echo "✓ FASTQ files created"
echo ""

# Create README
cat > "${OUTPUT}/README.txt" << EOF
HLA Reads Extraction Summary
============================

Sample: ${SAMPLE_ID}
Source BAM: ${BAM_FILE}
Extraction Date: $(date)
HLA Reads: ${HLA_READS}

Files Created:
--------------
${SAMPLE_ID}_R1.fastq.gz      - Forward reads (paired)
${SAMPLE_ID}_R2.fastq.gz      - Reverse reads (paired)
${SAMPLE_ID}_single.fastq.gz  - Unpaired reads
${SAMPLE_ID}_hla.sam          - SAM format

Use These Files With Any HLA Typing Tool:
-----------------------------------------

OptiType (HLA-A, B, C):
  OptiTypePipeline.py -i ${SAMPLE_ID}_R1.fastq.gz ${SAMPLE_ID}_R2.fastq.gz --dna -o output/

PHLAT (HLA-A, B, C, DRB1, DQB1, DPB1):
  python PHLAT.py -1 ${SAMPLE_ID}_R1.fastq.gz -2 ${SAMPLE_ID}_R2.fastq.gz -o output/

HLA-HD (High resolution, all HLA genes):
  hlahd.sh -t 4 -f freq_data/ ${SAMPLE_ID}_R1.fastq.gz ${SAMPLE_ID}_R2.fastq.gz dictionary/ ${SAMPLE_ID} output/

seq2HLA (HLA-A, B, C, DRB1, DQB1):
  seq2HLA -1 ${SAMPLE_ID}_R1.fastq.gz -2 ${SAMPLE_ID}_R2.fastq.gz -r output/${SAMPLE_ID}

Online Tools:
  - IMGTHLA (https://www.ebi.ac.uk/ipd/imgt/hla/)
  - Galaxy (usegalaxy.org) - has HLA typing workflows

EOF

# File sizes
echo "File sizes:"
ls -lh "${OUTPUT}/" | grep -E "\.fastq\.gz|\.sam" | awk '{print "  " $9 ": " $5}'
echo ""

# Summary
echo "=========================================="
echo "✓ Extraction Complete"
echo "=========================================="
echo ""
echo "Sample: ${SAMPLE_ID}"
echo "HLA reads: ${HLA_READS}"
echo "Output: ${OUTPUT}/"
echo ""
echo "These FASTQ files can now be analyzed with:"
echo "  - OptiType (recommended)"
echo "  - PHLAT"
echo "  - HLA-HD"
echo "  - seq2HLA"
echo "  - Any other HLA typing tool"
echo ""
echo "See README.txt in output directory for commands."
echo ""
echo "Completed: $(date)"
