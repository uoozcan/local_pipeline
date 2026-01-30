#!/bin/bash
#SBATCH --job-name=hla_pipeline
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/hla_pipeline_%j.out
#SBATCH --error=logs/hla_pipeline_%j.err

################################################################################
# HLA Typing Pipeline - SLURM Submission Script
# For CSC Puhti HPC Environment
#
# This script submits the HLA typing pipeline to SLURM with proper
# parameter handling. It fixes the issue where parameters were being
# executed as separate commands.
################################################################################

set -euo pipefail

# Create logs directory
mkdir -p logs

# Print job information
echo "=========================================="
echo "HLA Typing Pipeline - RNA-seq"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""

# Configuration - MODIFY THESE PARAMETERS
INPUT_DIR="/scratch/project_2008084/your_samples"
INPUT_TYPE="fastq"                    # bam, cram, or fastq
TOOLS="optitype,arcashla"             # Comma-separated: optitype,arcashla,spechla
OUTDIR="./results"
REFERENCE="hg38"                      # hg19 or hg38

# Tool-specific options
OPTITYPE_SEQ_TYPE="rna"               # dna or rna
ARCASHLA_GENES="A,B,C,DQA1,DQB1,DRB1"
SPECHLA_GENES="A,B,C,DQA1,DQB1,DRB1"
SPECHLA_EXON_ONLY=0                   # Set to 1 for exome data!

# Majority voting options
ENABLE_MAJORITY_VOTING="false"        # true or false
MV_MIN_TOOLS=2
MV_RESOLUTION=2
MV_GENES="A,B,C,DQA1,DQB1,DRB1"

# HPC resources
SLURM_ACCOUNT="project_2008084"
SLURM_PARTITION="small"
MAX_CPUS=40
MAX_MEMORY="180.GB"
MAX_TIME="24.h"

################################################################################
# DO NOT MODIFY BELOW THIS LINE (unless you know what you're doing)
################################################################################

# Load required modules
module purge
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0

echo "Modules loaded:"
module list
echo ""

# Count samples
if [ -d "$INPUT_DIR" ]; then
    if [ "$INPUT_TYPE" = "fastq" ]; then
        SAMPLE_COUNT=$(ls -1 "$INPUT_DIR"/*_R1*.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l || echo "0")
    else
        SAMPLE_COUNT=$(ls -1 "$INPUT_DIR"/*.{bam,cram} 2>/dev/null | wc -l || echo "0")
    fi
    echo "Processing $SAMPLE_COUNT samples"
    echo ""
fi

# Show samples (first 10)
if [ -d "$INPUT_DIR" ]; then
    echo "Samples:"
    if [ "$INPUT_TYPE" = "fastq" ]; then
        ls -1 "$INPUT_DIR"/*_R1*.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | head -10 | nl || true
    else
        ls -1 "$INPUT_DIR"/*.{bam,cram} 2>/dev/null | head -10 | nl || true
    fi
    echo ""
fi

# Run Nextflow pipeline
# CRITICAL: All parameters MUST be on the same line or properly escaped
# with backslashes. This is the fix for the "command not found" errors.

echo "Starting Nextflow pipeline..."
echo ""

nextflow run main.nf \
    --input "$INPUT_DIR" \
    --input_type "$INPUT_TYPE" \
    --tools "$TOOLS" \
    --reference_build "$REFERENCE" \
    --optitype_seq_type "$OPTITYPE_SEQ_TYPE" \
    --arcashla_genes "$ARCASHLA_GENES" \
    --spechla_genes "$SPECHLA_GENES" \
    --spechla_exon_only "$SPECHLA_EXON_ONLY" \
    --enable_majority_voting "$ENABLE_MAJORITY_VOTING" \
    --mv_min_tools "$MV_MIN_TOOLS" \
    --mv_resolution "$MV_RESOLUTION" \
    --mv_genes "$MV_GENES" \
    --slurm_account "$SLURM_ACCOUNT" \
    --slurm_partition "$SLURM_PARTITION" \
    --max_cpus "$MAX_CPUS" \
    --max_memory "$MAX_MEMORY" \
    --max_time "$MAX_TIME" \
    --outdir "$OUTDIR" \
    -profile puhti,singularity \
    -with-report "$OUTDIR/pipeline_info/report_${SLURM_JOB_ID}.html" \
    -with-timeline "$OUTDIR/pipeline_info/timeline_${SLURM_JOB_ID}.html" \
    -with-dag "$OUTDIR/pipeline_info/dag_${SLURM_JOB_ID}.html" \
    -resume

EXIT_CODE=$?

# Print completion status
echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline failed with exit code: $EXIT_CODE"
    echo "Check logs:"
    echo "  - logs/hla_pipeline_${SLURM_JOB_ID}.out"
    echo "  - logs/hla_pipeline_${SLURM_JOB_ID}.err"
    echo "  - .nextflow.log"
fi
echo "Completed: $(date)"
echo "=========================================="

exit $EXIT_CODE
