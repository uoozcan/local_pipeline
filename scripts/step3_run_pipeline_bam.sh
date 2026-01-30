#!/bin/bash
#SBATCH --job-name=hla_batch_bam
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --output=logs/hla_batch_bam_%j.log
#SBATCH --error=logs/hla_batch_bam_%j.err

# STEP 3: Run HLA Typing Pipeline on Batch BAM Files
# Tools: OptiType, ArcasHLA, SpecHLA, Majority Voting

set -e

echo "=========================================="
echo "HLA Typing Pipeline - Batch BAM Analysis"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Load required modules
echo "Loading modules..."
module load biokit
module load nextflow

echo "✓ Modules loaded"
echo ""

# Verify pipeline exists
if [ ! -f "${PIPELINE_DIR}/main.nf" ]; then
    echo "ERROR: Pipeline not found at ${PIPELINE_DIR}/main.nf"
    exit 1
fi

echo "✓ Pipeline found: ${PIPELINE_DIR}/main.nf"
echo ""

# Verify samplesheet exists
if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet not found: ${SAMPLE_SHEET}"
    echo "Please run step2_submit_array_job.sh first"
    exit 1
fi

sample_count=$(tail -n +2 ${SAMPLE_SHEET} | wc -l)
echo "✓ Sample sheet found: ${SAMPLE_SHEET}"
echo "  Samples to process: ${sample_count}"
echo ""

# Set up Nextflow environment
export NXF_OPTS="-Xms1G -Xmx6G"
export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"

# Create work directory
WORK_NF_DIR="${BASE_DIR}/work_batch_bam"
mkdir -p ${WORK_NF_DIR}

echo "=========================================="
echo "Pipeline Configuration"
echo "=========================================="
echo "Input type:       BAM"
echo "Batch:            ${BATCH_NAME}"
echo "Samples:          ${sample_count}"
echo "Sample sheet:     ${SAMPLE_SHEET}"
echo "Output directory: ${RESULTS_DIR}"
echo "Work directory:   ${WORK_NF_DIR}"
echo ""
echo "Tools:            ${TOOLS}"
echo "HLA genes:        ${HLA_GENES}"
echo "Seq type:         ${SEQ_TYPE}"
echo "Consensus:        ${ENABLE_CONSENSUS}"
echo "Min tools:        ${MIN_TOOLS_CONSENSUS}"
echo ""
echo "Resources:"
echo "  Max CPUs:       ${MAX_CPUS}"
echo "  Max Memory:     ${MAX_MEMORY}"
echo "  Max Time:       ${MAX_TIME}"
echo ""

# Pipeline command
PIPELINE_CMD="nextflow run ${PIPELINE_DIR}/main.nf \
    --input ${SAMPLE_SHEET} \
    --input_type csv \
    --tools ${TOOLS} \
    --hla_genes '${HLA_GENES}' \
    --seq_type ${SEQ_TYPE} \
    --enable_majority_voting \
    --mv_min_tools ${MIN_TOOLS_CONSENSUS} \
    --mv_resolution ${CONSENSUS_RESOLUTION} \
    --mv_genes '${HLA_GENES}' \
    --reference_genome ${REFERENCE_GENOME} \
    --slurm_account ${SLURM_ACCOUNT} \
    --slurm_partition ${SLURM_PARTITION} \
    --max_cpus ${MAX_CPUS} \
    --max_memory ${MAX_MEMORY} \
    --max_time ${MAX_TIME} \
    --outdir ${RESULTS_DIR} \
    -work-dir ${WORK_NF_DIR} \
    -profile singularity \
    -resume \
    -with-report ${RESULTS_DIR}/pipeline_report.html \
    -with-timeline ${RESULTS_DIR}/pipeline_timeline.html \
    -with-dag ${RESULTS_DIR}/pipeline_dag.svg"

echo "=========================================="
echo "Starting Pipeline Execution"
echo "=========================================="
echo ""
echo "Full command:"
echo "${PIPELINE_CMD}" | sed 's/ --/\n  --/g'
echo ""
echo "This may take several hours for ${sample_count} samples..."
echo "----------------------------------------"
echo ""

# Run pipeline
${PIPELINE_CMD}

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ Pipeline Completed Successfully!"
    echo "=========================================="
    echo ""
    echo "Results saved in: ${RESULTS_DIR}"
    echo ""
    echo "Summary files:"
    echo "  - Pipeline report:  ${RESULTS_DIR}/pipeline_report.html"
    echo "  - Timeline:         ${RESULTS_DIR}/pipeline_timeline.html"
    echo "  - DAG:              ${RESULTS_DIR}/pipeline_dag.svg"
    echo ""
    
    if [ -f "${RESULTS_DIR}/majority_voting/all_samples.consensus.tsv" ]; then
        echo "✓ Consensus results for all samples:"
        echo "  ${RESULTS_DIR}/majority_voting/all_samples.consensus.tsv"
        echo ""
        echo "Preview (first 10 samples):"
        head -11 "${RESULTS_DIR}/majority_voting/all_samples.consensus.tsv" || true
        echo ""
    fi
    
    echo "Individual sample results:"
    echo "  ${RESULTS_DIR}/"
    ls -d ${RESULTS_DIR}/*/ 2>/dev/null | head -5 | sed 's/^/    /' || true
    echo "  ... (${sample_count} samples total)"
    echo ""
    
    echo "NEXT STEP:"
    echo "  sbatch ${SCRIPTS_DIR}/step4_check_results.sh"
    echo ""
    echo "To download all results:"
    echo "  scp -r \${USER}@puhti.csc.fi:${RESULTS_DIR} ."
    
else
    echo "✗ Pipeline Failed!"
    echo "=========================================="
    echo ""
    echo "Exit code: ${EXIT_CODE}"
    echo ""
    echo "Check logs for errors:"
    echo "  - SLURM log:     ${LOGS_DIR}/hla_batch_bam_${SLURM_JOB_ID}.log"
    echo "  - Nextflow log:  ${WORK_NF_DIR}/.nextflow.log"
    echo ""
    echo "To resume pipeline:"
    echo "  sbatch ${SCRIPTS_DIR}/step3_run_pipeline_bam.sh"
fi

echo ""
echo "Completed: $(date)"
echo ""

exit ${EXIT_CODE}
