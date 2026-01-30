#!/usr/bin/env bash

#SBATCH --job-name=hla_typing_rnaseq

#SBATCH --account=project_2008084

#SBATCH --partition=small

#SBATCH --time=24:00:00

#SBATCH --mem=8G

#SBATCH --cpus-per-task=2

#SBATCH --output=../logs/hla_pipeline_%j.out

#SBATCH --error=../logs/hla_pipeline_%j.err



echo "=========================================="

echo "HLA Typing Pipeline - RNA-seq"

echo "=========================================="

echo "Job ID: ${SLURM_JOB_ID}"

echo "Started: $(date)"

echo ""



# Load modules

module purge

module load nextflow #/23.10.0

#module load singularity/3.11.4



# Set paths

BASE_DIR="/scratch/project_2008084/${USER}/hla_rnaseq_analysis"

PIPELINE_DIR="${BASE_DIR}/pipeline"

INPUT_DIR="${BASE_DIR}/pipeline_input"

RESULTS_DIR="${BASE_DIR}/results"

SINGULARITY_CACHE="${BASE_DIR}/singularity_cache"



# Create directories

mkdir -p ${RESULTS_DIR}

mkdir -p ${SINGULARITY_CACHE}



# Nextflow settings

export NXF_OPTS="-Xms512M -Xmx4G"

export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"



# Count samples

SAMPLE_COUNT=$(ls ${INPUT_DIR}/*_R1.fastq.gz 2>/dev/null | wc -l)

echo "Processing ${SAMPLE_COUNT} samples"

echo ""



# Show samples

echo "Samples:"

ls ${INPUT_DIR}/*_R1.fastq.gz | xargs -n1 basename | sed 's/_R1.fastq.gz//' | nl

echo ""



# Run pipeline

cd ${PIPELINE_DIR}



nextflow run main.nf \

    --input ${INPUT_DIR} \

    --input_type fastq \

    --tools optitype,arcashla,spechla \

    --optitype_seq_type rna \

    --arcashla_genes 'A,B,C,DRB1,DQB1,DPB1' \

    --spechla_genes 'A,B,C,DRB1,DQB1,DPB1' \

    --enable_majority_voting \

    --mv_min_tools 2 \

    --mv_resolution 2field \

    --mv_genes 'A,B,C,DRB1,DQB1,DPB1' \

    --slurm_account project_2008084 \

    --slurm_partition small \

    --max_cpus 40 \

    --max_memory 180.GB \

    --max_time 24.h \

    --outdir ${RESULTS_DIR} \

    -profile puhti,singularity \

    -resume \

    -with-report ${RESULTS_DIR}/pipeline_report.html \

    -with-timeline ${RESULTS_DIR}/pipeline_timeline.html \

    -with-dag ${RESULTS_DIR}/pipeline_dag.svg



EXIT_CODE=$?



echo ""

echo "=========================================="

if [ ${EXIT_CODE} -eq 0 ]; then

    echo "Pipeline completed successfully!"

    echo "Results: ${RESULTS_DIR}"

else

    echo "Pipeline failed with exit code: ${EXIT_CODE}"

    echo "Check logs:"

    echo "  - ${BASE_DIR}/logs/hla_pipeline_${SLURM_JOB_ID}.out"

    echo "  - ${BASE_DIR}/logs/hla_pipeline_${SLURM_JOB_ID}.err"

    echo "  - ${PIPELINE_DIR}/.nextflow.log"

fi



echo "Completed: $(date)"

exit ${EXIT_CODE}

