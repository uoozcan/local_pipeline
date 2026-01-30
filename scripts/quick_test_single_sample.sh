#!/bin/bash
#SBATCH --job-name=hla_quick_test
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=quick_test_%j.out
#SBATCH --error=quick_test_%j.err

set -e

echo "======================================================================"
echo "   HLA PIPELINE - QUICK SINGLE SAMPLE TEST"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""

# Load modules
#module load java/21
#module load biopython-env/3.10.6
module load biokit nextflow #/25.10.0

echo "Modules loaded:"
module list
echo ""

# Configuration
PIPELINE_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline"
INPUT_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/raw_fastq"
TEST_OUTPUT="./quick_test_results"
CONTAINER_CACHE="/scratch/project_2008084/hla_references/singularity_cache/containers"

# Find a single sample to test
echo "Finding test sample..."
TEST_R1=$(find $INPUT_DIR -name "*_R1_*.fastq.gz" | head -1)
TEST_R2=$(echo $TEST_R1 | sed 's/_R1_/_R2_/')

if [ -z "$TEST_R1" ] || [ ! -f "$TEST_R1" ]; then
    echo "ERROR: Could not find test sample in $INPUT_DIR"
    exit 1
fi

echo "Test sample R1: $TEST_R1"
echo "Test sample R2: $TEST_R2"
echo ""

# Extract sample name
SAMPLE_NAME=$(basename $TEST_R1 | sed 's/_R[12]_.*//')
echo "Sample name: $SAMPLE_NAME"
echo ""

# Create test input directory with only this sample
mkdir -p test_input
ln -sf $TEST_R1 test_input/
ln -sf $TEST_R2 test_input/

echo "======================================================================"
echo "TEST 1: ArcasHLA only"
echo "======================================================================"
echo ""

rm -rf ${TEST_OUTPUT}_arcashla
nextflow run ${PIPELINE_DIR}/main.nf \
    --input test_input \
    --input_type fastq \
    --tools arcashla \
    --outdir ${TEST_OUTPUT}_arcashla \
    --container_cache $CONTAINER_CACHE \
    --max_cpus 4 \
    --max_memory 32.GB \
    -profile puhti,singularity \
    -resume

echo ""
echo "ArcasHLA Results:"
find ${TEST_OUTPUT}_arcashla -name "*.json" -o -name "*.log" | while read f; do
    echo "--- $f ---"
    cat "$f" | head -50
    echo ""
done

echo ""
echo "======================================================================"
echo "TEST 2: OptiType only"
echo "======================================================================"
echo ""

rm -rf ${TEST_OUTPUT}_optitype
nextflow run ${PIPELINE_DIR}/main.nf \
    --input test_input \
    --input_type fastq \
    --tools optitype \
    --optitype_seq_type rna \
    --outdir ${TEST_OUTPUT}_optitype \
    --container_cache $CONTAINER_CACHE \
    --max_cpus 4 \
    --max_memory 32.GB \
    -profile puhti,singularity \
    -resume

echo ""
echo "OptiType Results:"
find ${TEST_OUTPUT}_optitype -name "*result.tsv" -o -name "*.log" | while read f; do
    echo "--- $f ---"
    cat "$f" | head -50
    echo ""
done

echo ""
echo "======================================================================"
echo "TEST 3: SpecHLA only"
echo "======================================================================"
echo ""

rm -rf ${TEST_OUTPUT}_spechla
nextflow run ${PIPELINE_DIR}/main.nf \
    --input test_input \
    --input_type fastq \
    --tools spechla \
    --spechla_exon_only 0 \
    --outdir ${TEST_OUTPUT}_spechla \
    --container_cache $CONTAINER_CACHE \
    --max_cpus 4 \
    --max_memory 32.GB \
    -profile puhti,singularity \
    -resume

echo ""
echo "SpecHLA Results:"
find ${TEST_OUTPUT}_spechla -name "*.tsv" -o -name "*.log" | while read f; do
    echo "--- $f ---"
    cat "$f" | head -50
    echo ""
done

echo ""
echo "======================================================================"
echo "SUMMARY"
echo "======================================================================"
echo ""

echo "Test completed: $(date)"
echo ""

echo "Result locations:"
echo "  ArcasHLA: ${TEST_OUTPUT}_arcashla"
echo "  OptiType: ${TEST_OUTPUT}_optitype"
echo "  SpecHLA: ${TEST_OUTPUT}_spechla"
echo ""

echo "Check results for:"
echo "1. Non-empty JSON/TSV files"
echo "2. Actual HLA alleles (not all '-')"
echo "3. Reasonable read counts in logs"
echo ""

# Cleanup
rm -rf test_input

echo "Quick test complete!"
