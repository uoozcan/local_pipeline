#!/bin/bash
#=============================================================================
# Standalone POLYSOLVER SLURM array for 1KGP calibration samples
#=============================================================================
# Runs POLYSOLVER v4 (Broad Institute) on 1KGP samples.
# Skips samples that already have a non-empty POLYSOLVER result.
#
# IMPORTANT: POLYSOLVER only types Class I (A, B, C).
#
# IMPORTANT: This script MUST run on Puhti (or any native Linux compute).
#   POLYSOLVER's bundled novoalign binary (2012 vintage) crashes on WSL2
#   kernels ≥ 5.15 with SIGSEGV due to the old binary's memory layout.
#   It works correctly on standard Linux (CentOS/RHEL/Ubuntu bare metal or VM).
#
# 1KGP CRAM compatibility:
#   1KGP CRAMs are hs37d5 (GRCh37). POLYSOLVER accepts hg19/hg38.
#   hs37d5 ≈ hg19 (same chr naming, no "chr" prefix) — set build=hg19.
#
# Usage (from login node):
#   bash run_polysolver_array.sh
#
# Outputs per sample:
#   1kgp_typing_results/{SAMPLE}/polysolver/{SAMPLE}_polysolver.txt
#   1kgp_typing_results/by_tool/polysolver/{SAMPLE}_polysolver.txt  (symlink)
#
# Prerequisites on Puhti:
#   - polysolver.sif at /projappl/project_2008084/containers/polysolver.sif
#     (pull with: singularity pull polysolver.sif docker://sachet/polysolver:v4)
#   - samtools loaded via module (for CRAM streaming)
#   - parse_polysolver_results.py in PIPELINE_BIN
#=============================================================================

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------
PROJECT_ID="project_2008084"
BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/polysolver_work"
HLA_TOOLS="/scratch/${PROJECT_ID}/hla_tools"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline_fresh/hla_typing_pipeline/bin"
POLYSOLVER_SIF="${HLA_TOOLS}/containers/polysolver.sif"
# 1KGP CRAMs are hs37d5 (hg37 = hg19 compatible)
POLYSOLVER_BUILD="hg19"
HLA_REGION="6:28000000-34000000"    # no "chr" prefix for hs37d5

# CRAM accessions
declare -A CRAM_ERR=(
    [NA19238]=ERR3239453
    [NA19239]=ERR3239454
    [NA18526]=ERR3239353
    [NA18542]=ERR3239356
    [HG00096]=ERR3240114
    [NA20502]=ERR3239785
)
EBI_BASE="ftp://ftp.sra.ebi.ac.uk/vol1/run"

#=============================================================================
# WRAPPER (runs on login node — discovers samples, submits array job)
#=============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    mkdir -p "$LOGS_DIR" "$WORK_BASE"
    mkdir -p "${RESULTS_DIR}/by_tool/polysolver"

    # Remove blank/invalid POLYSOLVER results so they get re-submitted
    echo "[INFO] Scanning for blank POLYSOLVER results to clean up..."
    CLEANED=0
    for F in "${RESULTS_DIR}/by_tool/polysolver/"*_polysolver.txt; do
        [[ -f "$F" ]] || [[ -L "$F" ]] || continue
        REAL_F="$(readlink -f "$F" 2>/dev/null || echo "$F")"
        if ! grep -qP '^[ABC]\t|^Gene\t' "$REAL_F" 2>/dev/null; then
            SAMPLE_C=$(basename "$F" _polysolver.txt)
            rm -f "$F"
            rm -f "${RESULTS_DIR}/${SAMPLE_C}/polysolver/${SAMPLE_C}_polysolver.txt"
            CLEANED=$(( CLEANED + 1 ))
        fi
    done
    echo "[INFO] Removed ${CLEANED} blank/invalid POLYSOLVER result(s)"

    # Discover samples typed by any other tool
    SAMPLE_LIST_FILE="${BASE}/conf/polysolver_pending.txt"
    mkdir -p "${BASE}/conf"
    : > "$SAMPLE_LIST_FILE"

    for TOOL in hlahd optitype arcashla spechla; do
        for F in "${RESULTS_DIR}"/*/; do
            SAMPLE=$(basename "$F")
            RESULT="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
            [[ -f "$RESULT" ]] || continue
            # Skip if POLYSOLVER already done
            DONE="${RESULTS_DIR}/by_tool/polysolver/${SAMPLE}_polysolver.txt"
            if [[ -s "$DONE" ]] && grep -qP '^[ABC]\t|^Gene\t' "$DONE" 2>/dev/null; then
                continue
            fi
            echo "$SAMPLE"
        done
    done | sort -u >> "$SAMPLE_LIST_FILE"

    N=$(wc -l < "$SAMPLE_LIST_FILE")

    if [[ "$N" -eq 0 ]]; then
        echo "[INFO] All samples already have POLYSOLVER results. Nothing to do."
        exit 0
    fi

    echo "[INFO] Samples to run: $N"
    cat "$SAMPLE_LIST_FILE"
    echo ""

    sbatch \
        --job-name=polysolver_1kgp \
        --account="${PROJECT_ID}" \
        --partition=small \
        --time=04:00:00 \
        --cpus-per-task=4 \
        --mem=8G \
        --array="1-${N}" \
        --output="${LOGS_DIR}/polysolver_%A_%a.out" \
        --error="${LOGS_DIR}/polysolver_%A_%a.err" \
        "$0"

    echo "[INFO] Array job submitted. Monitor with: squeue -u \$USER"
    exit 0
fi

#=============================================================================
# ARRAY TASK (runs inside SLURM)
#=============================================================================
set -euo pipefail

BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/polysolver_work"
HLA_TOOLS="/scratch/${PROJECT_ID}/hla_tools"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline_fresh/hla_typing_pipeline/bin"
POLYSOLVER_SIF="${HLA_TOOLS}/containers/polysolver.sif"
POLYSOLVER_BUILD="hg19"
HLA_REGION="6:28000000-34000000"
SAMPLE_LIST_FILE="${BASE}/conf/polysolver_pending.txt"

declare -A CRAM_ERR=(
    [NA19238]=ERR3239453
    [NA19239]=ERR3239454
    [NA18526]=ERR3239353
    [NA18542]=ERR3239356
    [HG00096]=ERR3240114
    [NA20502]=ERR3239785
)
EBI_BASE="ftp://ftp.sra.ebi.ac.uk/vol1/run"

#-----------------------------------------------------------------------------
# Resolve sample for this array task
#-----------------------------------------------------------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST_FILE")
if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: No sample at index ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST_FILE}"
    exit 1
fi

echo "=== POLYSOLVER: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Date: $(date)"

#-----------------------------------------------------------------------------
# Load modules
#-----------------------------------------------------------------------------
module purge
module load gcc
module load samtools/1.21
module load singularity

echo "samtools: $(samtools --version | head -1)"
echo "singularity: $(singularity --version)"

#-----------------------------------------------------------------------------
# Step 1: Stream HLA-region BAM from EBI CRAM
#   POLYSOLVER requires a coordinate-sorted BAM; stream region and sort.
#   CRAM is hs37d5 (hg19-like, no chr prefix).
#-----------------------------------------------------------------------------
WORKDIR="${WORK_BASE}/${SAMPLE}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

BAM="${WORKDIR}/${SAMPLE}.hla.bam"

ERR="${CRAM_ERR[$SAMPLE]:-}"
if [[ -z "$ERR" ]]; then
    echo "ERROR: No CRAM accession for ${SAMPLE}"
    exit 1
fi

PREFIX="${ERR:0:6}"
CRAM_URL="${EBI_BASE}/${PREFIX}/${ERR}/${SAMPLE}.final.cram"
echo "CRAM URL: ${CRAM_URL}"

echo "Streaming HLA region (${HLA_REGION})..."
samtools view -b -@ 4 "$CRAM_URL" "$HLA_REGION" \
    2>"${LOGS_DIR}/polysolver_stream_${SAMPLE}.log" \
    | samtools sort -@ 4 -o "$BAM" \
    2>>"${LOGS_DIR}/polysolver_stream_${SAMPLE}.log" \
    || { echo "ERROR: samtools stream/sort failed"; cat "${LOGS_DIR}/polysolver_stream_${SAMPLE}.log"; exit 1; }

samtools index "$BAM"
echo "[OK] BAM: $(du -sh $BAM | cut -f1) ($(samtools flagstat $BAM | head -1))"

#-----------------------------------------------------------------------------
# Step 2: Run POLYSOLVER
#   shell_call_hla_type args:
#     BAM  ethnicity  include_freq  build  format  insertCalc  outdir
#   - ethnicity=Unknown: use allele-frequency-free mode
#   - include_freq=1: normalise coverage (recommended for low-coverage samples)
#   - build=hg19: hs37d5 is hg19-compatible
#   - format=STDFQ: standard paired-end FASTQ-aligned BAM
#   - insertCalc=0: skip insert size calculation (saves time for small BAMs)
#-----------------------------------------------------------------------------
POLY_OUT="${WORKDIR}/polysolver_out"
mkdir -p "$POLY_OUT"

echo "[Step 2] Running POLYSOLVER (build=${POLYSOLVER_BUILD})..."

# SAMTOOLS_DIR must be set inside container for POLYSOLVER hg38/hg19 script path
singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$POLYSOLVER_SIF" \
    bash -c "
        export SAMTOOLS_DIR=/home/polysolver/binaries
        bash /home/polysolver/scripts/shell_call_hla_type \
            '${BAM}' \
            Unknown \
            1 \
            '${POLYSOLVER_BUILD}' \
            STDFQ \
            0 \
            '${POLY_OUT}'
    " 2>"${LOGS_DIR}/polysolver_run_${SAMPLE}.log" \
    || true   # POLYSOLVER exits non-zero when a gene is untyped

echo "POLYSOLVER finished."

# Show log tail for diagnostics
echo "--- POLYSOLVER log tail ---"
tail -20 "${LOGS_DIR}/polysolver_run_${SAMPLE}.log" 2>/dev/null || true
echo "---"

# Check for winners file
WINNERS="${POLY_OUT}/winners.hla.nofreq.txt"
if [[ ! -f "$WINNERS" ]]; then
    echo "ERROR: winners.hla.nofreq.txt not found in ${POLY_OUT}"
    echo "Directory contents:"
    ls -la "$POLY_OUT" 2>/dev/null || echo "(empty)"
    exit 1
fi

echo "[OK] Winners file found:"
cat "$WINNERS"

#-----------------------------------------------------------------------------
# Step 3: Parse winners.hla.nofreq.txt → standard pipeline TSV
#-----------------------------------------------------------------------------
RESULT_TSV="${WORKDIR}/${SAMPLE}_polysolver.txt"

python3 "${PIPELINE_BIN}/parse_polysolver_results.py" \
    --input "$WINNERS" \
    --sample "$SAMPLE" \
    --output "$RESULT_TSV" \
    || { echo "ERROR: parse_polysolver_results.py failed"; exit 1; }

echo "[OK] Parsed result: $RESULT_TSV"

# Preview output
echo "--- polysolver result preview ---"
head -8 "$RESULT_TSV" 2>/dev/null || echo "(empty)"
echo "---"

# Validate it contains gene calls
if ! grep -qP '^[ABC]\t' "$RESULT_TSV" 2>/dev/null; then
    echo "WARNING: No A/B/C allele calls in result for ${SAMPLE}"
    cat "$RESULT_TSV"
    # Do not exit — publish whatever was produced
fi

#-----------------------------------------------------------------------------
# Step 4: Publish result
#-----------------------------------------------------------------------------
OUT_DIR="${RESULTS_DIR}/${SAMPLE}/polysolver"
mkdir -p "$OUT_DIR"
cp "$RESULT_TSV" "${OUT_DIR}/${SAMPLE}_polysolver.txt"
echo "[OK] Published: ${OUT_DIR}/${SAMPLE}_polysolver.txt"

#-----------------------------------------------------------------------------
# Step 5: Symlink into by_tool layout for calibration
#-----------------------------------------------------------------------------
mkdir -p "${RESULTS_DIR}/by_tool/polysolver"
SRC="${OUT_DIR}/${SAMPLE}_polysolver.txt"
DST="${RESULTS_DIR}/by_tool/polysolver/${SAMPLE}_polysolver.txt"
[[ -e "$DST" ]] && rm -f "$DST"
ln -sf "$SRC" "$DST"
echo "[OK] Symlinked: ${DST}"

echo ""
echo "=== Done: ${SAMPLE} ==="
head -8 "${OUT_DIR}/${SAMPLE}_polysolver.txt"
