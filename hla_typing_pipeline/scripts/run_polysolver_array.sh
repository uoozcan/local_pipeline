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
FASTQ_DIR="${BASE}/1kgp_fastqs"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/polysolver_work"
HLA_TOOLS="/scratch/${PROJECT_ID}/hla_tools"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/bin"
# SIF location confirmed by user (singularity_cache layout)
POLYSOLVER_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/polysolver.sif"
# hs37d5 reference for FASTQ → BAM alignment (needed by POLYSOLVER)
# POLYSOLVER hg19 mode: extracts reads from chr6 without "chr" prefix
POLYSOLVER_BUILD="hg19"
HLA_REGION="6:28000000-34000000"    # no "chr" prefix for hs37d5
HG19_REF="/scratch/${PROJECT_ID}/references/hs37d5.fa"   # hg19/hs37d5 reference

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

    # polysolver_pending.txt persists between runs — only rebuilt when absent or --refresh passed.
    # After each batch is submitted, the batch samples are removed from this file so the next
    # run automatically picks up the next N samples (not the same ones again).
    SAMPLE_LIST_FILE="${BASE}/conf/polysolver_pending.txt"
    BATCH_LIST_FILE="${BASE}/conf/polysolver_batch.txt"
    mkdir -p "${BASE}/conf"

    if [[ ! -f "$SAMPLE_LIST_FILE" ]] || [[ "${1:-}" == "--refresh" ]]; then
        echo "[INFO] Building pending sample list..."
        : > "$SAMPLE_LIST_FILE"
        for TOOL in hlahd optitype arcashla spechla; do
            for F in "${RESULTS_DIR}"/*/; do
                SAMPLE=$(basename "$F")
                RESULT="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
                [[ -f "$RESULT" ]] || continue
                DONE="${RESULTS_DIR}/by_tool/polysolver/${SAMPLE}_polysolver.txt"
                if [[ -s "$DONE" ]] && grep -qP '^[ABC]\t|^Gene\t' "$DONE" 2>/dev/null; then
                    continue
                fi
                echo "$SAMPLE"
            done
        done | sort -u >> "$SAMPLE_LIST_FILE"
    fi

    N_TOTAL=$(wc -l < "$SAMPLE_LIST_FILE")

    if [[ "$N_TOTAL" -eq 0 ]]; then
        echo "[INFO] All samples done. Run with --refresh to rescan for any failures."
        rm -f "$SAMPLE_LIST_FILE"
        exit 0
    fi

    # Take first BATCH_SIZE samples and remove them from the pending list immediately
    # so the next run picks the next batch (not the same ones).
    BATCH_SIZE="${POLYSOLVER_BATCH_SIZE:-20}"
    N=$(( N_TOTAL < BATCH_SIZE ? N_TOTAL : BATCH_SIZE ))

    head -"$N" "$SAMPLE_LIST_FILE" > "$BATCH_LIST_FILE"
    tail -n +"$(( N + 1 ))" "$SAMPLE_LIST_FILE" > "${SAMPLE_LIST_FILE}.tmp"
    mv "${SAMPLE_LIST_FILE}.tmp" "$SAMPLE_LIST_FILE"

    N_REMAINING=$(wc -l < "$SAMPLE_LIST_FILE")
    echo "[INFO] Submitting batch of ${N} — ${N_REMAINING} samples remain for future batches"
    cat "$BATCH_LIST_FILE"
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
    if [[ "$N_REMAINING" -gt 0 ]]; then
        echo "[INFO] Re-run after this batch completes to submit the next ${BATCH_SIZE} samples (${N_REMAINING} remaining)."
    else
        echo "[INFO] This is the last batch. Run with --refresh after completion to check for failures."
    fi
    exit 0
fi

#=============================================================================
# ARRAY TASK (runs inside SLURM)
#=============================================================================
set -euo pipefail

BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
FASTQ_DIR="${BASE}/1kgp_fastqs"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/polysolver_work"
HLA_TOOLS="/scratch/${PROJECT_ID}/hla_tools"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/bin"
POLYSOLVER_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/polysolver.sif"
POLYSOLVER_BUILD="hg19"
HLA_REGION="6:28000000-34000000"
HG19_REF="/scratch/${PROJECT_ID}/references/hs37d5.fa"
SAMPLE_LIST_FILE="${BASE}/conf/polysolver_batch.txt"

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
# NOTE: singularity/apptainer is a system command on Puhti — no module needed
#-----------------------------------------------------------------------------
module purge
module load gcc
module load samtools/1.21
module load bwa   # needed for FASTQ → BAM alignment

echo "samtools:   $(samtools --version | head -1)"
echo "bwa:        $(bwa 2>&1 | head -1 || true)"
echo "singularity: $(singularity --version 2>/dev/null || apptainer --version 2>/dev/null || echo 'check path')"

#-----------------------------------------------------------------------------
# Step 1: Get coordinate-sorted HLA-region BAM for POLYSOLVER
#
# POLYSOLVER requires a coordinate-sorted BAM aligned to hg19 (hs37d5).
# Strategy:
#   a) If pre-existing FASTQs found in FASTQ_DIR → align to hs37d5 with bwa
#   b) If FASTQs missing but sample has known CRAM accession → stream from EBI
#-----------------------------------------------------------------------------
# Use Lustre scratch for work (TMPDIR on Puhti small-partition nodes is a tiny
# RAM-based /tmp — POLYSOLVER's ~4,400 temp BAMs overflow it).
# File quota is managed by aggressive cleanup after each sample (see below).
WORKDIR="${WORK_BASE}/${SAMPLE}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

BAM="${WORKDIR}/${SAMPLE}.hla.bam"
R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

if [[ -s "$R1" ]] && [[ -s "$R2" ]]; then
    echo "[Step 1] Aligning pre-existing FASTQs to hs37d5 → BAM..."
    if [[ ! -f "${HG19_REF}.bwt" ]]; then
        echo "ERROR: hs37d5 BWA index not found at ${HG19_REF}.bwt"
        echo "       Index with: bwa index ${HG19_REF}"
        exit 1
    fi
    bwa mem -t "${SLURM_CPUS_PER_TASK:-4}" "$HG19_REF" "$R1" "$R2" \
        2>"${LOGS_DIR}/polysolver_bwa_${SAMPLE}.log" \
        | samtools sort -@ 4 -o "$BAM" \
        2>>"${LOGS_DIR}/polysolver_bwa_${SAMPLE}.log" \
        || { echo "ERROR: bwa/samtools failed"; tail -5 "${LOGS_DIR}/polysolver_bwa_${SAMPLE}.log"; exit 1; }
    samtools index "$BAM"
    echo "[OK] BAM from FASTQs: $(du -sh $BAM | cut -f1) ($(samtools view -c $BAM) reads)"

else
    echo "[Step 1] FASTQs not found — streaming HLA region from EBI CRAM..."
    ERR="${CRAM_ERR[$SAMPLE]:-}"
    if [[ -z "$ERR" ]]; then
        echo "ERROR: No FASTQ in ${FASTQ_DIR} and no CRAM accession for ${SAMPLE}"
        exit 1
    fi
    PREFIX="${ERR:0:6}"
    CRAM_URL="${EBI_BASE}/${PREFIX}/${ERR}/${SAMPLE}.final.cram"
    echo "CRAM URL: ${CRAM_URL}"
    samtools view -b -@ 4 "$CRAM_URL" "$HLA_REGION" \
        2>"${LOGS_DIR}/polysolver_stream_${SAMPLE}.log" \
        | samtools sort -@ 4 -o "$BAM" \
        2>>"${LOGS_DIR}/polysolver_stream_${SAMPLE}.log" \
        || { echo "ERROR: samtools stream/sort failed"; cat "${LOGS_DIR}/polysolver_stream_${SAMPLE}.log"; exit 1; }
    samtools index "$BAM"
    echo "[OK] BAM from CRAM: $(du -sh $BAM | cut -f1) ($(samtools view -c $BAM) reads)"
fi

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
# include_freq=1 produces winners.hla.txt instead of winners.hla.nofreq.txt
[[ -f "$WINNERS" ]] || WINNERS="${POLY_OUT}/winners.hla.txt"
if [[ ! -f "$WINNERS" ]]; then
    echo "ERROR: neither winners.hla.nofreq.txt nor winners.hla.txt found in ${POLY_OUT}"
    echo "Directory contents:"
    ls -la "$POLY_OUT" 2>/dev/null || echo "(empty)"
    exit 1
fi

echo "[OK] Winners file found:"
cat "$WINNERS"

# Delete all temp BAMs/SAMs immediately — keeps ~4,400 files off scratch quota.
# Winners file is preserved; everything else in polysolver_out can go.
echo "[Cleanup] Removing temp BAMs from polysolver_out..."
find "$POLY_OUT" ! -name 'winners*' -type f -delete
echo "[Cleanup] Done ($(find "$POLY_OUT" | wc -l) files remaining in polysolver_out)"

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

# Clean up entire work dir — result is published, nothing needed here anymore
rm -rf "$WORKDIR"
echo "[Cleanup] Removed WORKDIR: ${WORKDIR}"

echo ""
echo "=== Done: ${SAMPLE} ==="
head -8 "${OUT_DIR}/${SAMPLE}_polysolver.txt"
