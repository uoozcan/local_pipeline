#!/bin/bash
#=============================================================================
# Standalone seq2HLA SLURM array for 1KGP calibration samples
#=============================================================================
# Runs seq2HLA v2.3 on all 1KGP samples already typed by other tools.
# Skips samples that already have a non-empty seq2HLA result.
#
# seq2HLA takes FASTQ input only (RNA-seq optimal; WGS will give lower
# confidence scores but will still type A/B/C/DRB1/DQA1/DQB1/DPA1/DPB1).
#
# Usage (from login node):
#   bash run_seq2hla_array.sh
#
# Outputs per sample:
#   1kgp_typing_results/{SAMPLE}/seq2hla/{SAMPLE}_seq2hla.txt
#   1kgp_typing_results/by_tool/seq2hla/{SAMPLE}_seq2hla.txt  (symlink)
#
# Prerequisites on Puhti:
#   - seq2hla.sif at /projappl/project_2008084/containers/seq2hla.sif
#   - Pre-existing 1kgp_fastqs/ OR EBI CRAM accessions (auto-stream fallback)
#   - Python 3 available on host for parse_seq2hla_results.py
#=============================================================================

#-----------------------------------------------------------------------------
# Configuration — edit these paths if your layout differs
#-----------------------------------------------------------------------------
PROJECT_ID="project_2008084"
BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
FASTQ_DIR="${BASE}/1kgp_fastqs"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/seq2hla_work"
HLA_TOOLS="/scratch/${PROJECT_ID}/hla_tools"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/hla_typing_pipeline/bin"
SEQ2HLA_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/seq2hla.sif"
HLA_REGION="chr6:28000000-34000000"

# CRAM accessions (fallback when FASTQ is missing)
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
    mkdir -p "${RESULTS_DIR}/by_tool/seq2hla"

    # Remove blank/invalid seq2HLA results so they get re-submitted
    echo "[INFO] Scanning for blank seq2HLA results to clean up..."
    CLEANED=0
    for F in "${RESULTS_DIR}/by_tool/seq2hla/"*_seq2hla.txt; do
        [[ -f "$F" ]] || [[ -L "$F" ]] || continue
        REAL_F="$(readlink -f "$F" 2>/dev/null || echo "$F")"
        if ! grep -qP 'HLA[*:]' "$REAL_F" 2>/dev/null; then
            SAMPLE_C=$(basename "$F" _seq2hla.txt)
            rm -f "$F"
            rm -f "${RESULTS_DIR}/${SAMPLE_C}/seq2hla/${SAMPLE_C}_seq2hla.txt"
            CLEANED=$(( CLEANED + 1 ))
        fi
    done
    echo "[INFO] Removed ${CLEANED} blank/invalid seq2HLA result(s)"

    SAMPLE_LIST_FILE="${BASE}/conf/seq2hla_pending.txt"
    BATCH_LIST_FILE="${BASE}/conf/seq2hla_batch.txt"
    mkdir -p "${BASE}/conf"

    if [[ ! -f "$SAMPLE_LIST_FILE" ]] || [[ "${1:-}" == "--refresh" ]]; then
        echo "[INFO] Building pending sample list..."
        : > "$SAMPLE_LIST_FILE"
        for TOOL in hlahd optitype arcashla spechla; do
            for F in "${RESULTS_DIR}"/*/; do
                SAMPLE=$(basename "$F")
                RESULT="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
                [[ -f "$RESULT" ]] || continue
                DONE="${RESULTS_DIR}/by_tool/seq2hla/${SAMPLE}_seq2hla.txt"
                if [[ -s "$DONE" ]] && grep -qP 'HLA[*:]' "$DONE" 2>/dev/null; then
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

    BATCH_SIZE="${SEQ2HLA_BATCH_SIZE:-20}"
    N=$(( N_TOTAL < BATCH_SIZE ? N_TOTAL : BATCH_SIZE ))

    head -"$N" "$SAMPLE_LIST_FILE" > "$BATCH_LIST_FILE"
    tail -n +"$(( N + 1 ))" "$SAMPLE_LIST_FILE" > "${SAMPLE_LIST_FILE}.tmp"
    mv "${SAMPLE_LIST_FILE}.tmp" "$SAMPLE_LIST_FILE"

    N_REMAINING=$(wc -l < "$SAMPLE_LIST_FILE")
    echo "[INFO] Submitting batch of ${N} — ${N_REMAINING} samples remain for future batches"
    cat "$BATCH_LIST_FILE"
    echo ""

    sbatch \
        --job-name=seq2hla_1kgp \
        --account="${PROJECT_ID}" \
        --partition=small \
        --time=02:00:00 \
        --cpus-per-task=4 \
        --mem=8G \
        --array="1-${N}" \
        --output="${LOGS_DIR}/seq2hla_%A_%a.out" \
        --error="${LOGS_DIR}/seq2hla_%A_%a.err" \
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
WORK_BASE="${BASE}/seq2hla_work"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/hla_typing_pipeline/bin"
SEQ2HLA_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/seq2hla.sif"
HLA_REGION="chr6:28000000-34000000"
SAMPLE_LIST_FILE="${BASE}/conf/seq2hla_batch.txt"

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

echo "=== seq2HLA: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Date: $(date)"

#-----------------------------------------------------------------------------
# Load modules
#-----------------------------------------------------------------------------
module purge
module load gcc
module load samtools/1.21
# NOTE: singularity/apptainer is a system command on Puhti — no module needed

#-----------------------------------------------------------------------------
# Step 1: Ensure R1/R2 FASTQs exist
#-----------------------------------------------------------------------------
R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

if [[ -s "$R1" ]] && [[ -s "$R2" ]]; then
    echo "[OK] Using existing FASTQs"
else
    echo "[INFO] FASTQs not found — streaming HLA region from EBI CRAM"

    ERR="${CRAM_ERR[$SAMPLE]:-}"
    if [[ -z "$ERR" ]]; then
        echo "ERROR: No CRAM accession for ${SAMPLE} and no FASTQ available"
        exit 1
    fi

    PREFIX="${ERR:0:6}"
    CRAM_URL="${EBI_BASE}/${PREFIX}/${ERR}/${SAMPLE}.final.cram"
    echo "CRAM URL: ${CRAM_URL}"

    TMPBAM="${WORK_BASE}/${SAMPLE}_hla.bam"
    NAMESORT="${WORK_BASE}/${SAMPLE}_namesort.bam"
    mkdir -p "$WORK_BASE"

    echo "Streaming HLA region..."
    samtools view -b -@ 4 -o "${TMPBAM}.tmp" "$CRAM_URL" "$HLA_REGION" \
        2>"${LOGS_DIR}/seq2hla_stream_${SAMPLE}.log" \
        || { echo "ERROR: samtools stream failed"; cat "${LOGS_DIR}/seq2hla_stream_${SAMPLE}.log"; exit 1; }
    mv "${TMPBAM}.tmp" "$TMPBAM"

    echo "Name-sorting..."
    samtools sort -n -@ 4 "$TMPBAM" -o "$NAMESORT" \
        2>>"${LOGS_DIR}/seq2hla_stream_${SAMPLE}.log"
    rm -f "$TMPBAM"

    echo "Converting to FASTQ..."
    mkdir -p "$FASTQ_DIR"
    samtools fastq -@ 4 \
        -1 "$R1" -2 "$R2" \
        -0 /dev/null -s /dev/null \
        "$NAMESORT" \
        2>>"${LOGS_DIR}/seq2hla_stream_${SAMPLE}.log"
    rm -f "$NAMESORT"

    echo "[OK] FASTQ written: $(du -sh $R1 | cut -f1) / $(du -sh $R2 | cut -f1)"
fi

#-----------------------------------------------------------------------------
# Step 2: Create workdir and run seq2HLA
#-----------------------------------------------------------------------------
WORKDIR="${WORK_BASE}/${SAMPLE}"
mkdir -p "$WORKDIR"

echo "[Step 2] Running seq2HLA..."
cd "$WORKDIR"

# seq2HLA writes output files relative to cwd; -r sets the prefix
singularity exec "$SEQ2HLA_SIF" \
    seq2HLA \
        -r "${SAMPLE}." \
        -p "${SLURM_CPUS_PER_TASK}" \
        -1 "$(realpath "$R1")" \
        -2 "$(realpath "$R2")" \
    || true   # seq2HLA can exit non-zero on low-confidence samples

echo "seq2HLA finished."

# List what was produced
echo "--- seq2HLA output files ---"
ls -la "${SAMPLE}"* 2>/dev/null || ls -la 2>/dev/null | grep "${SAMPLE}" || echo "(no output files found)"
echo "---"

#-----------------------------------------------------------------------------
# Step 3: Parse results → standard pipeline TSV
#  parse_seq2hla_results.py is Python 2/3 compatible; call via host python3
#  or fall back to python inside container if needed
#-----------------------------------------------------------------------------
RESULT_TSV="${WORKDIR}/${SAMPLE}_seq2hla.txt"

python3 "${PIPELINE_BIN}/parse_seq2hla_results.py" \
    --sample "$SAMPLE" \
    --prefix "${WORKDIR}/${SAMPLE}." \
    --output "$RESULT_TSV" \
    || { echo "ERROR: parse_seq2hla_results.py failed"; exit 1; }

echo "[OK] Parsed result: $RESULT_TSV"

# Preview output
echo "--- seq2hla result preview ---"
head -6 "$RESULT_TSV" 2>/dev/null || echo "(empty)"
echo "---"

# Validate it contains actual allele calls
if ! grep -qP 'HLA[*:]' "$RESULT_TSV" 2>/dev/null && \
   ! grep -qP '^[ABC]\t' "$RESULT_TSV" 2>/dev/null; then
    echo "WARNING: ${SAMPLE}_seq2hla.txt has no HLA allele calls — low confidence or WGS data"
    echo "File contents:"
    cat "$RESULT_TSV"
    # Do not exit — allow empty results to be published (calibration will skip this sample)
fi

#-----------------------------------------------------------------------------
# Step 4: Publish result
#-----------------------------------------------------------------------------
OUT_DIR="${RESULTS_DIR}/${SAMPLE}/seq2hla"
mkdir -p "$OUT_DIR"
cp "$RESULT_TSV" "${OUT_DIR}/${SAMPLE}_seq2hla.txt"
echo "[OK] Published: ${OUT_DIR}/${SAMPLE}_seq2hla.txt"

#-----------------------------------------------------------------------------
# Step 5: Symlink into by_tool layout for calibration
#-----------------------------------------------------------------------------
mkdir -p "${RESULTS_DIR}/by_tool/seq2hla"
SRC="${OUT_DIR}/${SAMPLE}_seq2hla.txt"
DST="${RESULTS_DIR}/by_tool/seq2hla/${SAMPLE}_seq2hla.txt"
[[ -e "$DST" ]] && rm -f "$DST"
ln -sf "$SRC" "$DST"
echo "[OK] Symlinked: ${DST}"

echo ""
echo "=== Done: ${SAMPLE} ==="
head -6 "${OUT_DIR}/${SAMPLE}_seq2hla.txt"
