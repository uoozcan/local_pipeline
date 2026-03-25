#!/bin/bash
#=============================================================================
# Standalone Kourami SLURM array for 1KGP calibration samples
#=============================================================================
# Runs Kourami v0.9.6 (assembly-graph HLA typing) on 1KGP samples.
# Skips samples that already have a non-empty Kourami result.
#
# Input strategy for 1KGP CRAMs (hs37d5 / GRCh37):
#   1KGP CRAMs use chromosome name "6" (no "chr" prefix).
#   Kourami's alignAndExtract script requires an hs38NoAltDH reference, but
#   the initial extraction is just a samtools region extraction.
#   We use the "direct samtools extraction" path in kourami.nf — no hs38 ref
#   needed. Reads are extracted from 6:28000000-34000000, converted to FASTQ,
#   then aligned to the Kourami HLA panel (which is IMGT-based, not chr-specific).
#   This works correctly regardless of whether the input was hg37 or hg38.
#
# Usage (from login node):
#   bash run_kourami_array.sh
#
# Outputs per sample:
#   1kgp_typing_results/{SAMPLE}/kourami/{SAMPLE}_kourami.txt
#   1kgp_typing_results/by_tool/kourami/{SAMPLE}_kourami.txt  (symlink)
#
# Prerequisites on Puhti:
#   - Kourami jar at /projappl/project_2008084/kourami/kourami-0.9.6/build/Kourami.jar
#   - Kourami DB at /projappl/project_2008084/kourami/kourami_db/
#     (must contain All_FINAL_with_Decoy.fa.gz + bwa index)
#   - samtools, bwa, java (≥11) loaded via modules
#   - parse_kourami_results.py in PIPELINE_BIN
#=============================================================================

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------
PROJECT_ID="project_2008084"
BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
FASTQ_DIR="${BASE}/1kgp_fastqs"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/kourami_work"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/bin"
KOURAMI_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/kourami.sif"
# JAR and DB are bundled inside the zlskidmore/kourami container
KOURAMI_JAR="/usr/local/bin/Kourami.jar"
KOURAMI_DB="/usr/local/bin/kourami-0.9.6/db"
HLA_REGION_NOPREFIX="6:28000000-34000000"    # hs37d5 uses no "chr" prefix
HLA_REGION_CHR="chr6:28000000-34000000"      # hg38 uses "chr" prefix

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
    mkdir -p "${RESULTS_DIR}/by_tool/kourami"

    # Remove blank/invalid Kourami results so they get re-submitted
    echo "[INFO] Scanning for blank Kourami results to clean up..."
    CLEANED=0
    for F in "${RESULTS_DIR}/by_tool/kourami/"*_kourami.txt; do
        [[ -f "$F" ]] || [[ -L "$F" ]] || continue
        REAL_F="$(readlink -f "$F" 2>/dev/null || echo "$F")"
        if ! grep -qP 'HLA[*:]|^[ABC]\t' "$REAL_F" 2>/dev/null; then
            SAMPLE_C=$(basename "$F" _kourami.txt)
            rm -f "$F"
            rm -f "${RESULTS_DIR}/${SAMPLE_C}/kourami/${SAMPLE_C}_kourami.txt"
            CLEANED=$(( CLEANED + 1 ))
        fi
    done
    echo "[INFO] Removed ${CLEANED} blank/invalid Kourami result(s)"

    # kourami_pending.txt persists between runs — only rebuilt when absent or --refresh passed.
    # After each batch is submitted, the batch samples are removed from this file so the next
    # run automatically picks up the next N samples (not the same ones again).
    SAMPLE_LIST_FILE="${BASE}/conf/kourami_pending.txt"
    BATCH_LIST_FILE="${BASE}/conf/kourami_batch.txt"
    mkdir -p "${BASE}/conf"

    if [[ ! -f "$SAMPLE_LIST_FILE" ]] || [[ "${1:-}" == "--refresh" ]]; then
        echo "[INFO] Building pending sample list..."
        : > "$SAMPLE_LIST_FILE"
        for TOOL in hlahd optitype arcashla spechla; do
            for F in "${RESULTS_DIR}"/*/; do
                SAMPLE=$(basename "$F")
                RESULT="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
                [[ -f "$RESULT" ]] || continue
                DONE="${RESULTS_DIR}/by_tool/kourami/${SAMPLE}_kourami.txt"
                if [[ -s "$DONE" ]] && grep -qP 'HLA[*:]|^[ABC]\t' "$DONE" 2>/dev/null; then
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

    BATCH_SIZE="${KOURAMI_BATCH_SIZE:-20}"
    N=$(( N_TOTAL < BATCH_SIZE ? N_TOTAL : BATCH_SIZE ))

    head -"$N" "$SAMPLE_LIST_FILE" > "$BATCH_LIST_FILE"
    tail -n +"$(( N + 1 ))" "$SAMPLE_LIST_FILE" > "${SAMPLE_LIST_FILE}.tmp"
    mv "${SAMPLE_LIST_FILE}.tmp" "$SAMPLE_LIST_FILE"

    N_REMAINING=$(wc -l < "$SAMPLE_LIST_FILE")
    echo "[INFO] Submitting batch of ${N} — ${N_REMAINING} samples remain for future batches"
    cat "$BATCH_LIST_FILE"
    echo ""

    sbatch \
        --job-name=kourami_1kgp \
        --account="${PROJECT_ID}" \
        --partition=small \
        --time=08:00:00 \
        --cpus-per-task=8 \
        --mem=14G \
        --array="1-${N}" \
        --output="${LOGS_DIR}/kourami_%A_%a.out" \
        --error="${LOGS_DIR}/kourami_%A_%a.err" \
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
WORK_BASE="${BASE}/kourami_work"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/bin"
KOURAMI_SIF="/scratch/${PROJECT_ID}/hla_references/singularity_cache/containers/kourami.sif"
KOURAMI_JAR="/usr/local/bin/Kourami.jar"
KOURAMI_DB="/usr/local/bin/kourami-0.9.6/db"
HLA_REGION_NOPREFIX="6:28000000-34000000"
HLA_REGION_CHR="chr6:28000000-34000000"
SAMPLE_LIST_FILE="${BASE}/conf/kourami_batch.txt"

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

echo "=== Kourami: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Date: $(date)"

#-----------------------------------------------------------------------------
# Load modules
# NOTE: bwa and java come from the kourami.sif container; only samtools is
# loaded from modules (container does not include samtools)
#-----------------------------------------------------------------------------
module purge
module load gcc
module load samtools/1.21

echo "samtools:    $(samtools --version | head -1)"
echo "singularity: $(singularity --version 2>/dev/null || apptainer --version 2>/dev/null || echo 'check path')"
echo "bwa (container):  $(singularity exec "$KOURAMI_SIF" bwa 2>&1 | head -1 || true)"
echo "java (container): $(singularity exec "$KOURAMI_SIF" java -version 2>&1 | head -1 || true)"

#-----------------------------------------------------------------------------
# Step 1: Create workdir; get HLA region FASTQ
#   - Use pre-existing FASTQs if available (reuse from other tools)
#   - Otherwise stream from EBI CRAM and extract HLA region
#-----------------------------------------------------------------------------
WORKDIR="${WORK_BASE}/${SAMPLE}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

R1_IN="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2_IN="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

# Kourami-specific HLA FASTQ (subsetted to chr6 HLA region)
R1="${WORKDIR}/${SAMPLE}._hla_1.fq.gz"
R2="${WORKDIR}/${SAMPLE}._hla_2.fq.gz"

if [[ -s "$R1" ]] && [[ -s "$R2" ]]; then
    echo "[OK] Using existing HLA FASTQ"
elif [[ -s "$R1_IN" ]] && [[ -s "$R2_IN" ]]; then
    echo "[INFO] Using full genome FASTQ from ${FASTQ_DIR} (no need to re-extract from CRAM)"
    # Kourami works best with HLA-region reads; use the full FASTQ as-is
    # (already ~HLA-enriched from the previous CRAM extraction step)
    R1="$R1_IN"
    R2="$R2_IN"
else
    echo "[INFO] Streaming HLA region from EBI CRAM"

    ERR="${CRAM_ERR[$SAMPLE]:-}"
    if [[ -z "$ERR" ]]; then
        echo "ERROR: No CRAM accession for ${SAMPLE} and no FASTQ available"
        exit 1
    fi

    PREFIX="${ERR:0:6}"
    CRAM_URL="${EBI_BASE}/${PREFIX}/${ERR}/${SAMPLE}.final.cram"
    echo "CRAM URL: ${CRAM_URL}"

    TMPBAM="${WORKDIR}/${SAMPLE}_hla.bam"
    NAMESORT="${WORKDIR}/${SAMPLE}_namesort.bam"

    # 1KGP CRAMs are hs37d5 (no "chr" prefix); try both naming styles
    echo "Detecting chromosome naming in CRAM..."
    if samtools view -H "$CRAM_URL" 2>/dev/null | grep -qP '^@SQ.*SN:chr6\t'; then
        HLA_REGION_USE="$HLA_REGION_CHR"
    else
        HLA_REGION_USE="$HLA_REGION_NOPREFIX"
    fi
    echo "Using region: ${HLA_REGION_USE}"

    echo "Streaming HLA region..."
    samtools view -b -@ 4 -o "${TMPBAM}.tmp" "$CRAM_URL" "$HLA_REGION_USE" \
        2>"${LOGS_DIR}/kourami_stream_${SAMPLE}.log" \
        || { echo "ERROR: samtools stream failed"; cat "${LOGS_DIR}/kourami_stream_${SAMPLE}.log"; exit 1; }
    mv "${TMPBAM}.tmp" "$TMPBAM"

    echo "Name-sorting..."
    samtools sort -n -@ 4 "$TMPBAM" -o "$NAMESORT" \
        2>>"${LOGS_DIR}/kourami_stream_${SAMPLE}.log"
    rm -f "$TMPBAM"

    echo "Converting to FASTQ..."
    mkdir -p "$FASTQ_DIR"
    samtools fastq -@ 4 \
        -1 "$R1" -2 "$R2" \
        -0 /dev/null -s /dev/null \
        "$NAMESORT" \
        2>>"${LOGS_DIR}/kourami_stream_${SAMPLE}.log"
    rm -f "$NAMESORT"

    echo "[OK] HLA FASTQ: $(du -sh $R1 | cut -f1) / $(du -sh $R2 | cut -f1)"
fi

#-----------------------------------------------------------------------------
# Step 2: Align HLA reads to Kourami panel (bwa inside container)
#-----------------------------------------------------------------------------
echo "[Step 2] Aligning to Kourami panel..."

singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    --bind "${FASTQ_DIR}:${FASTQ_DIR}" \
    "$KOURAMI_SIF" \
    bwa mem -t "${SLURM_CPUS_PER_TASK}" \
        "${KOURAMI_DB}/All_FINAL_with_Decoy.fa.gz" \
        "$R1" "$R2" \
    2>"${LOGS_DIR}/kourami_bwa_${SAMPLE}.log" \
    | samtools sort -@ "${SLURM_CPUS_PER_TASK}" \
        -o "${SAMPLE}.panel.bam"

samtools index "${SAMPLE}.panel.bam"
echo "[OK] Panel BAM: $(du -sh ${SAMPLE}.panel.bam | cut -f1)"

#-----------------------------------------------------------------------------
# Step 3: Kourami assembly-graph typing (java inside container)
#-----------------------------------------------------------------------------
echo "[Step 3] Running Kourami..."

singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$KOURAMI_SIF" \
    java -Xmx12g -jar "$KOURAMI_JAR" \
        -d "$KOURAMI_DB" \
        "${SAMPLE}.panel.bam" \
        -o "${SAMPLE}.kourami" \
    2>"${LOGS_DIR}/kourami_run_${SAMPLE}.log" \
    || true   # Kourami can exit non-zero if a gene fails

echo "Kourami finished."

# Preview .result file
echo "--- ${SAMPLE}.kourami.result preview ---"
head -10 "${SAMPLE}.kourami.result" 2>/dev/null || echo "(result file not found)"
echo "---"

if [[ ! -f "${SAMPLE}.kourami.result" ]]; then
    echo "ERROR: Kourami produced no .result file for ${SAMPLE}"
    echo "Kourami log tail:"
    tail -20 "${LOGS_DIR}/kourami_run_${SAMPLE}.log" 2>/dev/null || true
    exit 1
fi

#-----------------------------------------------------------------------------
# Step 4: Parse .result → standard pipeline TSV
#-----------------------------------------------------------------------------
RESULT_TSV="${WORKDIR}/${SAMPLE}_kourami.txt"

python3 "${PIPELINE_BIN}/parse_kourami_results.py" \
    --input "${SAMPLE}.kourami.result" \
    --sample "$SAMPLE" \
    --output "$RESULT_TSV" \
    || { echo "ERROR: parse_kourami_results.py failed"; exit 1; }

echo "[OK] Parsed result: $RESULT_TSV"

# Preview output
echo "--- kourami result preview ---"
head -8 "$RESULT_TSV" 2>/dev/null || echo "(empty)"
echo "---"

#-----------------------------------------------------------------------------
# Cleanup large intermediates
#-----------------------------------------------------------------------------
rm -f "${SAMPLE}.panel.bam" "${SAMPLE}.panel.bam.bai"
# Keep HLA FASTQs in WORK_BASE if they were newly generated (reusable)

#-----------------------------------------------------------------------------
# Step 5: Publish result
#-----------------------------------------------------------------------------
OUT_DIR="${RESULTS_DIR}/${SAMPLE}/kourami"
mkdir -p "$OUT_DIR"
cp "$RESULT_TSV" "${OUT_DIR}/${SAMPLE}_kourami.txt"
echo "[OK] Published: ${OUT_DIR}/${SAMPLE}_kourami.txt"

#-----------------------------------------------------------------------------
# Step 6: Symlink into by_tool layout for calibration
#-----------------------------------------------------------------------------
mkdir -p "${RESULTS_DIR}/by_tool/kourami"
SRC="${OUT_DIR}/${SAMPLE}_kourami.txt"
DST="${RESULTS_DIR}/by_tool/kourami/${SAMPLE}_kourami.txt"
[[ -e "$DST" ]] && rm -f "$DST"
ln -sf "$SRC" "$DST"
echo "[OK] Symlinked: ${DST}"

echo ""
echo "=== Done: ${SAMPLE} ==="
head -8 "${OUT_DIR}/${SAMPLE}_kourami.txt"
