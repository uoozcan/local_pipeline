#!/bin/bash
#=============================================================================
# Standalone SpecHLA SLURM array for 1KGP calibration samples
#=============================================================================
# Runs SpecHLA (local install at /projappl/project_2008084/SpecHLAx) on all
# samples already typed by other tools (hlahd/optitype/arcashla).
# Skips samples that already have a non-empty SpecHLA result.
#
# Usage (from login node):
#   bash run_spechla_array.sh
#
# Outputs per sample:
#   1kgp_typing_results/{SAMPLE}/spechla/{SAMPLE}_spechla.txt
#   1kgp_typing_results/by_tool/spechla/{SAMPLE}_spechla.txt  (symlink)
#=============================================================================

#-----------------------------------------------------------------------------
# Configuration — edit these paths if your layout differs
#-----------------------------------------------------------------------------
PROJECT_ID="project_2008084"
SPEC_HOME="/projappl/${PROJECT_ID}/SpecHLAx"
BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
FASTQ_DIR="${BASE}/1kgp_fastqs"
RESULTS_DIR="${BASE}/1kgp_typing_results"
LOGS_DIR="${BASE}/logs"
WORK_BASE="${BASE}/spechla_work"
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
    mkdir -p "${RESULTS_DIR}/by_tool/spechla"

    # Collect any existing results from spechla_work/ into the results layout
    echo "[INFO] Collecting existing SpecHLA results from ${WORK_BASE}..."
    COLLECTED=0
    for SAMPLE_DIR in "${WORK_BASE}"/*/; do
        [[ -d "$SAMPLE_DIR" ]] || continue
        SAMPLE=$(basename "$SAMPLE_DIR")
        # Find hla.result.txt — SpecHLA may nest it one or two levels deep
        RESULT_TXT=""
        for CANDIDATE in \
            "${SAMPLE_DIR}/${SAMPLE}/hla.result.txt" \
            "${SAMPLE_DIR}/${SAMPLE}/${SAMPLE}/hla.result.txt"; do
            if [[ -f "$CANDIDATE" ]] && [[ $(wc -l < "$CANDIDATE") -gt 1 ]]; then
                RESULT_TXT="$CANDIDATE"
                break
            fi
        done
        [[ -z "$RESULT_TXT" ]] && continue
        # Publish to results layout
        mkdir -p "${RESULTS_DIR}/${SAMPLE}/spechla"
        cp "$RESULT_TXT" "${RESULTS_DIR}/${SAMPLE}/spechla/${SAMPLE}_spechla.txt"
        # Symlink to by_tool
        DST="${RESULTS_DIR}/by_tool/spechla/${SAMPLE}_spechla.txt"
        [[ -e "$DST" ]] && rm -f "$DST"
        ln -sf "${RESULTS_DIR}/${SAMPLE}/spechla/${SAMPLE}_spechla.txt" "$DST"
        COLLECTED=$(( COLLECTED + 1 ))
    done
    echo "[INFO] Collected ${COLLECTED} existing SpecHLA results from spechla_work/"

    # Discover samples typed by any other tool
    SAMPLE_LIST_FILE="${BASE}/conf/spechla_pending.txt"
    : > "$SAMPLE_LIST_FILE"

    for TOOL in hlahd optitype arcashla; do
        for F in "${RESULTS_DIR}"/*/; do
            SAMPLE=$(basename "$F")
            RESULT="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
            [[ -f "$RESULT" ]] || continue
            # Skip if SpecHLA already done (non-empty result with >1 line)
            DONE="${RESULTS_DIR}/by_tool/spechla/${SAMPLE}_spechla.txt"
            if [[ -s "$DONE" ]] && [[ $(wc -l < "$DONE") -gt 1 ]]; then
                continue
            fi
            echo "$SAMPLE"
        done
    done | sort -u >> "$SAMPLE_LIST_FILE"

    N=$(wc -l < "$SAMPLE_LIST_FILE")

    if [[ "$N" -eq 0 ]]; then
        echo "[INFO] All samples already have SpecHLA results. Nothing to do."
        exit 0
    fi

    echo "[INFO] Samples to run: $N"
    cat "$SAMPLE_LIST_FILE"
    echo ""

    sbatch \
        --job-name=spechla_1kgp \
        --account="${PROJECT_ID}" \
        --partition=small \
        --time=24:00:00 \
        --cpus-per-task=8 \
        --mem=32G \
        --array="1-${N}" \
        --output="${LOGS_DIR}/spechla_%A_%a.out" \
        --error="${LOGS_DIR}/spechla_%A_%a.err" \
        "$0"

    echo "[INFO] Array job submitted. Monitor with: squeue -u \$USER"
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
WORK_BASE="${BASE}/spechla_work"
HLA_REGION="chr6:28000000-34000000"
SAMPLE_LIST_FILE="${BASE}/conf/spechla_pending.txt"

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

echo "=== SpecHLA: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Date: $(date)"

#-----------------------------------------------------------------------------
# Load modules + environment (from user's confirmed working setup)
#-----------------------------------------------------------------------------
module purge
module load gcc
module load intel-oneapi-mkl
module load samtools/1.21
module load bwa
module load bowtie2

export LD_LIBRARY_PATH=/appl/soft/bio/samtools/gcc_11.3.0/htslib/htslib-1.21/lib:/appl/spack/v018/install-tree/gcc-11.3.0/arpack-ng-3.8.0-mtifxr/lib64:${LD_LIBRARY_PATH:-}
export PATH="${SPEC_HOME}/spechla_env/bin:${SPEC_HOME}/script:$PATH"
export SPECHLA_HOME="${SPEC_HOME}"

echo "python: $(which python 2>/dev/null || echo 'not found')"
echo "SpecHLA: ${SPEC_HOME}/script/whole/SpecHLA.sh"

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
        2>"${LOGS_DIR}/spechla_stream_${SAMPLE}.log" \
        || { echo "ERROR: samtools stream failed"; cat "${LOGS_DIR}/spechla_stream_${SAMPLE}.log"; exit 1; }
    mv "${TMPBAM}.tmp" "$TMPBAM"

    echo "Name-sorting..."
    samtools sort -n -@ 4 "$TMPBAM" -o "$NAMESORT" \
        2>>"${LOGS_DIR}/spechla_stream_${SAMPLE}.log"
    rm -f "$TMPBAM"

    echo "Converting to FASTQ..."
    mkdir -p "$FASTQ_DIR"
    samtools fastq -@ 4 \
        -1 "$R1" -2 "$R2" \
        -0 /dev/null -s /dev/null \
        "$NAMESORT" \
        2>>"${LOGS_DIR}/spechla_stream_${SAMPLE}.log"
    rm -f "$NAMESORT"

    echo "[OK] FASTQ written: $(du -sh $R1 | cut -f1) / $(du -sh $R2 | cut -f1)"
fi

#-----------------------------------------------------------------------------
# Step 2: Create workdir; symlink FASTQs using absolute paths
#-----------------------------------------------------------------------------
WORKDIR="${WORK_BASE}/${SAMPLE}"
mkdir -p "${WORKDIR}/${SAMPLE}"

ln -sf "$(realpath "$R1")" "${WORKDIR}/R1.fastq.gz"
ln -sf "$(realpath "$R2")" "${WORKDIR}/R2.fastq.gz"

#-----------------------------------------------------------------------------
# Step 3: Run SpecHLA
#-----------------------------------------------------------------------------
echo "[Step 3] Running SpecHLA..."
cd "$WORKDIR"

bash "${SPEC_HOME}/script/whole/SpecHLA.sh" \
    -n "${SAMPLE}" \
    -1 R1.fastq.gz \
    -2 R2.fastq.gz \
    -o . \
    -j "${SLURM_CPUS_PER_TASK}" \
    -u 1

echo "SpecHLA finished."

#-----------------------------------------------------------------------------
# Step 4: Find hla.result.txt (SpecHLA may nest it)
#-----------------------------------------------------------------------------
RESULT_TXT=""
if [[ -f "${WORKDIR}/${SAMPLE}/hla.result.txt" ]]; then
    RESULT_TXT="${WORKDIR}/${SAMPLE}/hla.result.txt"
elif [[ -f "${WORKDIR}/${SAMPLE}/${SAMPLE}/hla.result.txt" ]]; then
    RESULT_TXT="${WORKDIR}/${SAMPLE}/${SAMPLE}/hla.result.txt"
else
    echo "ERROR: hla.result.txt not found under ${WORKDIR}/${SAMPLE}/"
    find "${WORKDIR}" -name "hla.result.txt" 2>/dev/null || true
    exit 1
fi

echo "[OK] Result: $RESULT_TXT"

# Validate it has allele data (more than just the header)
NLINES=$(wc -l < "$RESULT_TXT")
if [[ "$NLINES" -le 1 ]]; then
    echo "ERROR: hla.result.txt has only ${NLINES} line(s) — no allele calls produced"
    cat "$RESULT_TXT"
    exit 1
fi

#-----------------------------------------------------------------------------
# Step 5: Publish result
#-----------------------------------------------------------------------------
OUT_DIR="${RESULTS_DIR}/${SAMPLE}/spechla"
mkdir -p "$OUT_DIR"
cp "$RESULT_TXT" "${OUT_DIR}/${SAMPLE}_spechla.txt"
echo "[OK] Published: ${OUT_DIR}/${SAMPLE}_spechla.txt"

#-----------------------------------------------------------------------------
# Step 6: Symlink into by_tool layout for calibration
#-----------------------------------------------------------------------------
mkdir -p "${RESULTS_DIR}/by_tool/spechla"
SRC="${OUT_DIR}/${SAMPLE}_spechla.txt"
DST="${RESULTS_DIR}/by_tool/spechla/${SAMPLE}_spechla.txt"
[[ -e "$DST" ]] && rm -f "$DST"
ln -sf "$SRC" "$DST"
echo "[OK] Symlinked: ${DST}"

echo ""
echo "=== Done: ${SAMPLE} ==="
head -2 "${OUT_DIR}/${SAMPLE}_spechla.txt"
