#!/bin/bash
# =============================================================================
# run_wes_calibration_puhti.sh
# WES HLA calibration on CSC Puhti — 1000 Genomes Phase 3 exome BAMs
#
# Uses 1KGP Phase 3 WES BAMs (hg19/GRCh37) from EBI FTP.
# Ground truth: Gourraud et al. 2014 (same as WGS calibration).
# Produces: conf/tool_weights_wes_v1.json, conf/tool_accuracy_wes_v1.tsv
#
# Usage:
#   bash scripts/run_wes_calibration_puhti.sh [OPTIONS]
#
# Options:
#   --project PROJECT_ID   CSC project account (default: $SLURM_JOB_ACCOUNT or project_2008084)
#   --sample-list FILE     Sample list (default: conf/wes_samples_50.txt)
#   --batch-size N         Process next N unextracted samples per run (default: 10)
#   --tools TOOLS          Comma-separated tools (default: hlahd,spechla,arcashla,optitype)
#   --genes GENES          Comma-separated genes (default: A,B,C,DRB1,DQB1)
#   --skip-extract         Skip Phase 0 (FASTQs already present)
#   --skip-typing          Skip Phase 1 (results already present)
#   --calibrate-only       Skip to Phase 2 (calibrate from existing results)
#   --status               Show current progress and exit
#   --dry-run              Print commands without submitting
#   -h, --help             Show this help
#
# Batching: each run extracts and types the next --batch-size samples without FASTQs.
# Re-run after each batch completes. Calibration accumulates all results.
# Example (5 runs of 10):
#   bash scripts/run_wes_calibration_puhti.sh --project project_2008084 --batch-size 10
#
# Phase 0: SLURM array — stream HLA region from EBI WES BAM → paired FASTQs
# Phase 1: SLURM job  — run Nextflow typing batch (seq_type=wes)
# Phase 1b: SLURM job — collect results into by_tool/ layout
# Phase 2: SLURM job  — calibrate weights + compare voting strategies
# =============================================================================
set -euo pipefail

#-----------------------------------------------------------------------------
# Defaults
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="$(dirname "$SCRIPT_DIR")"   # hla_typing_pipeline/

PROJECT_ID="${SLURM_JOB_ACCOUNT:-project_2008084}"
TOOLS="hlahd,spechla,arcashla,optitype,seq2hla,kourami,polysolver"
GENES="A,B,C,DRB1,DQB1"
RESOLUTION="2-field"
SAMPLE_LIST_DEFAULT="${INSTALL_DIR}/conf/wes_samples_50.txt"
SAMPLE_LIST=""
BATCH_SIZE=10   # samples per run; 0 = all
SKIP_EXTRACT=false
SKIP_TYPING=false
CALIBRATE_ONLY=false
STATUS_ONLY=false
DRY_RUN=false

# EBI FTP paths
EBI_EXOME_INDEX="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.exome.sequence.index"
HLA_REGION="6:28000000-34000000"   # hg19/GRCh37 ENSEMBL (no chr prefix)

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)   PROJECT_ID="$2";       shift 2 ;;
        --sample-list) SAMPLE_LIST="$2";    shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";     shift 2 ;;
        --tools)     TOOLS="$2";            shift 2 ;;
        --genes)     GENES="$2";            shift 2 ;;
        --skip-extract)   SKIP_EXTRACT=true;   shift ;;
        --skip-typing)    SKIP_TYPING=true;    shift ;;
        --calibrate-only) CALIBRATE_ONLY=true; SKIP_EXTRACT=true; SKIP_TYPING=true; shift ;;
        --status)    STATUS_ONLY=true;      shift ;;
        --dry-run)   DRY_RUN=true;          shift ;;
        -h|--help)
            sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$SAMPLE_LIST" ]] && SAMPLE_LIST="$SAMPLE_LIST_DEFAULT"

#-----------------------------------------------------------------------------
# Derived paths
#-----------------------------------------------------------------------------
SCRATCH_BASE="/scratch/${PROJECT_ID}/hla_calibration/wes"
FASTQ_DIR="${SCRATCH_BASE}/fastqs"
RESULTS_DIR="${SCRATCH_BASE}/results"
BY_TOOL_DIR="${RESULTS_DIR}/by_tool"
INDEX_DIR="${SCRATCH_BASE}/index"
LOGS_DIR="${SCRATCH_BASE}/logs"
GT_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_hla_gt.tsv"
POP_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_populations.tsv"
WEIGHTS_OUT="${INSTALL_DIR}/conf/tool_weights_wes_v1.json"
TABLE_OUT="${INSTALL_DIR}/conf/tool_accuracy_wes_v1.tsv"

#-----------------------------------------------------------------------------
# Build effective sample list for this batch
# - Strip comments/blanks from master list
# - Filter to samples whose FASTQs do NOT yet exist (--batch-size 0 = all pending)
# - Take first BATCH_SIZE of those
#-----------------------------------------------------------------------------
EFFECTIVE_LIST="${SCRATCH_BASE}/wes_samples_effective.txt"
ALL_SAMPLES_LIST="${SCRATCH_BASE}/wes_samples_all.txt"

mkdir -p "${SCRATCH_BASE}" "${FASTQ_DIR}" "${RESULTS_DIR}" "${INDEX_DIR}" "${LOGS_DIR}"

# All samples (comments/blanks stripped)
grep -v '^#' "${SAMPLE_LIST}" | grep -v '^[[:space:]]*$' > "${ALL_SAMPLES_LIST}"
N_ALL=$(wc -l < "${ALL_SAMPLES_LIST}")

# Identify pending samples (no R1 FASTQ yet)
PENDING_LIST="${SCRATCH_BASE}/wes_samples_pending.txt"
> "${PENDING_LIST}"
while IFS= read -r S; do
    [[ ! -f "${FASTQ_DIR}/${S}_R1.fastq.gz" ]] && echo "$S" >> "${PENDING_LIST}"
done < "${ALL_SAMPLES_LIST}"
N_PENDING=$(wc -l < "${PENDING_LIST}")

# Apply batch size
if [[ "${BATCH_SIZE}" -gt 0 && "${N_PENDING}" -gt "${BATCH_SIZE}" ]]; then
    head -n "${BATCH_SIZE}" "${PENDING_LIST}" > "${EFFECTIVE_LIST}"
else
    cp "${PENDING_LIST}" "${EFFECTIVE_LIST}"
fi

N_TOTAL=$(wc -l < "${EFFECTIVE_LIST}")
N_DONE=$(( N_ALL - N_PENDING ))

#-----------------------------------------------------------------------------
# --status: show progress and exit
#-----------------------------------------------------------------------------
if [[ "$STATUS_ONLY" == "true" ]]; then
    echo "=== WES Calibration Status ==="
    echo "Project:     ${PROJECT_ID}"
    echo "Master list: ${N_ALL} samples total"
    echo ""
    echo "Phase 0 — FASTQs extracted:"
    N_FQ=$(find "${FASTQ_DIR}" -name "*_R1.fastq.gz" 2>/dev/null | wc -l || echo 0)
    echo "  ${N_FQ} / ${N_ALL} samples done   (${N_PENDING} pending)"
    echo ""
    echo "Phase 1 — Typing results:"
    for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
        N_RES=$(find "${BY_TOOL_DIR}/${TOOL}/" -name "*_${TOOL}.txt" 2>/dev/null | wc -l || echo 0)
        printf "  %-12s %3d samples\n" "$TOOL" "$N_RES"
    done
    echo ""
    echo "Phase 2 — Calibration weights:"
    [[ -f "$WEIGHTS_OUT" ]] && echo "  $WEIGHTS_OUT ($(python3 -c "import json; d=json.load(open('$WEIGHTS_OUT')); print(f\"n={d.get('n_samples','?')}\")" 2>/dev/null || echo "parse error"))" || echo "  NOT YET PRODUCED"
    exit 0
fi

# Nothing to extract?
if [[ "$SKIP_EXTRACT" == "false" && "$CALIBRATE_ONLY" == "false" && "${N_TOTAL}" -eq 0 ]]; then
    echo "[INFO] All ${N_ALL} samples already have FASTQs. Use --calibrate-only to run calibration."
    exit 0
fi

echo "==================================================================="
echo " WES HLA Calibration — CSC Puhti"
echo "==================================================================="
echo " Project:      ${PROJECT_ID}"
echo " Samples:      ${N_TOTAL} this batch  (${N_DONE}/${N_ALL} total done; ${N_PENDING} pending)"
echo " Tools:        ${TOOLS}"
echo " Genes:        ${GENES}"
echo " FASTQs:       ${FASTQ_DIR}"
echo " Results:      ${RESULTS_DIR}"
echo " Weights out:  ${WEIGHTS_OUT}"
echo " Dry-run:      ${DRY_RUN}"
echo "==================================================================="
echo ""

#-----------------------------------------------------------------------------
# Download exome sequence index (once; cached)
# Format: tab-sep, columns include SAMPLE_NAME and FASTQ_FILE (relative path)
# We use it to build exact BAM URLs per sample.
#-----------------------------------------------------------------------------
INDEX_FILE="${INDEX_DIR}/phase3_exome.sequence.index"
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "[INFO] Downloading 1KGP Phase 3 exome sequence index..."
    if [[ "$DRY_RUN" == "false" ]]; then
        wget -q -O "${INDEX_FILE}" "${EBI_EXOME_INDEX}" \
            || curl -s -o "${INDEX_FILE}" "${EBI_EXOME_INDEX}" \
            || { echo "[ERROR] Could not download exome sequence index from EBI FTP"; exit 1; }
        echo "[OK] Index saved: ${INDEX_FILE} ($(wc -l < "$INDEX_FILE") entries)"
    else
        echo "[DRY-RUN] wget -q -O ${INDEX_FILE} ${EBI_EXOME_INDEX}"
    fi
else
    echo "[INFO] Using cached index: ${INDEX_FILE}"
fi

#-----------------------------------------------------------------------------
# Build per-sample BAM URL map: SAMPLE -> full BAM FTP URL
# Columns in the index (header line starts with STUDY_ID):
#   col 1 = STUDY_ID, col 2 = SAMPLE_NAME, col 28 = FILE (relative path)
# The relative path starts with "data/" so FTP base = ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/
#-----------------------------------------------------------------------------
URL_MAP="${INDEX_DIR}/sample_bam_urls.tsv"
if [[ ! -f "$URL_MAP" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Building sample→BAM URL map..."
    python3 - << PYEOF
import sys, re

index_file = "${INDEX_FILE}"
url_map_out = "${URL_MAP}"
ftp_base = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/"

sample_urls = {}

with open(index_file) as fh:
    for line in fh:
        line = line.rstrip('\n')
        if not line or line.startswith('STUDY_ID'):
            continue
        parts = line.split('\t')
        if len(parts) < 28:
            continue
        sample = parts[9].strip()  # SAMPLE_NAME column (0-based index 9)
        filepath = parts[0].strip() # FILE column (0-based index 0 in some versions)
        # The index has FILE as first column starting with "data/"
        # Detect which column has the BAM path
        bam_col = None
        for i, p in enumerate(parts):
            if re.search(r'\.bam$', p) and 'exome' in p:
                bam_col = i
                break
        if bam_col is None:
            continue
        bam_path = parts[bam_col].strip()
        if not bam_path:
            continue
        # Strip leading slash or "vol1/ftp/phase3/" if already in path
        bam_path = re.sub(r'^/vol1/ftp/phase3/', '', bam_path)
        bam_path = re.sub(r'^phase3/', '', bam_path)
        url = ftp_base + bam_path.lstrip('/')
        # Keep the first BAM per sample (typically the primary alignment)
        if sample not in sample_urls:
            sample_urls[sample] = url

with open(url_map_out, 'w') as fout:
    for sample, url in sorted(sample_urls.items()):
        fout.write(f"{sample}\t{url}\n")

print(f"[OK] URL map written: {url_map_out} ({len(sample_urls)} samples)")
PYEOF
fi

#-----------------------------------------------------------------------------
# Phase 0: SLURM array — extract HLA FASTQs from WES BAMs
#-----------------------------------------------------------------------------
PHASE0_DEP=""

if [[ "$SKIP_EXTRACT" == "false" ]]; then
    echo "=== Phase 0: Submitting FASTQ extraction array (${N_TOTAL} tasks) ==="

    EXTRACT_SCRIPT="${SCRATCH_BASE}/phase0_extract.sh"
    cat > "$EXTRACT_SCRIPT" << 'EXTRACTEOF'
#!/bin/bash
#SBATCH --job-name=wes_extract
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=LOGS_PLACEHOLDER/extract_%A_%a.out
#SBATCH --error=LOGS_PLACEHOLDER/extract_%A_%a.err

set -euo pipefail
module purge
module load samtools 2>/dev/null || true

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "EFFECTIVE_LIST_PLACEHOLDER")
[[ -z "$SAMPLE" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

FASTQ_DIR="FASTQ_DIR_PLACEHOLDER"
URL_MAP="URL_MAP_PLACEHOLDER"
HLA_REGION="HLA_REGION_PLACEHOLDER"

R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Date:   $(date)"

# Skip if already done
if [[ -f "$R1" ]] && [[ -f "$R2" ]]; then
    echo "[SKIP] FASTQs already present for ${SAMPLE}"
    exit 0
fi

# Look up BAM URL
BAM_URL=$(grep -P "^${SAMPLE}\t" "$URL_MAP" | cut -f2 || true)
if [[ -z "$BAM_URL" ]]; then
    echo "[ERROR] No BAM URL found for sample ${SAMPLE} in URL map" >&2
    exit 1
fi
echo "[INFO] BAM URL: ${BAM_URL}"

# Stream HLA region from remote BAM → paired FASTQs
# samtools view requires the remote .bai index to be present at ${BAM_URL}.bai
mkdir -p "${FASTQ_DIR}"
TMP_BAM="${FASTQ_DIR}/${SAMPLE}_hla_tmp.bam"

samtools view -b -h -o "${TMP_BAM}" "${BAM_URL}" "${HLA_REGION}"
echo "[INFO] HLA reads: $(samtools view -c "$TMP_BAM")"

# Name-sort and convert to FASTQ (paired; singletons discarded)
samtools sort -n -@ 3 -m 2G -o "${TMP_BAM}.nsort.bam" "${TMP_BAM}"
samtools fastq \
    -1 "${R1}" \
    -2 "${R2}" \
    -s /dev/null \
    "${TMP_BAM}.nsort.bam"

rm -f "${TMP_BAM}" "${TMP_BAM}.nsort.bam"

echo "[OK] FASTQs:"
echo "  R1: ${R1} ($(wc -l < <(zcat "${R1}") / 4) reads)"
echo "  R2: ${R2}"
EXTRACTEOF

    # Substitute placeholders (heredoc can't expand variables inside 'EXTRACTEOF')
    sed -i \
        -e "s|LOGS_PLACEHOLDER|${LOGS_DIR}|g" \
        -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
        -e "s|FASTQ_DIR_PLACEHOLDER|${FASTQ_DIR}|g" \
        -e "s|URL_MAP_PLACEHOLDER|${URL_MAP}|g" \
        -e "s|HLA_REGION_PLACEHOLDER|${HLA_REGION}|g" \
        "$EXTRACT_SCRIPT"

    CONCURRENT=$([[ "${BATCH_SIZE}" -gt 0 ]] && echo "${BATCH_SIZE}" || echo "10")
    SBATCH_EXTRACT="sbatch --parsable \
        --account=${PROJECT_ID} \
        --array=1-${N_TOTAL}%${CONCURRENT} \
        ${EXTRACT_SCRIPT}"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_EXTRACT}"
        EXTRACT_JOB="DRY_RUN_JOB"
    else
        EXTRACT_JOB=$(eval "$SBATCH_EXTRACT")
        echo "[INFO] Phase 0 extraction array submitted: ${EXTRACT_JOB}"
        PHASE0_DEP="--dependency=afterok:${EXTRACT_JOB}"
    fi
else
    echo "[SKIP] Phase 0 — using existing FASTQs"
fi

#-----------------------------------------------------------------------------
# Phase 1: Build samplesheet + submit Nextflow batch typing job
#-----------------------------------------------------------------------------
PHASE1_DEP=""
TYPING_JOB=""

if [[ "$SKIP_TYPING" == "false" ]]; then
    echo ""
    echo "=== Phase 1: Submitting Nextflow typing job ==="

    # Generate samplesheet from extracted FASTQs
    SAMPLESHEET="${SCRATCH_BASE}/wes_samplesheet.csv"
    cat > "${SCRATCH_BASE}/gen_samplesheet.sh" << GENEOF
#!/bin/bash
echo "sample_id,fastq_1,fastq_2" > "${SAMPLESHEET}"
while IFS= read -r SAMPLE; do
    R1="${FASTQ_DIR}/\${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/\${SAMPLE}_R2.fastq.gz"
    if [[ -f "\$R1" ]] && [[ -f "\$R2" ]]; then
        echo "\${SAMPLE},\${R1},\${R2}" >> "${SAMPLESHEET}"
    else
        echo "[WARN] Missing FASTQs for \${SAMPLE} — skipped" >&2
    fi
done < "${EFFECTIVE_LIST}"
echo "[OK] Samplesheet: ${SAMPLESHEET} (\$(tail -n +2 ${SAMPLESHEET} | wc -l) samples)"
GENEOF

    # NOTE: samplesheet is generated inside the Phase 1 SLURM job (after Phase 0 FASTQs are ready)
    echo "[INFO] Samplesheet will be generated in Phase 1 job: ${SAMPLESHEET}"

    TYPING_SCRIPT="${SCRATCH_BASE}/phase1_typing.sh"
    cat > "$TYPING_SCRIPT" << TYPINGEOF
#!/bin/bash
#SBATCH --job-name=wes_typing
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=72:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=40
#SBATCH --output=${LOGS_DIR}/typing_%j.out
#SBATCH --error=${LOGS_DIR}/typing_%j.err

set -euo pipefail
module purge
module load nextflow 2>/dev/null || module load nextflow/23.10.0
module load apptainer 2>/dev/null || module load singularity

echo "=== Phase 1: HLA Typing (WES) ==="
echo "Date: \$(date)"
echo "Samplesheet: ${SAMPLESHEET}"

cd "${INSTALL_DIR}"

# Regenerate samplesheet in case Phase 0 completed after original run
bash "${SCRATCH_BASE}/gen_samplesheet.sh"

N_SAMPLES_SHEET=\$(tail -n +2 "${SAMPLESHEET}" | wc -l)
echo "Samples to type: \${N_SAMPLES_SHEET}"

nextflow run main.nf \\
    --input_samplesheet "${SAMPLESHEET}" \\
    --seq_type wes \\
    --reference hg38 \\
    --tools "${TOOLS}" \\
    --weighting equal \\
    --skip_qc true \\
    --outdir "${RESULTS_DIR}" \\
    -c conf/puhti.config \\
    -profile apptainer \\
    -resume \\
    -with-trace "${RESULTS_DIR}/pipeline_info/trace_wes.txt" \\
    -with-report "${RESULTS_DIR}/pipeline_info/report_wes.html" \\
    -work-dir "${SCRATCH_BASE}/work"

echo ""
echo "[OK] Phase 1 complete — results in ${RESULTS_DIR}"
TYPINGEOF

    SBATCH_TYPING="sbatch --parsable \
        ${PHASE0_DEP} \
        ${TYPING_SCRIPT}"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_TYPING}"
        TYPING_JOB="DRY_RUN_JOB"
    else
        TYPING_JOB=$(eval "$SBATCH_TYPING")
        echo "[INFO] Phase 1 typing job submitted: ${TYPING_JOB}"
        PHASE1_DEP="--dependency=afterok:${TYPING_JOB}"
    fi
else
    echo "[SKIP] Phase 1 — using existing typing results"
fi

#-----------------------------------------------------------------------------
# Phase 1b: Collect results into by_tool/ layout
#-----------------------------------------------------------------------------
COLLECT_SCRIPT="${SCRATCH_BASE}/phase1b_collect.sh"
cat > "$COLLECT_SCRIPT" << COLLECTEOF
#!/bin/bash
#SBATCH --job-name=wes_collect
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/collect_%j.out
#SBATCH --error=${LOGS_DIR}/collect_%j.err

set -euo pipefail
echo "=== Phase 1b: Collecting WES results ==="
echo "Date: \$(date)"

for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    mkdir -p "${BY_TOOL_DIR}/\${TOOL}"
done

TOTAL=0
while IFS= read -r SAMPLE; do
    for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
        SRC="${RESULTS_DIR}/\${SAMPLE}/\${TOOL}/\${SAMPLE}_\${TOOL}.txt"
        DST="${BY_TOOL_DIR}/\${TOOL}/\${SAMPLE}_\${TOOL}.txt"
        if [[ -f "\$SRC" ]] && [[ ! -e "\$DST" ]]; then
            ln -sf "\$SRC" "\$DST"
            TOTAL=\$((TOTAL+1))
        fi
    done
done < "${EFFECTIVE_LIST}"

echo ""
echo "Symlinks created: \$TOTAL"
echo ""
echo "Results by tool:"
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    N=\$(find "${BY_TOOL_DIR}/\${TOOL}/" -maxdepth 1 -name "*_\${TOOL}.txt" 2>/dev/null | wc -l)
    printf "  %-12s %3d samples\n" "\$TOOL" "\$N"
done
COLLECTEOF

SBATCH_COLLECT="sbatch --parsable \
    ${PHASE1_DEP} \
    ${COLLECT_SCRIPT}"

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] ${SBATCH_COLLECT}"
    COLLECT_JOB="DRY_RUN_JOB"
else
    COLLECT_JOB=$(eval "$SBATCH_COLLECT")
    echo "[INFO] Phase 1b collect job submitted: ${COLLECT_JOB}"
fi

#-----------------------------------------------------------------------------
# Phase 2: Calibrate weights + compare voting strategies
#-----------------------------------------------------------------------------
CALIBRATE_SCRIPT="${SCRATCH_BASE}/phase2_calibrate.sh"
COMPARE_OUT="${SCRATCH_BASE}/strategy_comparison_wes_calibrated.tsv"
COMPARE_RC_OUT="${SCRATCH_BASE}/strategy_comparison_wes_rc.tsv"

cat > "$CALIBRATE_SCRIPT" << CALEOF
#!/bin/bash
#SBATCH --job-name=wes_calibrate
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/calibrate_%j.out
#SBATCH --error=${LOGS_DIR}/calibrate_%j.err

set -euo pipefail
echo "=== Phase 2: WES HLA Tool Weight Calibration ==="
echo "Date: \$(date)"
echo ""

module purge
module load python-data 2>/dev/null || module load python/3.9 2>/dev/null || module load python

# Download ground truth if not present
if [[ ! -f "${GT_FILE}" ]]; then
    echo "[INFO] Downloading 1KGP Phase 1 ground truth..."
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" download-gt \
        --output "${GT_FILE}"
else
    echo "[INFO] Ground truth: ${GT_FILE}"
fi

# Download population file if not present
if [[ ! -f "${POP_FILE}" ]]; then
    echo "[INFO] Downloading 1KGP population panel..."
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" generate-population-file \
        --output "${POP_FILE}" \
        && echo "[OK] Population file: ${POP_FILE}" \
        || echo "[WARN] Could not download population file — continuing without stratification"
else
    echo "[INFO] Population file: ${POP_FILE}"
fi

# Calibrate
echo ""
echo "[INFO] Running WES calibration (genes: ${GENES})..."
POP_FILE_ARG=""
[[ -f "${POP_FILE}" ]] && POP_FILE_ARG="--population-file ${POP_FILE}"

# shellcheck disable=SC2086
python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" calibrate \\
    --ground-truth "${GT_FILE}" \\
    --results-dir "${BY_TOOL_DIR}" \\
    --data-type wes \\
    --genes "${GENES}" \\
    --resolution "${RESOLUTION}" \\
    --output-weights "${WEIGHTS_OUT}" \\
    --output-table  "${TABLE_OUT}" \\
    \$POP_FILE_ARG

echo ""
echo "[OK] Calibration outputs:"
echo "     Weights: ${WEIGHTS_OUT}"
echo "     Table:   ${TABLE_OUT}"

# Print weight summary
python3 - << PYEOF
import json, sys
try:
    with open('${WEIGHTS_OUT}') as f:
        w = json.load(f)
    n = w.get('n_samples', '?')
    print(f'\nWES calibrated weights (n={n} samples):')
    genes = list(w.get('genes', {}).keys())
    tools = sorted({t for g in w['genes'].values() for t in g.keys()})
    header = f"{'Gene':<8}" + "".join(f"  {t:<12}" for t in tools)
    print(header)
    print('-' * len(header))
    for gene in genes:
        row = f"{gene:<8}"
        for t in tools:
            row += f"  {w['genes'][gene].get(t, 0.0):<12.3f}"
        print(row)
except Exception as e:
    print(f'[WARN] Could not print weight summary: {e}', file=sys.stderr)
PYEOF

# Compare voting strategies
echo ""
echo "[INFO] Comparing strategies: equal vs calibrated (Wilcoxon test)..."
if [[ -f "${WEIGHTS_OUT}" ]]; then
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" compare-strategies \\
        --ground-truth "${GT_FILE}" \\
        --results-dir "${BY_TOOL_DIR}" \\
        --genes "${GENES}" \\
        --resolution "${RESOLUTION}" \\
        --weights-file "${WEIGHTS_OUT}" \\
        --ref-mode equal \\
        --test-mode calibrated \\
        --output "${COMPARE_OUT}" \\
        && echo "[OK] Strategy comparison: ${COMPARE_OUT}" \\
        || echo "[WARN] Strategy comparison failed (scipy may not be available)"

    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" compare-strategies \\
        --ground-truth "${GT_FILE}" \\
        --results-dir "${BY_TOOL_DIR}" \\
        --genes "${GENES}" \\
        --resolution "${RESOLUTION}" \\
        --ref-mode equal \\
        --test-mode read_confidence \\
        --output "${COMPARE_RC_OUT}" \\
        && echo "[OK] RC comparison: ${COMPARE_RC_OUT}" \\
        || echo "[WARN] RC comparison failed"
else
    echo "[WARN] Weights file not found — skipping strategy comparison"
fi

echo ""
echo "================================================================="
echo " WES Calibration Complete"
echo "================================================================="
echo " Weights:      ${WEIGHTS_OUT}"
echo " Accuracy:     ${TABLE_OUT}"
[[ -f "${COMPARE_OUT}" ]]    && echo " Strat(cal):   ${COMPARE_OUT}"
[[ -f "${COMPARE_RC_OUT}" ]] && echo " Strat(rc):    ${COMPARE_RC_OUT}"
echo ""
echo " Next: push updated weights to GitHub and rerun pipeline"
echo "   --weighting calibrated --weights_file conf/tool_weights_wes_v1.json"
echo "================================================================="
CALEOF

SBATCH_CAL="sbatch --parsable \
    --dependency=afterok:${COLLECT_JOB} \
    ${CALIBRATE_SCRIPT}"

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] ${SBATCH_CAL}"
else
    CALIBRATE_JOB=$(eval "$SBATCH_CAL")
    echo "[INFO] Phase 2 calibration job submitted: ${CALIBRATE_JOB}"
fi

echo ""
echo "==================================================================="
echo " All jobs submitted."
echo ""
echo " Monitor: squeue -u \$USER"
echo " Status:  bash scripts/run_wes_calibration_puhti.sh --status"
echo " Logs:    ${LOGS_DIR}/"
echo ""
echo " Expected outputs:"
echo "   ${WEIGHTS_OUT}"
echo "   ${TABLE_OUT}"
echo "==================================================================="
