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
BAM_DIR="${SCRATCH_BASE}/bams"
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

mkdir -p "${SCRATCH_BASE}" "${BAM_DIR}" "${RESULTS_DIR}" "${INDEX_DIR}" "${LOGS_DIR}"

# All samples (comments/blanks stripped)
grep -v '^#' "${SAMPLE_LIST}" | grep -v '^[[:space:]]*$' > "${ALL_SAMPLES_LIST}"
N_ALL=$(wc -l < "${ALL_SAMPLES_LIST}")

# Identify pending samples (no HLA BAM yet)
PENDING_LIST="${SCRATCH_BASE}/wes_samples_pending.txt"
> "${PENDING_LIST}"
while IFS= read -r S; do
    [[ ! -f "${BAM_DIR}/${S}_hla.bam" ]] && echo "$S" >> "${PENDING_LIST}"
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
# GT filter: only keep samples present in the Gourraud 2014 ground truth
# Downloads GT if not present (needed for calibration anyway).
#-----------------------------------------------------------------------------
if [[ ! -f "${GT_FILE}" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Downloading 1KGP ground truth (Gourraud 2014)..."
    module load python-data 2>/dev/null || module load python/3.9 2>/dev/null || true
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" download-gt \
        --output "${GT_FILE}" \
        && echo "[OK] GT saved: ${GT_FILE}" \
        || { echo "[ERROR] Could not download GT — aborting"; exit 1; }
fi

if [[ -f "${GT_FILE}" ]]; then
    FILTERED_LIST="${SCRATCH_BASE}/wes_samples_gt_confirmed.txt"
    python3 - << GTEOF
import sys
gt_samples = set()
with open("${GT_FILE}") as f:
    for line in f:
        s = line.split('\t')[0].strip() if line.strip() else ''
        if s and s != 'Sample':
            gt_samples.add(s)
kept, skipped = [], []
with open("${EFFECTIVE_LIST}") as f:
    for line in f:
        s = line.strip()
        if s:
            (kept if s in gt_samples else skipped).append(s)
with open("${FILTERED_LIST}", "w") as f:
    f.write("\n".join(kept) + ("\n" if kept else ""))
if skipped:
    print(f"[WARN] {len(skipped)} sample(s) not in GT, excluded: {', '.join(skipped)}", file=sys.stderr)
print(f"[INFO] GT-confirmed: {len(kept)} / {len(kept)+len(skipped)} samples in this batch")
GTEOF
    EFFECTIVE_LIST="${FILTERED_LIST}"
    N_TOTAL=$(wc -l < "${EFFECTIVE_LIST}")
fi

#-----------------------------------------------------------------------------
# --status: show progress and exit
#-----------------------------------------------------------------------------
if [[ "$STATUS_ONLY" == "true" ]]; then
    echo "=== WES Calibration Status ==="
    echo "Project:     ${PROJECT_ID}"
    echo "Master list: ${N_ALL} samples total"
    echo ""
    echo "Phase 0 — HLA BAMs downloaded:"
    N_DONE=$(( N_ALL - N_PENDING ))
    echo "  ${N_DONE} / ${N_ALL} samples done   (${N_PENDING} pending)"
    echo ""
    echo "Phase 1 — Typing results:"
    for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
        N_RES=$(find "${BY_TOOL_DIR}/${TOOL}/" -name "*_${TOOL}.txt" 2>/dev/null | wc -l || true)
        N_RES=$(echo "${N_RES}" | tr -d '[:space:]')
        N_RES=${N_RES:-0}
        printf "  %-12s %3d samples\n" "$TOOL" "$N_RES"
    done
    echo ""
    echo "Phase 2 — Calibration weights:"
    [[ -f "$WEIGHTS_OUT" ]] && echo "  $WEIGHTS_OUT ($(python3 -c "import json; d=json.load(open('$WEIGHTS_OUT')); print(f\"n={d.get('n_samples','?')}\")" 2>/dev/null || echo "parse error"))" || echo "  NOT YET PRODUCED"
    exit 0
fi

# Nothing to extract?
if [[ "$SKIP_EXTRACT" == "false" && "$CALIBRATE_ONLY" == "false" && "${N_TOTAL}" -eq 0 ]]; then
    echo "[INFO] All ${N_ALL} samples already have HLA BAMs. Use --calibrate-only to run calibration."
    exit 0
fi

echo "==================================================================="
echo " WES HLA Calibration — CSC Puhti"
echo "==================================================================="
echo " Project:      ${PROJECT_ID}"
echo " Samples:      ${N_TOTAL} this batch  (${N_DONE}/${N_ALL} total done; ${N_PENDING} pending)"
echo " Tools:        ${TOOLS}"
echo " Genes:        ${GENES}"
echo " BAMs:         ${BAM_DIR}"
echo " Results:      ${RESULTS_DIR}"
echo " Weights out:  ${WEIGHTS_OUT}"
echo " Dry-run:      ${DRY_RUN}"
echo "==================================================================="
echo ""

#-----------------------------------------------------------------------------
# Build per-sample BAM URL map: SAMPLE -> full BAM FTP URL
#
# Root cause of earlier failure: 20130502.phase3.exome.sequence.index lists
# raw FASTQ files, not BAM alignments. BAM URLs must be constructed from the
# 1KGP population panel which maps sample_id -> population code (e.g. YRI).
#
# BAM URL pattern:
#   ftp://.../phase3/data/{SAMPLE}/exome_alignment/
#     {SAMPLE}.mapped.ILLUMINA.bwa.{POP}.exome.{DATE}.bam
# Dates tried in order: 20121211, 20130415, 20120522
#-----------------------------------------------------------------------------
PANEL_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
PANEL_FILE="${INDEX_DIR}/1kgp_panel.tsv"
URL_MAP="${INDEX_DIR}/sample_bam_urls.tsv"

# Download population panel (cached)
if [[ ! -f "$PANEL_FILE" ]]; then
    echo "[INFO] Downloading 1KGP population panel..."
    if [[ "$DRY_RUN" == "false" ]]; then
        wget -q -O "${PANEL_FILE}" "${PANEL_URL}" \
            || curl -s -o "${PANEL_FILE}" "${PANEL_URL}" \
            || { echo "[ERROR] Could not download population panel from EBI"; exit 1; }
        echo "[OK] Panel saved: ${PANEL_FILE} ($(wc -l < "$PANEL_FILE") entries)"
    else
        echo "[DRY-RUN] wget -q -O ${PANEL_FILE} ${PANEL_URL}"
    fi
else
    echo "[INFO] Using cached population panel: ${PANEL_FILE}"
fi

# (Re)build URL map whenever the panel changes or map is missing
# Delete stale map built from the wrong (FASTQ) index
if [[ -f "$URL_MAP" ]]; then
    N_MAPPED=$(wc -l < "$URL_MAP")
    if [[ "$N_MAPPED" -eq 0 ]]; then
        echo "[INFO] Stale empty URL map found — rebuilding from population panel"
        rm -f "$URL_MAP"
    fi
fi

if [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Building sample→BAM URL map from population panel..."
    python3 - << PYEOF
import sys, os

panel_file  = "${PANEL_FILE}"
sample_list = "${ALL_SAMPLES_LIST}"
url_map_out = "${URL_MAP}"
ftp_base    = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
# Dates to try in order (most samples use 20121211)
DATES = ["20121211", "20130415", "20120522"]

# Load population panel: sample_id \t pop \t super_pop \t gender
pop_map = {}
with open(panel_file) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("sample"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            pop_map[parts[0]] = parts[1]   # sample_id -> pop (e.g. YRI)

# Load sample list
with open(sample_list) as fh:
    samples = [l.strip() for l in fh if l.strip() and not l.startswith("#")]

found, missing = 0, []
with open(url_map_out, "w") as fout:
    for sample in samples:
        pop = pop_map.get(sample)
        if not pop:
            missing.append(sample)
            continue
        # Use the first date as primary URL; the extraction script tries all dates
        url = f"{ftp_base}/{sample}/exome_alignment/{sample}.mapped.ILLUMINA.bwa.{pop}.exome.{DATES[0]}.bam"
        # Also write alt-date URLs as fallbacks (tab-separated after main URL)
        alts = [f"{ftp_base}/{sample}/exome_alignment/{sample}.mapped.ILLUMINA.bwa.{pop}.exome.{d}.bam"
                for d in DATES[1:]]
        fout.write(f"{sample}\t{url}\t{'|'.join(alts)}\n")
        found += 1

print(f"[OK] URL map: {found} samples mapped, {len(missing)} missing from panel")
if missing:
    print(f"[WARN] Not in panel: {', '.join(missing)}", file=sys.stderr)
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

BAM_DIR="BAM_DIR_PLACEHOLDER"
URL_MAP="URL_MAP_PLACEHOLDER"
HLA_REGION="HLA_REGION_PLACEHOLDER"

BAM_OUT="${BAM_DIR}/${SAMPLE}_hla.bam"

echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Date:   $(date)"

# Skip if already done (BAM + index both present)
if [[ -f "${BAM_OUT}" ]] && [[ -f "${BAM_OUT}.bai" ]]; then
    echo "[SKIP] HLA BAM already present for ${SAMPLE}"
    exit 0
fi

# Look up BAM URL (col 2 = primary; col 3 = pipe-separated alt-date fallbacks)
LINE=$(grep -P "^${SAMPLE}\t" "$URL_MAP" || true)
if [[ -z "$LINE" ]]; then
    echo "[ERROR] No BAM URL found for sample ${SAMPLE} in URL map" >&2
    exit 1
fi
PRIMARY_URL=$(echo "$LINE" | cut -f2)
ALT_URLS=$(echo "$LINE" | cut -f3 | tr '|' ' ')

# Try primary URL, then alt-date fallbacks — stream HLA region directly as BAM
mkdir -p "${BAM_DIR}"
TMP_BAM="${BAM_DIR}/${SAMPLE}_hla_tmp.bam"
BAM_URL=""
for URL in ${PRIMARY_URL} ${ALT_URLS}; do
    echo "[INFO] Trying: ${URL}"
    if samtools view -b -h -o "${TMP_BAM}" "${URL}" "${HLA_REGION}" 2>/dev/null; then
        NREADS=$(samtools view -c "${TMP_BAM}" 2>/dev/null || echo 0)
        if [[ "${NREADS}" -gt 0 ]]; then
            BAM_URL="${URL}"
            echo "[OK] HLA reads: ${NREADS} from ${URL}"
            break
        else
            echo "[WARN] 0 reads from ${URL} — trying next"
            rm -f "${TMP_BAM}"
        fi
    else
        echo "[WARN] samtools failed for ${URL} — trying next"
        rm -f "${TMP_BAM}"
    fi
done

if [[ -z "$BAM_URL" ]]; then
    echo "[ERROR] Could not stream HLA reads for ${SAMPLE} from any URL" >&2
    exit 1
fi

# Coordinate-sort (required for BAM index) and index
samtools sort -@ 3 -m 2G -o "${BAM_OUT}" "${TMP_BAM}"
samtools index "${BAM_OUT}"
rm -f "${TMP_BAM}"

NREADS_FINAL=$(samtools view -c "${BAM_OUT}")
echo "[OK] HLA BAM: ${BAM_OUT} (${NREADS_FINAL} reads, indexed)"
EXTRACTEOF

    # Substitute placeholders (heredoc can't expand variables inside 'EXTRACTEOF')
    sed -i \
        -e "s|LOGS_PLACEHOLDER|${LOGS_DIR}|g" \
        -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
        -e "s|BAM_DIR_PLACEHOLDER|${BAM_DIR}|g" \
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
    echo "[SKIP] Phase 0 — using existing HLA BAMs"
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
    SAMPLESHEET="${SCRATCH_BASE}/wes_bam_samplesheet.csv"
    cat > "${SCRATCH_BASE}/gen_samplesheet.sh" << GENEOF
#!/bin/bash
echo "sample_id,bam_path" > "${SAMPLESHEET}"
while IFS= read -r SAMPLE; do
    BAM="${BAM_DIR}/\${SAMPLE}_hla.bam"
    if [[ -f "\$BAM" ]] && [[ -f "\${BAM}.bai" ]]; then
        echo "\${SAMPLE},\${BAM}" >> "${SAMPLESHEET}"
    else
        echo "[WARN] Missing HLA BAM for \${SAMPLE} — skipped" >&2
    fi
done < "${ALL_SAMPLES_LIST}"
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
module load apptainer 2>/dev/null || true

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
    --reference hg19 \\
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
