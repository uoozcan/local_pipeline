#!/bin/bash
# =============================================================================
# run_rna_calibration_puhti.sh
# RNA-seq HLA calibration on CSC Puhti — Geuvadis ERP001942 paired FASTQs
#
# Uses Geuvadis ERP001942 RNA-seq FASTQs (LCL cell lines, 1KGP Phase 1 samples).
# Ground truth: Gourraud et al. 2014 (same as WGS/WES calibration).
# Produces: conf/tool_weights_rna_v1.json, conf/tool_accuracy_rna_v1.tsv
#
# Usage:
#   bash scripts/run_rna_calibration_puhti.sh [OPTIONS]
#
# Options:
#   --project PROJECT_ID   CSC project account (default: $SLURM_JOB_ACCOUNT or project_2008084)
#   --sample-list FILE     Sample list (default: conf/rna_samples_50.txt)
#   --batch-size N         Process next N undownloaded samples per run (default: 10)
#   --tools TOOLS          Comma-separated tools (default: spechla,hlahd,arcashla,optitype,t1k,seq2hla)
#   --genes GENES          Comma-separated genes (default: A,B,C,DRB1,DQB1)
#   --skip-download        Skip Phase 0 (FASTQs already present)
#   --skip-typing          Skip Phase 1 (results already present)
#   --calibrate-only       Skip to Phase 2 (calibrate from existing results)
#   --keep-inputs          Do NOT delete FASTQs after Phase 1b (default: delete to save space)
#   --status               Show current progress and exit
#   --dry-run              Print commands without submitting
#   -h, --help             Show this help
#
# Batching: each run downloads and types the next --batch-size samples without FASTQs.
# Re-run after each batch completes. Calibration accumulates all results.
# Example (5 runs of 10):
#   bash scripts/run_rna_calibration_puhti.sh --project project_2008084 --batch-size 10
#
# Phase 0: SLURM array — download paired FASTQs from ENA FTP (ERP001942)
# Phase 1: SLURM job  — run Nextflow typing batch (seq_type=rna)
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
TOOLS="spechla,hlahd,arcashla,optitype,t1k,seq2hla"
GENES="A,B,C,DRB1,DQB1"
RESOLUTION="2-field"
SAMPLE_LIST_DEFAULT="${INSTALL_DIR}/conf/rna_samples_50.txt"
SAMPLE_LIST=""
BATCH_SIZE=10   # samples per run; 0 = all
SKIP_DOWNLOAD=false
SKIP_TYPING=false
CALIBRATE_ONLY=false
STATUS_ONLY=false
DRY_RUN=false
KEEP_INPUTS=false   # set true to retain downloaded FASTQs after Phase 1b

# ENA project for Geuvadis RNA-seq
ENA_PROJECT="ERP001942"
ENA_META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=sample_alias,sample_accession,run_accession,fastq_ftp&format=tsv"

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)     PROJECT_ID="$2";       shift 2 ;;
        --sample-list) SAMPLE_LIST="$2";      shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";       shift 2 ;;
        --tools)       TOOLS="$2";            shift 2 ;;
        --genes)       GENES="$2";            shift 2 ;;
        --skip-download)  SKIP_DOWNLOAD=true; shift ;;
        --skip-typing)    SKIP_TYPING=true;   shift ;;
        --calibrate-only) CALIBRATE_ONLY=true; SKIP_DOWNLOAD=true; SKIP_TYPING=true; shift ;;
        --keep-inputs)    KEEP_INPUTS=true;   shift ;;
        --status)      STATUS_ONLY=true;      shift ;;
        --dry-run)     DRY_RUN=true;          shift ;;
        -h|--help)
            sed -n '2,31p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$SAMPLE_LIST" ]] && SAMPLE_LIST="$SAMPLE_LIST_DEFAULT"

#-----------------------------------------------------------------------------
# Derived paths
#-----------------------------------------------------------------------------
SCRATCH_BASE="/scratch/${PROJECT_ID}/hla_calibration/rna"
FASTQ_DIR="${SCRATCH_BASE}/fastqs"
RESULTS_DIR="${SCRATCH_BASE}/results"
BY_TOOL_DIR="${RESULTS_DIR}/by_tool"
INDEX_DIR="${SCRATCH_BASE}/index"
LOGS_DIR="${SCRATCH_BASE}/logs"
GT_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_hla_gt.tsv"
POP_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_populations.tsv"
WEIGHTS_OUT="${INSTALL_DIR}/conf/tool_weights_rna_v1.json"
TABLE_OUT="${INSTALL_DIR}/conf/tool_accuracy_rna_v1.tsv"
ENA_META="${INDEX_DIR}/geuvadis_run_table.tsv"
URL_MAP="${INDEX_DIR}/sample_fastq_urls.tsv"
PROCESSED_DIR="${SCRATCH_BASE}/processed"   # sentinel files: ${SAMPLE}.done after cleanup

#-----------------------------------------------------------------------------
# Build effective sample list for this batch
# - Strip comments/blanks from master list
# - Filter to samples whose FASTQs do NOT yet exist (--batch-size 0 = all pending)
# - Take first BATCH_SIZE of those
#-----------------------------------------------------------------------------
EFFECTIVE_LIST="${SCRATCH_BASE}/rna_samples_effective.txt"
ALL_SAMPLES_LIST="${SCRATCH_BASE}/rna_samples_all.txt"

mkdir -p "${SCRATCH_BASE}" "${FASTQ_DIR}" "${RESULTS_DIR}" "${INDEX_DIR}" "${LOGS_DIR}" "${PROCESSED_DIR}"

# All samples (comments/blanks stripped)
grep -v '^#' "${SAMPLE_LIST}" | grep -v '^[[:space:]]*$' > "${ALL_SAMPLES_LIST}"
N_ALL=$(wc -l < "${ALL_SAMPLES_LIST}")

#-----------------------------------------------------------------------------
# GT filter: only keep samples present in the Gourraud 2014 ground truth
# Applied to master list so N_ALL/N_PENDING/N_DONE all reflect GT-only counts
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
    FILTERED_LIST="${SCRATCH_BASE}/rna_samples_gt_confirmed.txt"
    python3 - << GTEOF
import sys
gt_samples = set()
with open("${GT_FILE}") as f:
    for line in f:
        s = line.split('\t')[0].strip() if line.strip() else ''
        if s and s != 'Sample':
            gt_samples.add(s)
kept, skipped = [], []
with open("${ALL_SAMPLES_LIST}") as f:
    for line in f:
        s = line.strip()
        if s:
            (kept if s in gt_samples else skipped).append(s)
with open("${FILTERED_LIST}", "w") as f:
    f.write("\n".join(kept) + ("\n" if kept else ""))
if skipped:
    print(f"[WARN] {len(skipped)} sample(s) not in GT, excluded: {', '.join(skipped)}", file=sys.stderr)
print(f"[INFO] GT-confirmed: {len(kept)} / {len(kept)+len(skipped)} samples in master list")
GTEOF
    ALL_SAMPLES_LIST="${FILTERED_LIST}"
    N_ALL=$(wc -l < "${ALL_SAMPLES_LIST}")
fi

# Identify pending samples: (R1/R2 missing/empty) AND no processed sentinel
PENDING_LIST="${SCRATCH_BASE}/rna_samples_pending.txt"
> "${PENDING_LIST}"
while IFS= read -r S; do
    R1="${FASTQ_DIR}/${S}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${S}_R2.fastq.gz"
    FASTQS_MISSING=false
    { [[ ! -f "${R1}" ]] || [[ ! -s "${R1}" ]] || [[ ! -f "${R2}" ]] || [[ ! -s "${R2}" ]]; } \
        && FASTQS_MISSING=true
    if [[ "${FASTQS_MISSING}" == "true" ]] && [[ ! -f "${PROCESSED_DIR}/${S}.done" ]]; then
        echo "$S" >> "${PENDING_LIST}"
    fi
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
    echo "=== RNA-seq Calibration Status ==="
    echo "Project:     ${PROJECT_ID}"
    echo "Master list: ${N_ALL} samples total"
    echo ""
    echo "Phase 0 — FASTQs downloaded / processed:"
    N_PROCESSED=$(find "${PROCESSED_DIR}" -maxdepth 1 -name "*.done" 2>/dev/null | wc -l || true)
    N_PROCESSED=$(echo "${N_PROCESSED}" | tr -d '[:space:]'); N_PROCESSED=${N_PROCESSED:-0}
    N_FASTQ_PRESENT=$(find "${FASTQ_DIR}" -maxdepth 1 -name "*_R1.fastq.gz" 2>/dev/null | wc -l || true)
    N_FASTQ_PRESENT=$(echo "${N_FASTQ_PRESENT}" | tr -d '[:space:]'); N_FASTQ_PRESENT=${N_FASTQ_PRESENT:-0}
    echo "  ${N_DONE} / ${N_ALL} samples ready  (${N_PENDING} pending)"
    echo "  ${N_FASTQ_PRESENT} FASTQ pairs on disk   ${N_PROCESSED} cleaned (sentinel)"
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

# Nothing to download?
if [[ "$SKIP_DOWNLOAD" == "false" && "$CALIBRATE_ONLY" == "false" && "${N_TOTAL}" -eq 0 ]]; then
    echo "[INFO] All ${N_ALL} samples already processed (FASTQs present or cleaned). Use --calibrate-only to run calibration."
    exit 0
fi

echo "==================================================================="
echo " RNA-seq HLA Calibration — CSC Puhti"
echo "==================================================================="
echo " Project:      ${PROJECT_ID}"
echo " Samples:      ${N_TOTAL} this batch  (${N_DONE}/${N_ALL} total done; ${N_PENDING} pending)"
echo " Tools:        ${TOOLS}"
echo " Genes:        ${GENES}"
echo " FASTQs:       ${FASTQ_DIR}"
echo " Results:      ${RESULTS_DIR}"
echo " Weights out:  ${WEIGHTS_OUT}"
echo " Keep inputs:  ${KEEP_INPUTS}  (false = delete FASTQs after Phase 1b)"
echo " Dry-run:      ${DRY_RUN}"
echo "==================================================================="
echo ""

#-----------------------------------------------------------------------------
# Download Geuvadis ENA run table and build FASTQ URL map
#
# ENA metadata maps sample_accession (NA/HG ID) → ERR accession → fastq_ftp
# fastq_ftp field: semicolon-separated FTP paths, e.g.:
#   ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188022/ERR188022_1.fastq.gz;...
# Samples with multiple runs are concatenated by wget (each downloaded separately
# then combined, or just take the first run per sample for simplicity).
#-----------------------------------------------------------------------------
if [[ ! -f "${ENA_META}" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Downloading Geuvadis ENA run table..."
    curl -s "${ENA_META_URL}" > "${ENA_META}" \
        || wget -q -O "${ENA_META}" "${ENA_META_URL}" \
        || { echo "[ERROR] Could not download ENA run table"; exit 1; }
    echo "[OK] ENA metadata: ${ENA_META} ($(wc -l < "${ENA_META}") entries)"
else
    [[ "$DRY_RUN" == "false" ]] && echo "[INFO] Using cached ENA metadata: ${ENA_META}"
fi

if [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Building sample→FASTQ URL map from ENA metadata..."
    python3 - << PYEOF
import sys, os

meta_file   = "${ENA_META}"
sample_list = "${ALL_SAMPLES_LIST}"
url_map_out = "${URL_MAP}"

# Load sample list
with open(sample_list) as fh:
    samples = set(l.strip() for l in fh if l.strip() and not l.startswith("#"))

# Parse ENA metadata TSV: sample_accession, run_accession, fastq_ftp
# fastq_ftp: semicolon-separated paths (e.g. R1;R2 or R1;R2;unpaired)
# Use first run per sample (some samples have multiple runs; first is largest)
sample_runs = {}   # sample_id -> [(r1_url, r2_url), ...]
with open(meta_file) as fh:
    header = None
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            # Find column indices
            try:
                i_sample = header.index("sample_alias")   # 1KGP ID (NA18501 etc.), not SAMEA
                i_run    = header.index("run_accession")
                i_ftp    = header.index("fastq_ftp")
            except ValueError:
                print(f"[ERROR] Unexpected ENA header: {header}", file=sys.stderr)
                sys.exit(1)
            continue
        if len(parts) <= max(i_sample, i_run, i_ftp):
            continue
        sample_id = parts[i_sample].strip()
        # Strip project prefix if present (e.g. "GEUV:NA18501" -> "NA18501")
        if ':' in sample_id:
            sample_id = sample_id.split(':', 1)[1]
        if sample_id not in samples:
            continue
        ftp_paths = [p.strip() for p in parts[i_ftp].split(";") if p.strip()]
        # Expect paired FASTQs: _1.fastq.gz and _2.fastq.gz
        r1 = next((p for p in ftp_paths if "_1.fastq.gz" in p), None)
        r2 = next((p for p in ftp_paths if "_2.fastq.gz" in p), None)
        if r1 and r2:
            if sample_id not in sample_runs:
                sample_runs[sample_id] = []
            sample_runs[sample_id].append((r1, r2))

found, missing, multi = 0, [], []
with open(url_map_out, "w") as fout:
    for sample in sorted(samples):
        runs = sample_runs.get(sample)
        if not runs:
            missing.append(sample)
            continue
        if len(runs) > 1:
            multi.append(f"{sample}({len(runs)} runs)")
        # Write all runs as pipe-separated alternates after first run
        r1_primary, r2_primary = runs[0]
        r1_alts = "|".join(r for r, _ in runs[1:]) if len(runs) > 1 else ""
        r2_alts = "|".join(r for _, r in runs[1:]) if len(runs) > 1 else ""
        fout.write(f"{sample}\t{r1_primary}\t{r2_primary}\t{r1_alts}\t{r2_alts}\n")
        found += 1

print(f"[OK] URL map: {found} samples mapped, {len(missing)} missing from ENA metadata")
if missing:
    print(f"[WARN] Not in ENA metadata: {', '.join(missing)}", file=sys.stderr)
if multi:
    print(f"[INFO] Multiple runs (first used): {', '.join(multi)}")
PYEOF
fi

#-----------------------------------------------------------------------------
# Phase 0: SLURM array — download paired FASTQs from ENA FTP
#-----------------------------------------------------------------------------
PHASE0_DEP=""

if [[ "$SKIP_DOWNLOAD" == "false" ]]; then
    echo "=== Phase 0: Submitting FASTQ download array (${N_TOTAL} tasks) ==="

    DOWNLOAD_SCRIPT="${SCRATCH_BASE}/phase0_download.sh"
    cat > "$DOWNLOAD_SCRIPT" << 'DLEOF'
#!/bin/bash
#SBATCH --job-name=rna_download
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=LOGS_PLACEHOLDER/download_%A_%a.out
#SBATCH --error=LOGS_PLACEHOLDER/download_%A_%a.err

set -euo pipefail

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "EFFECTIVE_LIST_PLACEHOLDER")
[[ -z "$SAMPLE" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

FASTQ_DIR="FASTQ_DIR_PLACEHOLDER"
URL_MAP="URL_MAP_PLACEHOLDER"

R1_OUT="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2_OUT="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Date:   $(date)"

# Skip if both FASTQs already present and non-empty
if [[ -f "${R1_OUT}" ]] && [[ -s "${R1_OUT}" ]] && [[ -f "${R2_OUT}" ]] && [[ -s "${R2_OUT}" ]]; then
    echo "[SKIP] FASTQs already present for ${SAMPLE}"
    exit 0
fi

# Look up FASTQ URLs (col 2 = R1 primary, col 3 = R2 primary, col 4 = R1 alts, col 5 = R2 alts)
LINE=$(grep -P "^${SAMPLE}\t" "$URL_MAP" || true)
if [[ -z "$LINE" ]]; then
    echo "[ERROR] No FASTQ URL found for sample ${SAMPLE} in URL map" >&2
    exit 1
fi
R1_URL=$(echo "$LINE" | cut -f2)
R2_URL=$(echo "$LINE" | cut -f3)
R1_ALTS=$(echo "$LINE" | cut -f4 | tr '|' ' ')
R2_ALTS=$(echo "$LINE" | cut -f5 | tr '|' ' ')

mkdir -p "${FASTQ_DIR}"

# Download R1 — try primary URL, then alts
download_file() {
    local URL="$1"
    local OUT="$2"
    # ENA FTP paths don't include ftp:// prefix; prepend it
    local FULL_URL="ftp://${URL}"
    echo "[INFO] Downloading: ${FULL_URL}"
    wget -q --tries=3 --timeout=120 -O "${OUT}.tmp" "${FULL_URL}" && mv "${OUT}.tmp" "${OUT}"
}

R1_OK=false
for URL in ${R1_URL} ${R1_ALTS}; do
    if download_file "${URL}" "${R1_OUT}"; then
        echo "[OK] R1: ${R1_OUT} ($(du -sh "${R1_OUT}" | cut -f1))"
        R1_OK=true
        break
    else
        echo "[WARN] Failed: ${URL} — trying next"
        rm -f "${R1_OUT}.tmp"
    fi
done

R2_OK=false
for URL in ${R2_URL} ${R2_ALTS}; do
    if download_file "${URL}" "${R2_OUT}"; then
        echo "[OK] R2: ${R2_OUT} ($(du -sh "${R2_OUT}" | cut -f1))"
        R2_OK=true
        break
    else
        echo "[WARN] Failed: ${URL} — trying next"
        rm -f "${R2_OUT}.tmp"
    fi
done

if [[ "$R1_OK" == "false" ]] || [[ "$R2_OK" == "false" ]]; then
    echo "[ERROR] Could not download FASTQs for ${SAMPLE}" >&2
    exit 1
fi

echo "[OK] ${SAMPLE}: paired FASTQs downloaded"
DLEOF

    # Substitute placeholders
    sed -i \
        -e "s|LOGS_PLACEHOLDER|${LOGS_DIR}|g" \
        -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
        -e "s|FASTQ_DIR_PLACEHOLDER|${FASTQ_DIR}|g" \
        -e "s|URL_MAP_PLACEHOLDER|${URL_MAP}|g" \
        "$DOWNLOAD_SCRIPT"

    CONCURRENT=$([[ "${BATCH_SIZE}" -gt 0 ]] && echo "${BATCH_SIZE}" || echo "10")
    SBATCH_DOWNLOAD="sbatch --parsable \
        --account=${PROJECT_ID} \
        --array=1-${N_TOTAL}%${CONCURRENT} \
        ${DOWNLOAD_SCRIPT}"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_DOWNLOAD}"
        DOWNLOAD_JOB="DRY_RUN_JOB"
    else
        DOWNLOAD_JOB=$(eval "$SBATCH_DOWNLOAD")
        echo "[INFO] Phase 0 download array submitted: ${DOWNLOAD_JOB}"
        PHASE0_DEP="--dependency=afterok:${DOWNLOAD_JOB}"
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

    SAMPLESHEET="${SCRATCH_BASE}/rna_fastq_samplesheet.csv"

    # gen_samplesheet.sh: iterates ALL_SAMPLES_LIST so all downloaded FASTQs are included
    cat > "${SCRATCH_BASE}/gen_samplesheet.sh" << GENEOF
#!/bin/bash
echo "sample_id,fastq_1,fastq_2" > "${SAMPLESHEET}"
while IFS= read -r SAMPLE; do
    R1="${FASTQ_DIR}/\${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/\${SAMPLE}_R2.fastq.gz"
    if [[ -f "\$R1" ]] && [[ -s "\$R1" ]] && [[ -f "\$R2" ]] && [[ -s "\$R2" ]]; then
        echo "\${SAMPLE},\${R1},\${R2}" >> "${SAMPLESHEET}"
    else
        echo "[WARN] Missing FASTQs for \${SAMPLE} — skipped" >&2
    fi
done < "${ALL_SAMPLES_LIST}"
echo "[OK] Samplesheet: ${SAMPLESHEET} (\$(tail -n +2 ${SAMPLESHEET} | wc -l) samples)"
GENEOF

    echo "[INFO] Samplesheet will be generated in Phase 1 job: ${SAMPLESHEET}"

    TYPING_SCRIPT="${SCRATCH_BASE}/phase1_typing.sh"
    cat > "$TYPING_SCRIPT" << TYPINGEOF
#!/bin/bash
#SBATCH --job-name=rna_typing
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

echo "=== Phase 1: HLA Typing (RNA-seq) ==="
echo "Date: \$(date)"
echo "Samplesheet: ${SAMPLESHEET}"

cd "${INSTALL_DIR}"

# Regenerate samplesheet (captures all FASTQs downloaded in Phase 0)
bash "${SCRATCH_BASE}/gen_samplesheet.sh"

N_SAMPLES_SHEET=\$(tail -n +2 "${SAMPLESHEET}" | wc -l)
echo "Samples to type: \${N_SAMPLES_SHEET}"

nextflow run main.nf \\
    --input_samplesheet "${SAMPLESHEET}" \\
    --seq_type rna \\
    --tools "${TOOLS}" \\
    --weighting equal \\
    --skip_qc true \\
    --project "${PROJECT_ID}" \\
    --outdir "${RESULTS_DIR}" \\
    -c conf/puhti.config \\
    -profile apptainer \\
    -resume \\
    -with-trace "${RESULTS_DIR}/pipeline_info/trace_rna.txt" \\
    -with-report "${RESULTS_DIR}/pipeline_info/report_rna.html" \\
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
#SBATCH --job-name=rna_collect
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/collect_%j.out
#SBATCH --error=${LOGS_DIR}/collect_%j.err

set -euo pipefail
echo "=== Phase 1b: Collecting RNA-seq results ==="
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
done < "${ALL_SAMPLES_LIST}"

echo ""
echo "Symlinks created: \$TOTAL"
echo ""
echo "Results by tool:"
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    N=\$(find "${BY_TOOL_DIR}/\${TOOL}/" -maxdepth 1 -name "*_\${TOOL}.txt" 2>/dev/null | wc -l)
    printf "  %-12s %3d samples\n" "\$TOOL" "\$N"
done

# Cleanup: delete FASTQs to free scratch space (skip if --keep-inputs)
# Geuvadis FASTQs are 3-8 GB/pair; 50 samples could use 150-400 GB without cleanup
if [[ "KEEP_INPUTS_PLACEHOLDER" == "false" ]]; then
    echo ""
    echo "=== Cleanup: deleting FASTQs to free scratch space ==="
    CLEANED=0
    SKIPPED_CLEAN=0
    mkdir -p "PROCESSED_DIR_PLACEHOLDER"
    while IFS= read -r SAMPLE; do
        N_RES=\$(find "BY_TOOL_DIR_PLACEHOLDER/" -maxdepth 2 \
                      -name "\${SAMPLE}_*.txt" 2>/dev/null | wc -l || true)
        N_RES=\${N_RES:-0}
        if [[ "\${N_RES}" -gt 0 ]]; then
            rm -f "FASTQ_DIR_PLACEHOLDER/\${SAMPLE}_R1.fastq.gz" \
                  "FASTQ_DIR_PLACEHOLDER/\${SAMPLE}_R2.fastq.gz"
            touch "PROCESSED_DIR_PLACEHOLDER/\${SAMPLE}.done"
            echo "[OK] Cleaned: \${SAMPLE} (\${N_RES} tool results)"
            CLEANED=\$((CLEANED+1))
        else
            echo "[WARN] No results for \${SAMPLE} — keeping FASTQs" >&2
            SKIPPED_CLEAN=\$((SKIPPED_CLEAN+1))
        fi
    done < "EFFECTIVE_LIST_PLACEHOLDER"
    echo "Cleaned: \${CLEANED} FASTQ pairs deleted, \${SKIPPED_CLEAN} kept (no results yet)"
else
    echo ""
    echo "[INFO] --keep-inputs set: FASTQs retained"
fi
COLLECTEOF

sed -i \
    -e "s|KEEP_INPUTS_PLACEHOLDER|${KEEP_INPUTS}|g" \
    -e "s|PROCESSED_DIR_PLACEHOLDER|${PROCESSED_DIR}|g" \
    -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
    -e "s|BY_TOOL_DIR_PLACEHOLDER|${BY_TOOL_DIR}|g" \
    -e "s|FASTQ_DIR_PLACEHOLDER|${FASTQ_DIR}|g" \
    "${COLLECT_SCRIPT}"

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
COMPARE_OUT="${SCRATCH_BASE}/strategy_comparison_rna_calibrated.tsv"
COMPARE_RC_OUT="${SCRATCH_BASE}/strategy_comparison_rna_rc.tsv"

cat > "$CALIBRATE_SCRIPT" << CALEOF
#!/bin/bash
#SBATCH --job-name=rna_calibrate
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/calibrate_%j.out
#SBATCH --error=${LOGS_DIR}/calibrate_%j.err

set -euo pipefail
echo "=== Phase 2: RNA-seq HLA Tool Weight Calibration ==="
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
echo "[INFO] Running RNA-seq calibration (genes: ${GENES})..."
POP_FILE_ARG=""
[[ -f "${POP_FILE}" ]] && POP_FILE_ARG="--population-file ${POP_FILE}"

# shellcheck disable=SC2086
python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" calibrate \\
    --ground-truth "${GT_FILE}" \\
    --results-dir "${BY_TOOL_DIR}" \\
    --data-type rna \\
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
    print(f'\nRNA-seq calibrated weights (n={n} samples):')
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
echo " RNA-seq Calibration Complete"
echo "================================================================="
echo " Weights:      ${WEIGHTS_OUT}"
echo " Accuracy:     ${TABLE_OUT}"
[[ -f "${COMPARE_OUT}" ]]    && echo " Strat(cal):   ${COMPARE_OUT}"
[[ -f "${COMPARE_RC_OUT}" ]] && echo " Strat(rc):    ${COMPARE_RC_OUT}"
echo ""
echo " Next: push updated weights to GitHub and rerun pipeline"
echo "   --weighting calibrated --weights_file conf/tool_weights_rna_v1.json"
echo "================================================================="
CALEOF

SBATCH_CAL="sbatch --parsable \
    --dependency=afterok:${COLLECT_JOB} \
    ${CALIBRATE_SCRIPT}"

REMAINING=$(( N_PENDING - N_TOTAL ))
if [[ "$CALIBRATE_ONLY" == "true" ]] || [[ "$REMAINING" -eq 0 ]]; then
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_CAL}"
    else
        CALIBRATE_JOB=$(eval "$SBATCH_CAL")
        echo "[INFO] Phase 2 calibration job submitted: ${CALIBRATE_JOB}"
    fi
else
    echo "[INFO] ${REMAINING} sample(s) still pending — skipping Phase 2 calibration."
    echo "[INFO] Re-run script for next batch. Use --calibrate-only when all batches complete."
fi

echo ""
echo "==================================================================="
echo " All jobs submitted."
echo ""
echo " Monitor: squeue -u \$USER"
echo " Status:  bash scripts/run_rna_calibration_puhti.sh --status"
echo " Logs:    ${LOGS_DIR}/"
echo ""
echo " Expected outputs:"
echo "   ${WEIGHTS_OUT}"
echo "   ${TABLE_OUT}"
echo "==================================================================="
