#!/bin/bash
# =============================================================================
# run_wgs_calibration_puhti.sh
# WGS HLA calibration on CSC Puhti — 1KGP NYGC 30x CRAMs (GRCh38)
#
# Uses 1KGP NYGC 30x high-coverage CRAMs (PRJEB31736) from EBI SRA.
# Ground truth: Gourraud et al. 2014 (same as WES calibration).
# Produces: conf/tool_weights_wgs_puhti_v1.json, conf/tool_accuracy_wgs_puhti_v1.tsv
#
# Usage:
#   bash scripts/run_wgs_calibration_puhti.sh [OPTIONS]
#
# Options:
#   --project PROJECT_ID   CSC project account (default: $SLURM_JOB_ACCOUNT or project_2008084)
#   --sample-list FILE     Sample list (default: conf/wgs_samples_50.txt)
#   --batch-size N         Process next N unextracted samples per run (default: 10)
#   --tools TOOLS          Comma-separated tools (default: hlahd,spechla,arcashla,optitype,xhla,kourami,polysolver)
#   --genes GENES          Comma-separated genes (default: A,B,C,DRB1,DQB1)
#   --skip-extract         Skip Phase 0 (BAMs already present)
#   --skip-typing          Skip Phase 1 (results already present)
#   --calibrate-only       Skip to Phase 2 (calibrate from existing results)
#   --keep-inputs          Do NOT delete HLA BAMs after Phase 1b (default: delete to save space)
#   --status               Show current progress and exit
#   --dry-run              Print commands without submitting
#   -h, --help             Show this help
#
# Batching: each run extracts and types the next --batch-size samples without BAMs.
# Re-run after each batch completes. Calibration accumulates all results.
# Example (5 runs of 10):
#   bash scripts/run_wgs_calibration_puhti.sh --project project_2008084 --batch-size 10
#
# Phase 0: SLURM array — stream HLA region from EBI 30x CRAM → sorted BAM
# Phase 1: SLURM job  — run Nextflow typing batch (seq_type=dna, all BAM tools)
# Phase 1b: SLURM job — collect results into by_tool/ layout + delete BAMs
# Phase 2: SLURM job  — calibrate weights + compare voting strategies
# =============================================================================
set -euo pipefail

#-----------------------------------------------------------------------------
# Defaults
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="$(dirname "$SCRIPT_DIR")"   # hla_typing_pipeline/

PROJECT_ID="${SLURM_JOB_ACCOUNT:-project_2008084}"
# hlala excluded: needs full-genome PRG graph (not HLA-region BAM)
# seq2hla excluded: RNA-seq optimised, unreliable on WGS
TOOLS="hlahd,spechla,arcashla,optitype,kourami,polysolver"
GENES="A,B,C,DRB1,DQB1"
RESOLUTION="2-field"
SAMPLE_LIST_DEFAULT="${INSTALL_DIR}/conf/wgs_samples_50.txt"
SAMPLE_LIST=""
BATCH_SIZE=10   # samples per run; 0 = all
SKIP_EXTRACT=false
SKIP_TYPING=false
CALIBRATE_ONLY=false
STATUS_ONLY=false
DRY_RUN=false
KEEP_INPUTS=false   # set true to retain HLA BAMs after Phase 1b

# ENA project for 1KGP NYGC 30x WGS CRAMs
ENA_PROJECT="PRJEB31736"
ENA_META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=run_accession,submitted_ftp&format=tsv"

# HLA region (GRCh38/hg38 with UCSC chr prefix)
HLA_REGION="chr6:28000000-34000000"

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)    PROJECT_ID="$2";      shift 2 ;;
        --sample-list) SAMPLE_LIST="$2";   shift 2 ;;
        --batch-size) BATCH_SIZE="$2";     shift 2 ;;
        --tools)      TOOLS="$2";          shift 2 ;;
        --genes)      GENES="$2";          shift 2 ;;
        --skip-extract)   SKIP_EXTRACT=true;   shift ;;
        --skip-typing)    SKIP_TYPING=true;    shift ;;
        --calibrate-only) CALIBRATE_ONLY=true; SKIP_EXTRACT=true; SKIP_TYPING=true; shift ;;
        --keep-inputs)    KEEP_INPUTS=true;    shift ;;
        --status)     STATUS_ONLY=true;    shift ;;
        --dry-run)    DRY_RUN=true;        shift ;;
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
SCRATCH_BASE="/scratch/${PROJECT_ID}/hla_calibration/wgs"
BAM_DIR="${SCRATCH_BASE}/bams"
RESULTS_DIR="${SCRATCH_BASE}/results"
BY_TOOL_DIR="${RESULTS_DIR}/by_tool"
INDEX_DIR="${SCRATCH_BASE}/index"
LOGS_DIR="${SCRATCH_BASE}/logs"
CRAM_REF_CACHE="${SCRATCH_BASE}/cram_ref_cache"
PROCESSED_DIR="${SCRATCH_BASE}/processed"   # sentinel files: ${SAMPLE}.done after cleanup
GT_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_hla_gt.tsv"
POP_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_populations.tsv"
WEIGHTS_OUT="${INSTALL_DIR}/conf/tool_weights_wgs_puhti_v1.json"
TABLE_OUT="${INSTALL_DIR}/conf/tool_accuracy_wgs_puhti_v1.tsv"
ENA_META="${INDEX_DIR}/nygc_30x_run_table.tsv"
URL_MAP="${INDEX_DIR}/sample_cram_urls.tsv"

#-----------------------------------------------------------------------------
# Build effective sample list for this batch
#-----------------------------------------------------------------------------
EFFECTIVE_LIST="${SCRATCH_BASE}/wgs_samples_effective.txt"
ALL_SAMPLES_LIST="${SCRATCH_BASE}/wgs_samples_all.txt"

mkdir -p "${SCRATCH_BASE}" "${BAM_DIR}" "${RESULTS_DIR}" "${INDEX_DIR}" \
         "${LOGS_DIR}" "${CRAM_REF_CACHE}" "${PROCESSED_DIR}"

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
    FILTERED_LIST="${SCRATCH_BASE}/wgs_samples_gt_confirmed.txt"
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

# Identify pending samples: no HLA BAM AND no processed sentinel
PENDING_LIST="${SCRATCH_BASE}/wgs_samples_pending.txt"
> "${PENDING_LIST}"
while IFS= read -r S; do
    if [[ ! -f "${BAM_DIR}/${S}_hla.bam" ]] && [[ ! -f "${PROCESSED_DIR}/${S}.done" ]]; then
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
    echo "=== WGS Calibration Status ==="
    echo "Project:     ${PROJECT_ID}"
    echo "Master list: ${N_ALL} samples total"
    echo ""
    echo "Phase 0 — HLA BAMs / processing:"
    N_DONE=$(( N_ALL - N_PENDING ))
    N_PROCESSED=$(find "${PROCESSED_DIR}" -maxdepth 1 -name "*.done" 2>/dev/null | wc -l || true)
    N_PROCESSED=$(echo "${N_PROCESSED}" | tr -d '[:space:]'); N_PROCESSED=${N_PROCESSED:-0}
    N_BAM_PRESENT=$(find "${BAM_DIR}" -maxdepth 1 -name "*_hla.bam" 2>/dev/null | wc -l || true)
    N_BAM_PRESENT=$(echo "${N_BAM_PRESENT}" | tr -d '[:space:]'); N_BAM_PRESENT=${N_BAM_PRESENT:-0}
    echo "  ${N_DONE} / ${N_ALL} samples ready  (${N_PENDING} pending)"
    echo "  ${N_BAM_PRESENT} HLA BAMs on disk   ${N_PROCESSED} cleaned (sentinel)"
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
    echo "[INFO] All ${N_ALL} samples already processed (BAMs present or cleaned). Use --calibrate-only to run calibration."
    exit 0
fi

echo "==================================================================="
echo " WGS HLA Calibration — CSC Puhti"
echo "==================================================================="
echo " Project:      ${PROJECT_ID}"
echo " Samples:      ${N_TOTAL} this batch  (${N_DONE}/${N_ALL} total done; ${N_PENDING} pending)"
echo " Tools:        ${TOOLS}"
echo " Genes:        ${GENES}"
echo " BAMs:         ${BAM_DIR}"
echo " Results:      ${RESULTS_DIR}"
echo " Weights out:  ${WEIGHTS_OUT}"
echo " Keep inputs:  ${KEEP_INPUTS}  (false = delete BAMs after Phase 1b)"
echo " Dry-run:      ${DRY_RUN}"
echo "==================================================================="
echo ""

#-----------------------------------------------------------------------------
# Download ENA run table for PRJEB31736 and build CRAM URL map
#
# ENA run table: run_accession, submitted_ftp
# submitted_ftp contains the CRAM FTP path, e.g.:
#   ftp.sra.ebi.ac.uk/.../NA19238.final.cram;ftp.sra.ebi.ac.uk/.../NA19238.final.cram.crai
# CRAM URL: https://ftp.sra.ebi.ac.uk/vol1/run/{ERR[0:6]}/{ERR}/{SAMPLE}.final.cram
# (using HTTPS rather than ftp:// — Puhti samtools 1.21 supports libcurl/HTTPS but not FTP)
#-----------------------------------------------------------------------------
if [[ ! -f "${ENA_META}" ]] && [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Downloading NYGC 30x ENA run table (PRJEB31736)..."
    curl -s "${ENA_META_URL}" > "${ENA_META}" \
        || wget -q -O "${ENA_META}" "${ENA_META_URL}" \
        || { echo "[ERROR] Could not download ENA run table"; exit 1; }
    echo "[OK] ENA metadata: ${ENA_META} ($(wc -l < "${ENA_META}") entries)"
else
    [[ "$DRY_RUN" == "false" ]] && echo "[INFO] Using cached ENA metadata: ${ENA_META}"
fi

if [[ "$DRY_RUN" == "false" ]]; then
    echo "[INFO] Building sample→CRAM URL map from ENA metadata..."
    python3 - << PYEOF
import sys, os, re

meta_file   = "${ENA_META}"
sample_list = "${ALL_SAMPLES_LIST}"
url_map_out = "${URL_MAP}"

# Load sample list
with open(sample_list) as fh:
    samples = set(l.strip() for l in fh if l.strip() and not l.startswith("#"))

# Parse ENA metadata: run_accession, submitted_ftp
# Extract sample ID from CRAM filename (e.g. NA19238.final.cram)
found, missing = 0, []
sample_map = {}   # sample_id -> (err_accession, cram_url)

with open(meta_file) as fh:
    header = None
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            try:
                i_run = header.index("run_accession")
                i_ftp = header.index("submitted_ftp")
            except ValueError:
                print(f"[ERROR] Unexpected ENA header: {header}", file=sys.stderr)
                sys.exit(1)
            continue
        if len(parts) <= max(i_run, i_ftp):
            continue
        err = parts[i_run].strip()
        ftp = parts[i_ftp].strip()
        # Extract sample ID from CRAM filename: {SAMPLE}.final.cram
        m = re.search(r'/([A-Z0-9]+)\.final\.cram(?:;|$)', ftp)
        if not m:
            continue
        s = m.group(1)
        if s not in samples:
            continue
        # Build HTTPS URL: https://ftp.sra.ebi.ac.uk/vol1/run/{ERR[0:6]}/{ERR}/{SAMPLE}.final.cram
        cram_url = f"https://ftp.sra.ebi.ac.uk/vol1/run/{err[:6]}/{err}/{s}.final.cram"
        if s not in sample_map:   # keep first ERR per sample
            sample_map[s] = (err, cram_url)

with open(url_map_out, "w") as fout:
    for s in sorted(samples):
        if s in sample_map:
            err, url = sample_map[s]
            fout.write(f"{s}\t{err}\t{url}\n")
            found += 1
        else:
            missing.append(s)

print(f"[OK] URL map: {found} samples mapped, {len(missing)} not found in ENA metadata")
if missing:
    print(f"[WARN] Not in PRJEB31736: {', '.join(missing)}", file=sys.stderr)
PYEOF
fi

#-----------------------------------------------------------------------------
# Phase 0: SLURM array — extract HLA region from 30x CRAM → sorted BAM
#-----------------------------------------------------------------------------
PHASE0_DEP=""

if [[ "$SKIP_EXTRACT" == "false" ]]; then
    echo "=== Phase 0: Submitting CRAM extraction array (${N_TOTAL} tasks) ==="

    EXTRACT_SCRIPT="${SCRATCH_BASE}/phase0_extract.sh"
    cat > "$EXTRACT_SCRIPT" << 'EXTRACTEOF'
#!/bin/bash
#SBATCH --job-name=wgs_extract
#SBATCH --partition=small
#SBATCH --time=03:00:00
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
REF_CACHE="REF_CACHE_PLACEHOLDER"

BAM_OUT="${BAM_DIR}/${SAMPLE}_hla.bam"

echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Date:   $(date)"

# Skip if already done
if [[ -f "${BAM_OUT}" ]] && [[ -f "${BAM_OUT}.bai" ]]; then
    echo "[SKIP] HLA BAM already present for ${SAMPLE}"
    exit 0
fi

# Look up CRAM URL
LINE=$(grep -P "^${SAMPLE}\t" "$URL_MAP" || true)
if [[ -z "$LINE" ]]; then
    echo "[WARN] No CRAM URL found for sample ${SAMPLE} — skipping" >&2
    exit 0
fi
CRAM_URL=$(echo "$LINE" | cut -f3)
echo "[INFO] CRAM URL: ${CRAM_URL}"

# CRAM MD5 reference: use ENA reference server (requires internet on compute nodes)
# REF_CACHE stores decoded reference chunks to avoid repeated downloads
export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${REF_CACHE}/%2s/%2s/%s"
mkdir -p "${REF_CACHE}"

mkdir -p "${BAM_DIR}"
TMP_BAM="${BAM_DIR}/${SAMPLE}_hla_tmp.bam"

echo "[INFO] Streaming HLA region from CRAM: ${HLA_REGION}"
if samtools view -b -h -o "${TMP_BAM}" "${CRAM_URL}" "${HLA_REGION}" 2>/dev/null; then
    NREADS=$(samtools view -c "${TMP_BAM}" 2>/dev/null || echo 0)
    if [[ "${NREADS}" -gt 0 ]]; then
        echo "[OK] HLA reads: ${NREADS}"
    else
        echo "[ERROR] 0 reads extracted from ${CRAM_URL}" >&2
        rm -f "${TMP_BAM}"
        exit 1
    fi
else
    echo "[ERROR] samtools failed for ${CRAM_URL}" >&2
    rm -f "${TMP_BAM}"
    exit 1
fi

# Coordinate-sort + index
samtools sort -@ 3 -m 2G -o "${BAM_OUT}" "${TMP_BAM}"
samtools index "${BAM_OUT}"
rm -f "${TMP_BAM}"

NREADS_FINAL=$(samtools view -c "${BAM_OUT}")
echo "[OK] HLA BAM: ${BAM_OUT} (${NREADS_FINAL} reads, indexed)"
EXTRACTEOF

    # Substitute placeholders
    sed -i \
        -e "s|LOGS_PLACEHOLDER|${LOGS_DIR}|g" \
        -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
        -e "s|BAM_DIR_PLACEHOLDER|${BAM_DIR}|g" \
        -e "s|URL_MAP_PLACEHOLDER|${URL_MAP}|g" \
        -e "s|HLA_REGION_PLACEHOLDER|${HLA_REGION}|g" \
        -e "s|REF_CACHE_PLACEHOLDER|${CRAM_REF_CACHE}|g" \
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

    SAMPLESHEET="${SCRATCH_BASE}/wgs_bam_samplesheet.csv"

    # gen_samplesheet.sh iterates ALL_SAMPLES_LIST so all available BAMs are included
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

    echo "[INFO] Samplesheet will be generated in Phase 1 job: ${SAMPLESHEET}"

    TYPING_SCRIPT="${SCRATCH_BASE}/phase1_typing.sh"
    cat > "$TYPING_SCRIPT" << TYPINGEOF
#!/bin/bash
#SBATCH --job-name=wgs_typing
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

echo "=== Phase 1: HLA Typing (WGS 30x) ==="
echo "Date: \$(date)"
echo "Samplesheet: ${SAMPLESHEET}"

cd "${INSTALL_DIR}"

# Regenerate samplesheet (captures all BAMs from Phase 0)
bash "${SCRATCH_BASE}/gen_samplesheet.sh"

N_SAMPLES_SHEET=\$(tail -n +2 "${SAMPLESHEET}" | wc -l)
echo "Samples to type: \${N_SAMPLES_SHEET}"

nextflow run main.nf \\
    --input_samplesheet "${SAMPLESHEET}" \\
    --seq_type dna \\
    --reference hg38 \\
    --tools "${TOOLS}" \\
    --weighting equal \\
    --skip_qc true \\
    --project "${PROJECT_ID}" \\
    --outdir "${RESULTS_DIR}" \\
    -c conf/puhti.config \\
    -profile apptainer \\
    -resume \\
    -with-trace "${RESULTS_DIR}/pipeline_info/trace_wgs.txt" \\
    -with-report "${RESULTS_DIR}/pipeline_info/report_wgs.html" \\
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
# Phase 1b: Collect results + cleanup BAMs
#-----------------------------------------------------------------------------
COLLECT_SCRIPT="${SCRATCH_BASE}/phase1b_collect.sh"
cat > "$COLLECT_SCRIPT" << COLLECTEOF
#!/bin/bash
#SBATCH --job-name=wgs_collect
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/collect_%j.out
#SBATCH --error=${LOGS_DIR}/collect_%j.err

set -euo pipefail
echo "=== Phase 1b: Collecting WGS results ==="
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

# Cleanup: delete HLA BAMs to free scratch space (skip if --keep-inputs)
if [[ "KEEP_INPUTS_PLACEHOLDER" == "false" ]]; then
    echo ""
    echo "=== Cleanup: deleting HLA BAMs to free scratch space ==="
    CLEANED=0
    SKIPPED_CLEAN=0
    mkdir -p "PROCESSED_DIR_PLACEHOLDER"
    while IFS= read -r SAMPLE; do
        N_RES=\$(find "BY_TOOL_DIR_PLACEHOLDER/" -maxdepth 2 \
                      -name "\${SAMPLE}_*.txt" 2>/dev/null | wc -l || true)
        N_RES=\${N_RES:-0}
        if [[ "\${N_RES}" -gt 0 ]]; then
            rm -f "BAM_DIR_PLACEHOLDER/\${SAMPLE}_hla.bam" \
                  "BAM_DIR_PLACEHOLDER/\${SAMPLE}_hla.bam.bai"
            touch "PROCESSED_DIR_PLACEHOLDER/\${SAMPLE}.done"
            echo "[OK] Cleaned: \${SAMPLE} (\${N_RES} tool results)"
            CLEANED=\$((CLEANED+1))
        else
            echo "[WARN] No results for \${SAMPLE} — keeping BAM" >&2
            SKIPPED_CLEAN=\$((SKIPPED_CLEAN+1))
        fi
    done < "EFFECTIVE_LIST_PLACEHOLDER"
    echo "Cleaned: \${CLEANED} BAMs deleted, \${SKIPPED_CLEAN} kept (no results yet)"
else
    echo ""
    echo "[INFO] --keep-inputs set: HLA BAMs retained"
fi
COLLECTEOF

sed -i \
    -e "s|KEEP_INPUTS_PLACEHOLDER|${KEEP_INPUTS}|g" \
    -e "s|PROCESSED_DIR_PLACEHOLDER|${PROCESSED_DIR}|g" \
    -e "s|EFFECTIVE_LIST_PLACEHOLDER|${EFFECTIVE_LIST}|g" \
    -e "s|BY_TOOL_DIR_PLACEHOLDER|${BY_TOOL_DIR}|g" \
    -e "s|BAM_DIR_PLACEHOLDER|${BAM_DIR}|g" \
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
COMPARE_OUT="${SCRATCH_BASE}/strategy_comparison_wgs_calibrated.tsv"
COMPARE_RC_OUT="${SCRATCH_BASE}/strategy_comparison_wgs_rc.tsv"

cat > "$CALIBRATE_SCRIPT" << CALEOF
#!/bin/bash
#SBATCH --job-name=wgs_calibrate
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/calibrate_%j.out
#SBATCH --error=${LOGS_DIR}/calibrate_%j.err

set -euo pipefail
echo "=== Phase 2: WGS HLA Tool Weight Calibration ==="
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

echo ""
echo "[INFO] Running WGS calibration (genes: ${GENES})..."
POP_FILE_ARG=""
[[ -f "${POP_FILE}" ]] && POP_FILE_ARG="--population-file ${POP_FILE}"

# shellcheck disable=SC2086
python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" calibrate \\
    --ground-truth "${GT_FILE}" \\
    --results-dir "${BY_TOOL_DIR}" \\
    --data-type wgs \\
    --genes "${GENES}" \\
    --resolution "${RESOLUTION}" \\
    --output-weights "${WEIGHTS_OUT}" \\
    --output-table  "${TABLE_OUT}" \\
    \$POP_FILE_ARG

echo ""
echo "[OK] Calibration outputs:"
echo "     Weights: ${WEIGHTS_OUT}"
echo "     Table:   ${TABLE_OUT}"

python3 - << PYEOF
import json, sys
try:
    with open('${WEIGHTS_OUT}') as f:
        w = json.load(f)
    n = w.get('n_samples', '?')
    print(f'\nWGS calibrated weights (n={n} samples):')
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
echo " WGS Calibration Complete"
echo "================================================================="
echo " Weights:      ${WEIGHTS_OUT}"
echo " Accuracy:     ${TABLE_OUT}"
[[ -f "${COMPARE_OUT}" ]]    && echo " Strat(cal):   ${COMPARE_OUT}"
[[ -f "${COMPARE_RC_OUT}" ]] && echo " Strat(rc):    ${COMPARE_RC_OUT}"
echo ""
echo " Next: push updated weights to GitHub and rerun pipeline"
echo "   --weighting calibrated --weights_file conf/tool_weights_wgs_puhti_v1.json"
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
echo " Status:  bash scripts/run_wgs_calibration_puhti.sh --status"
echo " Logs:    ${LOGS_DIR}/"
echo ""
echo " Expected outputs:"
echo "   ${WEIGHTS_OUT}"
echo "   ${TABLE_OUT}"
echo "==================================================================="
