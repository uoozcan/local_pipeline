#!/bin/bash
#=============================================================================
# HLA Tool Weight Calibration — CSC Puhti
#=============================================================================
# Runs the full empirical calibration workflow on Puhti as dependent SLURM jobs:
#
#   Phase 0 (extract array):  Stream HLA-region FASTQs from EBI 30x CRAMs
#   Phase 1 (typing  array):  Run all HLA tools on 1KGP WGS FASTQ pairs
#   Phase 1b (collect job):   Symlink outputs into by_tool/ layout
#   Phase 2 (calibrate job):  Download GT + population file + calibrate + Wilcoxon
#
# Usage
# -----
#   # Most common: stream+type+calibrate 50 samples from EBI 30x CRAMs
#   bash submit_calibration_puhti.sh --project project_2008084
#
#   # Use pre-extracted FASTQs (skip Phase 0)
#   bash submit_calibration_puhti.sh \
#       --project project_2008084 \
#       --fastq-dir /scratch/project_2008084/1kgp_fastqs
#
#   # Skip typing entirely (only run calibration on existing results)
#   bash submit_calibration_puhti.sh \
#       --project project_2008084 \
#       --skip-typing
#
#   # Fetch ENA accessions first (needed once before first run)
#   bash submit_calibration_puhti.sh \
#       --project project_2008084 \
#       --fetch-accessions
#
# Output (after all jobs complete)
# ---------------------------------
#   conf/tool_weights_wgs_v2.json      calibrated weights (5 genes, 30x WGS)
#   conf/tool_accuracy_wgs_v2.tsv      accuracy table — Supplementary Table
#   conf/tool_accuracy_wgs_v2_*.tsv    per-population accuracy tables
#   conf/strategy_comparison_*.tsv     Wilcoxon test: equal vs calibrated
#
# Notes
# -----
#   - 1KGP 30x NYGC CRAMs use GRCh38 with UCSC chr names (chr6:...)
#   - Converting CRAMs to FASTQ first avoids chr naming issues inside containers
#   - EBI → Puhti download speed is typically fast (~100-200 MB/s)
#   - FASTQ extraction: ~2h per sample  |  Typing: ~6-8h per sample
#   - Phase 0 + Phase 1 throttled to 10 concurrent (--array=1-N%10)
#=============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

#-----------------------------------------------------------------------------
# Defaults
#-----------------------------------------------------------------------------
PROJECT_ID=""
INSTALL_DIR=""
FASTQ_DIR=""
N_SAMPLES=50
BATCH_SIZE=5            # samples per run (set lower to save scratch space)
TOOLS="hlahd,spechla,arcashla,optitype"
GENES="A,B,C,DRB1,DQB1"
REFERENCE="hg38"
RESOLUTION="2-field"
SOURCE="30x"            # 30x (NYGC GRCh38) or phase3 (4x hg19)
SKIP_TYPING=false
FETCH_ACCESSIONS=false
CLEANUP_WORK=true       # delete Nextflow work dir after collecting (saves ~10-20 GB/batch)
USER_GT_FILE=""         # path to user-provided GT CSV (overrides Gourraud 2014 download)
SPECHLA_SIF=""          # path to existing spechla SIF (optional; default: container_dir/spechla_with_spechap.sif)
HLA_REGION_HG38="chr6:28000000-34000000"
HLA_REGION_HG19="6:28000000-34000000"

#-----------------------------------------------------------------------------
# Parse arguments
#-----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $(basename "$0") --project <project_id> [options]

Required:
  --project ID          CSC project ID (e.g. project_2008084)

Options:
  --source SOURCE       30x (default, NYGC GRCh38) or phase3 (4x hg19)
  --fastq-dir DIR       Use pre-extracted FASTQs; skip Phase 0 extraction
  --n-samples N         Total 1KGP samples to calibrate (default: 50)
  --batch-size N        Samples to process per run (default: 5, saves scratch space)
                        Run script 10 times to process all 50 samples.
                        Each batch uses ~15-30 GB scratch; work dirs cleaned after.
  --tools TOOLS         Comma-separated tools (default: hlahd,spechla,arcashla,optitype)
  --genes GENES         Genes for calibration (default: A,B,C,DRB1,DQB1)
  --gt-file PATH        Use this CSV as ground truth (overrides downloading 1KGP Phase 1 GT).
                        Columns: sample,A1,A2,B1,B2[,...] (comma-separated, header required).
                        SAMPLE_LIST is built from samples in this file; genes auto-detected
                        from column headers. Ideal for custom/local GT datasets.
  --spechla-sif PATH    Path to an existing spechla Singularity SIF file.
                        Skips building/pulling a new one. Useful when the SIF was already
                        built during a previous pipeline installation.
                        Default: \${INSTALL_DIR}/containers/spechla_with_spechap.sif
  --skip-typing         Skip Phase 0+1 — run calibration on existing results
  --no-cleanup          Keep Nextflow work dirs after collection (for debugging)
  --fetch-accessions    Download ENA run report for all 1KGP 30x samples, then exit
  --install-dir DIR     Override install dir (default: /projappl/<project>/hla_typing)
  --help                Show this help

Space usage per batch (5 samples):
  FASTQs:       ~1-2 GB  (kept across batches, reused if already present)
  Nextflow work: ~10-20 GB (deleted after results are collected)
  Results:      ~0.5 GB  (kept permanently for calibration)
  Total peak:   ~15-25 GB per batch
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --project)          PROJECT_ID="$2";       shift 2 ;;
        --source)           SOURCE="$2";            shift 2 ;;
        --fastq-dir)        FASTQ_DIR="$2";         shift 2 ;;
        --n-samples)        N_SAMPLES="$2";         shift 2 ;;
        --tools)            TOOLS="$2";             shift 2 ;;
        --genes)            GENES="$2";             shift 2 ;;
        --gt-file)          USER_GT_FILE="$2";      shift 2 ;;
        --spechla-sif)      SPECHLA_SIF="$2";       shift 2 ;;
        --batch-size)       BATCH_SIZE="$2";        shift 2 ;;
        --no-cleanup)       CLEANUP_WORK=false;     shift ;;
        --skip-typing)      SKIP_TYPING=true;       shift ;;
        --fetch-accessions) FETCH_ACCESSIONS=true;  shift ;;
        --install-dir)      INSTALL_DIR="$2";       shift 2 ;;
        --help|-h)          usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

if [[ -z "$PROJECT_ID" ]]; then
    PROJECT_ID="${CSC_PROJECT:-}"
    [[ -z "$PROJECT_ID" ]] && echo "ERROR: --project is required" && usage
fi

#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
# Auto-detect install directory (try projappl first, then common scratch paths)
if [[ -z "$INSTALL_DIR" ]]; then
    for _candidate in \
        "/projappl/${PROJECT_ID}/hla_typing" \
        "/projappl/${PROJECT_ID}/hla_typing_pipeline" \
        "/scratch/${PROJECT_ID}/${USER}/hla_typing_pipeline" \
        "/scratch/${PROJECT_ID}/${USER}/local_pipeline/hla_typing_pipeline"; do
        if [[ -f "${_candidate}/main.nf" ]]; then
            INSTALL_DIR="$_candidate"
            break
        fi
    done
fi
INSTALL_DIR="${INSTALL_DIR:-/projappl/${PROJECT_ID}/hla_typing}"
if [[ ! -f "${INSTALL_DIR}/main.nf" ]]; then
    echo "[ERROR] Pipeline not found at: $INSTALL_DIR"
    echo "        Pass --install-dir /path/to/hla_typing_pipeline to override."
    exit 1
fi

SCRATCH_BASE="/scratch/${PROJECT_ID}/${USER}/hla_calibration"
FASTQ_DIR="${FASTQ_DIR:-${SCRATCH_BASE}/1kgp_fastqs}"
RESULTS_DIR="${SCRATCH_BASE}/1kgp_typing_results"
GT_FILE="${SCRATCH_BASE}/conf/1kgp_hla_gt.tsv"
POP_FILE="${SCRATCH_BASE}/conf/1kgp_populations.tsv"
ACC_FILE="${SCRATCH_BASE}/conf/1kgp_30x_accessions.tsv"
WEIGHTS_OUT="${INSTALL_DIR}/conf/tool_weights_wgs_v2.json"
TABLE_OUT="${INSTALL_DIR}/conf/tool_accuracy_wgs_v2.tsv"
LOGS_DIR="${SCRATCH_BASE}/logs"
SAMPLE_LIST="${SCRATCH_BASE}/conf/1kgp_sample_list.txt"

[[ "$SOURCE" == "30x" ]] && HLA_REGION="$HLA_REGION_HG38" || HLA_REGION="$HLA_REGION_HG19"

mkdir -p "$SCRATCH_BASE" "$FASTQ_DIR" "$RESULTS_DIR" "$LOGS_DIR" "${SCRATCH_BASE}/conf"

# Validate and stage user-provided GT file
if [[ -n "$USER_GT_FILE" ]]; then
    if [[ ! -f "$USER_GT_FILE" ]]; then
        echo "[ERROR] --gt-file not found: $USER_GT_FILE"
        exit 1
    fi
    GT_FILE="${SCRATCH_BASE}/conf/ground_truth_data.csv"
    cp "$USER_GT_FILE" "$GT_FILE"
    echo "[INFO] Using provided GT file: $USER_GT_FILE"
fi

echo "============================================================"
echo " HLA Calibration Workflow — CSC Puhti"
echo "============================================================"
echo " Project:    $PROJECT_ID"
echo " Install:    $INSTALL_DIR"
echo " Scratch:    $SCRATCH_BASE"
echo " FASTQ dir:  $FASTQ_DIR"
echo " Results:    $RESULTS_DIR"
echo " Source:     $SOURCE ($HLA_REGION)"
echo " Tools:      $TOOLS"
echo " Genes:      $GENES"
echo " N samples:  $N_SAMPLES (batch size: $BATCH_SIZE)"
echo " Cleanup:    $CLEANUP_WORK (delete work dirs after collect)"
[[ -n "$USER_GT_FILE" ]] && echo " GT file:    $USER_GT_FILE"
echo "============================================================"
echo ""

#-----------------------------------------------------------------------------
# Build sample list — either from user-provided GT CSV or default 50-sample list
#-----------------------------------------------------------------------------
if [[ -n "$USER_GT_FILE" ]]; then
    # Extract sample IDs from first column of provided GT CSV (skip header)
    python3 - << PYEOF
import csv, sys
with open('${GT_FILE}') as f:
    for i, row in enumerate(csv.reader(f)):
        if i == 0:
            continue  # skip header
        s = row[0].strip()
        if s:
            print(s)
PYEOF
    # Redirect output to SAMPLE_LIST
    python3 -c "
import csv
rows = []
with open('${GT_FILE}') as f:
    for i, row in enumerate(csv.reader(f)):
        if i == 0: continue
        s = row[0].strip()
        if s: rows.append(s)
with open('${SAMPLE_LIST}', 'w') as out:
    out.write('\n'.join(rows) + '\n')
"
    # Auto-detect genes from CSV header columns (e.g. A1,A2,B1,B2 → A,B)
    GENES=$(python3 -c "
import csv, re
with open('${GT_FILE}') as f:
    header = next(csv.reader(f))
seen, genes = set(), []
for col in header[1:]:
    m = re.match(r'^(DRB1|DQA1|DQB1|DPA1|DPB1|[A-Z])[_\d]', col.strip())
    if m and m.group(1) not in seen:
        genes.append(m.group(1)); seen.add(m.group(1))
print(','.join(genes))
")
    echo "[INFO] Genes auto-detected from GT file: $GENES"
else
    # Default: recommended 50-sample list spanning 5 superpopulations (1KGP Phase 1 GT)
    # NOTE: NA12878/NA12889 NOT in Phase 1 GT. YRI trio + CHB + GBR + TSI preferred.
    cat > "${SCRATCH_BASE}/conf/1kgp_recommended_samples.txt" << 'SAMPLES'
NA19238
NA19239
NA19240
NA19209
NA19210
NA19129
NA19130
NA19131
NA19152
NA19153
NA18526
NA18542
NA18524
NA18529
NA18532
NA18537
NA18561
NA18562
NA18563
NA18564
HG00096
HG00097
HG00099
HG00100
HG00101
HG00102
HG00103
HG00105
HG00106
HG00107
NA20502
NA20503
NA20504
NA20505
NA20506
NA20507
NA20508
NA20509
NA20510
NA20511
HG01565
HG01566
HG01567
HG01568
HG01570
HG01571
HG01572
HG01573
NA18939
NA18940
SAMPLES
    head -n "$N_SAMPLES" "${SCRATCH_BASE}/conf/1kgp_recommended_samples.txt" > "$SAMPLE_LIST"
fi

ACTUAL_N=$(wc -l < "$SAMPLE_LIST")
echo "[INFO] Sample list: $ACTUAL_N samples"

#-----------------------------------------------------------------------------
# --fetch-accessions mode: download ENA run report for PRJEB31736, then exit
#-----------------------------------------------------------------------------
if [[ "$FETCH_ACCESSIONS" == "true" ]]; then
    echo ""
    echo "=== Fetching ENA accessions for 1KGP 30x NYGC CRAMs (PRJEB31736) ==="

    # Download ground truth first (needed to check has_gt)
    if [[ ! -f "$GT_FILE" ]]; then
        echo "[INFO] Downloading 1KGP Phase 1 ground truth..."
        module load python-data 2>/dev/null || module load python/3.9 2>/dev/null || module load python
        python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" download-gt \
            --output "$GT_FILE"
    fi

    ENA_TSV="${SCRATCH_BASE}/conf/ena_run_report.tsv"
    echo "[INFO] Fetching ENA run report (this may take 1-2 minutes)..."
    curl -sL \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB31736&result=read_run&fields=run_accession,submitted_ftp&format=tsv" \
        > "$ENA_TSV"
    NROWS=$(tail -n +2 "$ENA_TSV" | wc -l || echo 0)
    echo "[INFO] Downloaded $NROWS run entries"

    python3 - << PYEOF
import csv, re, os

gt_file = "${GT_FILE}"
ena_tsv = "${ENA_TSV}"
out_tsv = "${ACC_FILE}"

gt = set()
if os.path.exists(gt_file):
    with open(gt_file) as f:
        next(f, None)
        for line in f:
            row = line.strip().split('\t')
            if row: gt.add(row[0].strip())
print(f"[INFO] GT samples: {len(gt)}")

written = 0
with open(ena_tsv) as f, open(out_tsv, 'w') as out:
    out.write("sample_id\terr_accession\thas_gt\n")
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ftp = row.get('submitted_ftp', '').strip()
        err = row.get('run_accession', '').strip()
        if not ftp or not err: continue
        m = re.search(r'/([A-Z0-9]+)\.final\.cram(?:;|\$)', ftp)
        if not m: continue
        s = m.group(1)
        has_gt = 'yes' if s in gt else 'no'
        out.write(f"{s}\t{err}\t{has_gt}\n")
        written += 1

print(f"[INFO] Wrote {written} accessions")
gt_count = sum(1 for ln in open(out_tsv).readlines()[1:] if ln.split('\t')[2].strip() == 'yes')
print(f"[INFO] Samples with ground truth: {gt_count}")
print(f"[INFO] Accession file: ${ACC_FILE}")
PYEOF

    echo ""
    echo "[DONE] Accessions saved to: $ACC_FILE"
    echo ""
    echo "Next: run the calibration:"
    echo "  bash $(basename "$0") --project $PROJECT_ID --source $SOURCE --n-samples $N_SAMPLES"
    exit 0
fi

#-----------------------------------------------------------------------------
# Load ENA accessions (needed for Phase 0 FASTQ extraction)
#-----------------------------------------------------------------------------
declare -A CRAM_ERR=()

# Hardcoded fallback — 7 samples already confirmed working locally
CRAM_ERR[NA19238]=ERR3239453
CRAM_ERR[NA19239]=ERR3239454
# NA19240: ERR3239455 returns 404 at EBI — run --fetch-accessions to get correct accession
CRAM_ERR[NA18526]=ERR3239353
CRAM_ERR[NA18542]=ERR3239356
CRAM_ERR[HG00096]=ERR3240114
CRAM_ERR[NA20502]=ERR3239785

if [[ -f "$ACC_FILE" ]]; then
    while IFS=$'\t' read -r sample err has_gt; do
        [[ "$has_gt" == "yes" ]] && CRAM_ERR["$sample"]="$err"
    done < <(tail -n +2 "$ACC_FILE")
    echo "[INFO] Loaded ${#CRAM_ERR[@]} samples with 30x accessions + ground truth"
else
    echo "[INFO] Using hardcoded 7-sample fallback (run --fetch-accessions for full 50+)"
fi

# Rebuild SAMPLE_LIST to only include GT-confirmed samples (present in CRAM_ERR)
# This ensures ACTUAL_N reflects real typing workload and avoids empty array tasks.
GT_CONFIRMED_LIST="${SCRATCH_BASE}/conf/1kgp_gt_confirmed_samples.txt"
: > "$GT_CONFIRMED_LIST"
while IFS= read -r S; do
    [[ -n "${CRAM_ERR[$S]:-}" ]] && echo "$S" >> "$GT_CONFIRMED_LIST"
done < "$SAMPLE_LIST"

GT_N=$(wc -l < "$GT_CONFIRMED_LIST")
CANDIDATE_N=$ACTUAL_N
SKIPPED=$(( CANDIDATE_N - GT_N ))

if [[ $SKIPPED -gt 0 ]]; then
    echo "[INFO] $SKIPPED candidate sample(s) have no 30x accession with ground truth — excluded"
fi
echo "[INFO] GT-confirmed samples to type: ${GT_N} / ${CANDIDATE_N}"

cp "$GT_CONFIRMED_LIST" "$SAMPLE_LIST"
ACTUAL_N=$GT_N

if [[ $ACTUAL_N -eq 0 ]]; then
    echo "[ERROR] No GT-confirmed samples found."
    echo "        Run --fetch-accessions first to populate the accession table, then retry."
    exit 1
fi

#-----------------------------------------------------------------------------
# Batch filtering: trim SAMPLE_LIST to the next BATCH_SIZE untyped samples
# (samples that are missing at least one tool's result file)
#-----------------------------------------------------------------------------
if [[ "$SKIP_TYPING" != "true" ]]; then
    BATCH_LIST="${SCRATCH_BASE}/conf/batch_current.txt"
    : > "$BATCH_LIST"
    ALREADY_COMPLETE=0
    while IFS= read -r _S; do
        _ALL_DONE=true
        for _T in $(echo "$TOOLS" | tr ',' ' '); do
            [[ ! -f "${RESULTS_DIR}/${_S}/${_T}/${_S}_${_T}.txt" ]] && _ALL_DONE=false && break
        done
        if [[ "$_ALL_DONE" == "true" ]]; then
            ALREADY_COMPLETE=$(( ALREADY_COMPLETE + 1 ))
        elif [[ $(wc -l < "$BATCH_LIST") -lt $BATCH_SIZE ]]; then
            echo "$_S" >> "$BATCH_LIST"
        fi
    done < "$SAMPLE_LIST"

    BATCH_N=$(wc -l < "$BATCH_LIST")
    TOTAL_REMAINING=$(( ACTUAL_N - ALREADY_COMPLETE ))
    BATCHES_LEFT=$(( (TOTAL_REMAINING + BATCH_SIZE - 1) / BATCH_SIZE ))

    echo "[INFO] Already typed:   $ALREADY_COMPLETE / $ACTUAL_N"
    echo "[INFO] Remaining:       $TOTAL_REMAINING samples (~$BATCHES_LEFT batch(es) of $BATCH_SIZE)"
    echo "[INFO] This batch:      $BATCH_N sample(s)"

    if [[ $BATCH_N -eq 0 ]]; then
        echo ""
        echo "[DONE] All $ACTUAL_N GT-confirmed samples are already typed."
        echo "       Run with --skip-typing to re-run calibration only."
        exit 0
    fi

    cp "$BATCH_LIST" "$SAMPLE_LIST"
    ACTUAL_N=$BATCH_N
fi

#-----------------------------------------------------------------------------
# Phase 0: FASTQ extraction array job (stream HLA region from EBI 30x CRAMs)
#-----------------------------------------------------------------------------
PHASE0_JOB=""
PHASE0_DEP=""

if [[ "$SKIP_TYPING" != "true" ]] && [[ "$SOURCE" == "30x" ]]; then
    # Build extraction sample list (only samples with accession and missing FASTQ)
    EXTRACT_LIST="${SCRATCH_BASE}/conf/extract_sample_list.txt"
    : > "$EXTRACT_LIST"
    while IFS= read -r SAMPLE; do
        R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
        if [[ -f "$R1" ]] && [[ -f "$R2" ]]; then
            echo "[SKIP] $SAMPLE — FASTQ already present"
        elif [[ -n "${CRAM_ERR[$SAMPLE]:-}" ]]; then
            echo "$SAMPLE" >> "$EXTRACT_LIST"
        else
            echo "[WARN] $SAMPLE — no ENA accession (run --fetch-accessions)"
        fi
    done < "$SAMPLE_LIST"

    EXTRACT_N=$(wc -l < "$EXTRACT_LIST")

    if [[ "$EXTRACT_N" -gt 0 ]]; then
        echo ""
        echo "=== Phase 0: Submitting FASTQ extraction array ($EXTRACT_N samples) ==="

        PHASE0_SCRIPT="${SCRATCH_BASE}/phase0_extract.sh"
        cat > "$PHASE0_SCRIPT" << PHASE0EOF
#!/bin/bash
#SBATCH --job-name=1kgp_extract
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/extract_%A_%a.out
#SBATCH --error=${LOGS_DIR}/extract_%A_%a.err

set -euo pipefail

SAMPLE_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${EXTRACT_LIST}")
[[ -z "\$SAMPLE_ID" ]] && echo "No sample for task \${SLURM_ARRAY_TASK_ID}" && exit 1

echo "--- Task \${SLURM_ARRAY_TASK_ID}: \$SAMPLE_ID ---"
echo "Date: \$(date)"

# Load samtools (biokit provides htslib + CRAM support with EBI MD5 reference cache)
module load biokit 2>/dev/null || module load samtools 2>/dev/null || \
    { echo "ERROR: samtools not available"; exit 1; }

# CRAM reference cache (EBI hosted MD5 refs)
export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${SCRATCH_BASE}/ref_cache/%2s/%2s/%s"
mkdir -p "${SCRATCH_BASE}/ref_cache"

# Map sample → ERR accession (embedded from CRAM_ERR map at generation time)
declare -A CRAM_ERR
$(for s in "${!CRAM_ERR[@]}"; do echo "CRAM_ERR[$s]=${CRAM_ERR[$s]}"; done)

ERR="\${CRAM_ERR[\$SAMPLE_ID]:-}"
if [[ -z "\$ERR" ]]; then
    echo "ERROR: No ERR accession for \$SAMPLE_ID" && exit 1
fi

PREFIX3="\${ERR:0:6}"
CRAM_URL="ftp://ftp.sra.ebi.ac.uk/vol1/run/\${PREFIX3}/\${ERR}/\${SAMPLE_ID}.final.cram"
echo "CRAM URL: \$CRAM_URL"
echo "HLA region: ${HLA_REGION_HG38}"

TMPBAM="${FASTQ_DIR}/\${SAMPLE_ID}_hla_tmp.bam"
NAMESORT="${FASTQ_DIR}/\${SAMPLE_ID}_namesorted.bam"
R1="${FASTQ_DIR}/\${SAMPLE_ID}_R1.fastq.gz"
R2="${FASTQ_DIR}/\${SAMPLE_ID}_R2.fastq.gz"

mkdir -p "${FASTQ_DIR}"

# Skip if FASTQ already exists
if [[ -f "\$R1" ]] && [[ -f "\$R2" ]]; then
    echo "[SKIP] FASTQ already present: \$R1"
    exit 0
fi

# Step 1: stream HLA region BAM from EBI CRAM
echo "Streaming HLA region from EBI..."
samtools view -b -@ 3 -o "\${TMPBAM}.tmp" "\$CRAM_URL" "${HLA_REGION_HG38}" \
    2>"${LOGS_DIR}/extract_\${SAMPLE_ID}_samtools.log" \
    || { echo "ERROR: samtools view failed"; cat "${LOGS_DIR}/extract_\${SAMPLE_ID}_samtools.log"; exit 1; }
mv "\${TMPBAM}.tmp" "\$TMPBAM"

READ_COUNT=\$(samtools view -c "\$TMPBAM" 2>/dev/null || echo "?")
SIZE=\$(du -sh "\$TMPBAM" | cut -f1)
echo "[OK] HLA BAM: \$SIZE, \$READ_COUNT reads"

if [[ "\$READ_COUNT" != "?" ]] && [[ "\$READ_COUNT" -lt 5000 ]]; then
    echo "[WARN] Only \$READ_COUNT reads — this sample may have low HLA coverage"
fi

# Step 2: name-sort (required for samtools fastq to produce proper pairs)
echo "Name-sorting..."
samtools sort -n -@ 3 "\$TMPBAM" -o "\$NAMESORT" \
    2>>"${LOGS_DIR}/extract_\${SAMPLE_ID}_samtools.log"
rm -f "\$TMPBAM"

# Step 3: convert to paired FASTQ
echo "Converting to FASTQ..."
samtools fastq -@ 3 \
    -1 "\$R1" -2 "\$R2" \
    -0 /dev/null -s /dev/null \
    "\$NAMESORT" \
    2>>"${LOGS_DIR}/extract_\${SAMPLE_ID}_samtools.log"
rm -f "\$NAMESORT"

R1_READS=\$(zcat "\$R1" | awk 'NR%4==1' | wc -l || echo "?")
R1_SIZE=\$(du -sh "\$R1" | cut -f1)
echo "[OK] \${SAMPLE_ID}: R1=\${R1_SIZE} (\${R1_READS} reads), R2=\$(du -sh \$R2 | cut -f1)"
PHASE0EOF

        PHASE0_JOB=$(sbatch --parsable \
            --array=1-${EXTRACT_N}%10 \
            "$PHASE0_SCRIPT")
        echo "[INFO] Phase 0 array job submitted: ${PHASE0_JOB} (array 1-${EXTRACT_N}%10)"
        PHASE0_DEP="--dependency=afterany:${PHASE0_JOB}"
    else
        echo "[INFO] All FASTQs already present — skipping Phase 0"
    fi
fi

#-----------------------------------------------------------------------------
# Phase 1: Typing array job (FASTQ → HLA calls)
#-----------------------------------------------------------------------------
TYPING_JOB=""

if [[ "$SKIP_TYPING" != "true" ]]; then
    echo ""
    echo "=== Phase 1: Submitting HLA typing array ==="

    TYPING_ARRAY_SCRIPT="${SCRATCH_BASE}/phase1_typing_array.sh"
    cat > "$TYPING_ARRAY_SCRIPT" << ARRAYEOF
#!/bin/bash
#SBATCH --job-name=1kgp_hla_typing
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --output=${LOGS_DIR}/typing_%A_%a.out
#SBATCH --error=${LOGS_DIR}/typing_%A_%a.err

set -euo pipefail

SAMPLE_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
[[ -z "\$SAMPLE_ID" ]] && echo "No sample for task \${SLURM_ARRAY_TASK_ID}" && exit 1
echo "--- Task \${SLURM_ARRAY_TASK_ID}: \$SAMPLE_ID ---"

R1="${FASTQ_DIR}/\${SAMPLE_ID}_R1.fastq.gz"
R2="${FASTQ_DIR}/\${SAMPLE_ID}_R2.fastq.gz"

if [[ ! -f "\$R1" ]] || [[ ! -f "\$R2" ]]; then
    echo "[WARN] FASTQ not found for \$SAMPLE_ID — skipping"
    exit 0
fi

# Check if all tools already typed for this sample
ALREADY_DONE=true
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    [[ ! -f "${RESULTS_DIR}/\${SAMPLE_ID}/\${TOOL}/\${SAMPLE_ID}_\${TOOL}.txt" ]] && ALREADY_DONE=false && break
done
if [[ "\$ALREADY_DONE" == "true" ]]; then
    echo "[SKIP] \$SAMPLE_ID — all tools already typed"
    exit 0
fi

module purge
module load nextflow/23.10.0 2>/dev/null || module load nextflow
module load singularity 2>/dev/null || true

export SINGULARITY_CACHEDIR="${INSTALL_DIR}/singularity_cache"
export NXF_SINGULARITY_CACHEDIR="\$SINGULARITY_CACHEDIR"
export NXF_ANSI_LOG=false
# Per-sample NXF_HOME prevents session lock contention between array tasks
export NXF_HOME="${SCRATCH_BASE}/nxf_home/\${SAMPLE_ID}"
mkdir -p "\$NXF_HOME"

# Write a per-sample samplesheet (FASTQ mode — avoids chr naming issues in containers)
SAMPLE_CSV="${SCRATCH_BASE}/samplesheets/\${SAMPLE_ID}.csv"
mkdir -p "${SCRATCH_BASE}/samplesheets"
echo "sample_id,fastq_1,fastq_2" > "\$SAMPLE_CSV"
echo "\${SAMPLE_ID},\${R1},\${R2}" >> "\$SAMPLE_CSV"

WORK_DIR="${SCRATCH_BASE}/work/\${SAMPLE_ID}_\${SLURM_JOB_ID}"
mkdir -p "\$WORK_DIR"

# cd to pipeline dir so Nextflow auto-loads nextflow.config (params defaults)
cd "${INSTALL_DIR}"

nextflow run main.nf \
    --input_samplesheet "\$SAMPLE_CSV" \
    --outdir "${RESULTS_DIR}" \
    --tools "${TOOLS}" \
    --seq_type dna \
    --reference "${REFERENCE}" \
    --resolution "${RESOLUTION}" \
    --max_cpus 20 \
    --max_memory "60GB" \
    --output_format text \
    --hla_genes 'HLA-A,HLA-B,HLA-C,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DPA1,HLA-DPB1' \
    --min_hla_reads 1000 \
    --min_read_length 50 \
    --min_tools 1 \
    --install_dir "${INSTALL_DIR}" \
    ${SPECHLA_SIF:+--spechla_sif "${SPECHLA_SIF}"} \
    -profile singularity \
    -work-dir "\$WORK_DIR" \
    -c "${INSTALL_DIR}/conf/puhti.config" 2>&1

echo "[DONE] \$SAMPLE_ID"
ARRAYEOF

    TYPING_JOB=$(sbatch --parsable \
        ${PHASE0_DEP} \
        --array=1-${ACTUAL_N}%10 \
        "$TYPING_ARRAY_SCRIPT")
    echo "[INFO] Phase 1 array job submitted: ${TYPING_JOB} (array 1-${ACTUAL_N}%10)"
    COLLECT_DEP="--dependency=afterany:${TYPING_JOB}"
else
    echo "[INFO] Skipping Phase 0+1 (--skip-typing)"
    COLLECT_DEP=""
fi

#-----------------------------------------------------------------------------
# Phase 1b: Collect results into by_tool/ layout
#-----------------------------------------------------------------------------
COLLECT_JOB=""
COLLECT_SCRIPT="${SCRATCH_BASE}/phase1b_collect.sh"
cat > "$COLLECT_SCRIPT" << COLLECTEOF
#!/bin/bash
#SBATCH --job-name=1kgp_collect
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/collect_%j.out
#SBATCH --error=${LOGS_DIR}/collect_%j.err

set -euo pipefail
echo "=== Phase 1b: Collecting results ==="
echo "Date: \$(date)"

for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    mkdir -p "${RESULTS_DIR}/by_tool/\${TOOL}"
done

TOTAL=0
while IFS= read -r SAMPLE; do
    for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
        SRC="${RESULTS_DIR}/\${SAMPLE}/\${TOOL}/\${SAMPLE}_\${TOOL}.txt"
        DST="${RESULTS_DIR}/by_tool/\${TOOL}/\${SAMPLE}_\${TOOL}.txt"
        if [[ -f "\$SRC" ]] && [[ ! -e "\$DST" ]]; then
            ln -sf "\$SRC" "\$DST"
            TOTAL=\$((TOTAL+1))
        fi
    done
done < "${SAMPLE_LIST}"

echo ""
echo "Symlinks created: \$TOTAL"
echo ""
echo "Results by tool:"
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    N=\$(find "${RESULTS_DIR}/by_tool/\${TOOL}/" -maxdepth 1 -name "*_\${TOOL}.txt" 2>/dev/null | wc -l)
    printf "  %-12s %3d samples\n" "\$TOOL" "\$N"
done

# Clean up Nextflow work directories to free scratch space (enabled by default)
if [[ "${CLEANUP_WORK}" == "true" ]]; then
    echo ""
    echo "Cleaning up Nextflow work directories..."
    WORK_BASE="${SCRATCH_BASE}/work"
    if [[ -d "\$WORK_BASE" ]]; then
        BEFORE=\$(du -sh "\$WORK_BASE" 2>/dev/null | cut -f1 || echo "?")
        rm -rf "\${WORK_BASE:?}"/*
        echo "[OK] Removed work dirs (\${BEFORE} freed); results preserved in ${RESULTS_DIR}"
    fi
else
    echo ""
    echo "[INFO] Work dirs kept for debugging: ${SCRATCH_BASE}/work/"
fi

echo ""
echo "Remaining scratch usage:"
du -sh "${SCRATCH_BASE}" 2>/dev/null || true
COLLECTEOF

COLLECT_JOB=$(sbatch --parsable \
    ${COLLECT_DEP} \
    "$COLLECT_SCRIPT")
echo "[INFO] Phase 1b collect job submitted: ${COLLECT_JOB}"

#-----------------------------------------------------------------------------
# Phase 2: Calibration + population stratification + Wilcoxon comparison
#-----------------------------------------------------------------------------
echo ""
echo "=== Phase 2: Submitting calibration job ==="

CALIBRATE_SCRIPT="${SCRATCH_BASE}/phase2_calibrate.sh"
COMPARE_OUT="${SCRATCH_BASE}/conf/strategy_comparison_calibrated.tsv"
COMPARE_RC_OUT="${SCRATCH_BASE}/conf/strategy_comparison_rc.tsv"

cat > "$CALIBRATE_SCRIPT" << CALEOF
#!/bin/bash
#SBATCH --job-name=hla_calibrate
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/calibrate_%j.out
#SBATCH --error=${LOGS_DIR}/calibrate_%j.err

set -euo pipefail

echo "=== Phase 2: HLA Tool Weight Calibration ==="
echo "Date: \$(date)"
echo ""

module purge
module load python-data 2>/dev/null || \
    module load python/3.9  2>/dev/null || \
    module load python

# Step 1: Download 1KGP Phase 1 ground truth (Gourraud et al. 2014)
if [[ ! -f "${GT_FILE}" ]]; then
    echo "[INFO] Downloading 1KGP Phase 1 ground truth..."
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" download-gt \
        --output "${GT_FILE}"
else
    echo "[INFO] Ground truth: ${GT_FILE}"
fi

# Step 2: Generate population file (sample → superpopulation TSV)
if [[ ! -f "${POP_FILE}" ]]; then
    echo ""
    echo "[INFO] Downloading 1KGP population panel..."
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" generate-population-file \
        --output "${POP_FILE}" \
        && echo "[OK] Population file: ${POP_FILE}" \
        || echo "[WARN] Could not download population file — continuing without stratification"
else
    echo "[INFO] Population file: ${POP_FILE}"
fi

# Step 3: Run calibration (5 genes, with per-population tables)
echo ""
echo "[INFO] Running calibration (genes: ${GENES})..."

POP_FILE_ARG=""
[[ -f "${POP_FILE}" ]] && POP_FILE_ARG="--population-file ${POP_FILE}"

# shellcheck disable=SC2086
python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" calibrate \
    --ground-truth "${GT_FILE}" \
    --results-dir "${RESULTS_DIR}/by_tool" \
    --data-type wgs \
    --genes "${GENES}" \
    --resolution "${RESOLUTION}" \
    --output-weights "${WEIGHTS_OUT}" \
    --output-table "${TABLE_OUT}" \
    \$POP_FILE_ARG

echo ""
echo "[OK] Calibration outputs:"
echo "     Weights: ${WEIGHTS_OUT}"
echo "     Table:   ${TABLE_OUT}"
[[ -f "${POP_FILE}" ]] && echo "     Per-pop: \$(dirname ${TABLE_OUT})/tool_accuracy_wgs_v2_*.tsv"

# Print weight summary
python3 - << PYEOF
import json, sys
try:
    with open('${WEIGHTS_OUT}') as f:
        w = json.load(f)
    n = w.get('n_samples', '?')
    print(f'\nCalibrated weights (n={n} samples):')
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

# Step 4: Compare voting strategies (equal vs calibrated) — Wilcoxon signed-rank
echo ""
echo "[INFO] Comparing strategies: equal vs calibrated (Wilcoxon test)..."
if [[ -f "${WEIGHTS_OUT}" ]]; then
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" compare-strategies \
        --ground-truth "${GT_FILE}" \
        --results-dir "${RESULTS_DIR}/by_tool" \
        --genes "${GENES}" \
        --resolution "${RESOLUTION}" \
        --weights-file "${WEIGHTS_OUT}" \
        --ref-mode equal \
        --test-mode calibrated \
        --output "${COMPARE_OUT}" \
        && echo "[OK] Strategy comparison: ${COMPARE_OUT}" \
        || echo "[WARN] Strategy comparison failed (scipy may not be available)"

    # Also compare equal vs read_confidence
    echo ""
    echo "[INFO] Comparing strategies: equal vs read_confidence..."
    python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" compare-strategies \
        --ground-truth "${GT_FILE}" \
        --results-dir "${RESULTS_DIR}/by_tool" \
        --genes "${GENES}" \
        --resolution "${RESOLUTION}" \
        --ref-mode equal \
        --test-mode read_confidence \
        --output "${COMPARE_RC_OUT}" \
        && echo "[OK] RC comparison: ${COMPARE_RC_OUT}" \
        || echo "[WARN] RC comparison failed"
else
    echo "[WARN] Weights file not found — skipping strategy comparison"
fi

echo ""
echo "============================================================"
echo " Phase 2 complete"
echo "============================================================"
echo " Weights:           ${WEIGHTS_OUT}"
echo " Accuracy table:    ${TABLE_OUT}"
[[ -f "${COMPARE_OUT}" ]] && echo " Strategy (equal vs cal): ${COMPARE_OUT}"
[[ -f "${COMPARE_RC_OUT}" ]] && echo " Strategy (equal vs rc):  ${COMPARE_RC_OUT}"
echo ""
echo " To use calibrated weights in pipeline runs:"
echo "   nextflow run ${INSTALL_DIR}/main.nf \\"
echo "       --weighting calibrated \\"
echo "       --weights-file ${WEIGHTS_OUT} \\"
echo "       ..."
echo "============================================================"
CALEOF

CALIBRATE_JOB=$(sbatch --parsable \
    --dependency=afterok:${COLLECT_JOB} \
    "$CALIBRATE_SCRIPT")
echo "[INFO] Phase 2 calibration job submitted: ${CALIBRATE_JOB}"

#-----------------------------------------------------------------------------
# Final summary
#-----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo " All jobs submitted"
echo "============================================================"
[[ -n "$PHASE0_JOB" ]] && echo " Phase 0 (extract):  ${PHASE0_JOB} (array 1-${EXTRACT_N:-?}%10)"
[[ -n "$TYPING_JOB" ]] && echo " Phase 1 (typing):   ${TYPING_JOB} (array 1-${ACTUAL_N}%10)"
echo " Phase 1b (collect): ${COLLECT_JOB}"
echo " Phase 2 (calibrate): ${CALIBRATE_JOB}"
echo ""
echo " Monitor:"
echo "   watch squeue -u \$USER"
echo "   tail -f ${LOGS_DIR}/calibrate_*.out"
echo ""
echo " Expected outputs after completion:"
echo "   ${WEIGHTS_OUT}"
echo "   ${TABLE_OUT}"
echo "   ${COMPARE_OUT}"
echo ""
echo " Next batch (auto-skips already-typed samples):"
echo "   bash $(basename "$0") --project ${PROJECT_ID}"
echo ""
echo " Re-run calibration only (all typing done):"
echo "   bash $(basename "$0") --project ${PROJECT_ID} --skip-typing"
echo ""
echo " Check Puhti billing usage:"
echo "   csc-projects"
echo "============================================================"
