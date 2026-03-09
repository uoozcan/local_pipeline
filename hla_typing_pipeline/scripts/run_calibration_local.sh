#!/bin/bash
#=============================================================================
# HLA Tool Weight Calibration — Local (laptop/workstation)
#=============================================================================
# Runs the empirical calibration workflow locally (no SLURM).
#
# Two data sources are supported (--source flag):
#
#   phase3 (default): 1KGP Phase 3 low-coverage (~4x) BAMs from EBI.
#       Fast to stream (~10 min/sample), but too low coverage for B gene.
#       HLA region uses ENSEMBL chr names: "6:28000000-34000000"
#
#   30x: 1KGP NYGC 30x high-coverage CRAMs from EBI SRA (GRCh38).
#       ~15-30 min to stream; reliable for all 8 genes including Class II.
#       Requires REF_PATH for CRAM MD5 reference cache.
#       HLA region uses UCSC chr names: "chr6:28000000-34000000"
#
# Usage
# -----
#   # Fetch ENA accessions for all 1KGP 30x samples with GT (~1 min)
#   bash scripts/run_calibration_local.sh --source 30x --fetch-accessions \
#       --work-dir /home/umut/hla_calibration
#
#   # Check current status
#   bash scripts/run_calibration_local.sh --source 30x --status \
#       --work-dir /home/umut/hla_calibration
#
#   # Run batch (extract FASTQs + type with all tools; repeat until 50 samples)
#   bash scripts/run_calibration_local.sh --source 30x \
#       --n-samples 15 \
#       --tools "hlahd,spechla,optitype,arcashla" \
#       --work-dir /home/umut/hla_calibration \
#       --parallel-extract 2 \
#       --target 50
#
# Estimated runtime per sample:
#   optitype:  ~30 min
#   arcashla:  ~60 min
#   spechla:   ~90 min
#   hlahd:     ~60 min
#   xhla:      ~60-120 min (align.pl, faster on small HLA BAM)
#   Total:     ~5-6 hours per sample, sequential
#=============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

#-----------------------------------------------------------------------------
# Defaults — edit these or pass via arguments
#-----------------------------------------------------------------------------
N_SAMPLES=5
TOOLS="hlahd,spechla,arcashla,optitype,xhla"
WORK_DIR="/mnt/d/hla_calibration"
SOURCE="phase3"   # phase3 | 30x
PARALLEL_EXTRACT=2
TARGET_N=50
STATUS_ONLY=false
FETCH_ACCESSIONS=false
USER_GT_FILE=""         # path to user-provided GT CSV (overrides Gourraud 2014 download)
PREPARE_BATCHES=false
BATCH_SIZE=5
BATCH_FILE=""

# Use the existing run.config — it already has correct container paths,
# samtools wrapper binds, memory limits, and errorStrategy for all tools
RUN_CONFIG="${PIPELINE_DIR}/run.config"

# EBI base URL for 1KGP phase 3 low-coverage alignments (HTTPS)
EBI_BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"

# 30x CRAM ERR accessions (NYGC, GRCh38, 2019) — hardcoded fallback (7 samples)
declare -A CRAM_ERR=(
    [NA19238]=ERR3239453 [NA19239]=ERR3239454 [NA19240]=ERR3239455
    [NA18526]=ERR3239353 [NA18542]=ERR3239356
    [HG00096]=ERR3240114 [NA20502]=ERR3239785
)

#-----------------------------------------------------------------------------
# Parse arguments
#-----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: $(basename "$0") [options]

Options:
  --source SOURCE         Data source: phase3 (4x, default) or 30x (NYGC high-cov)
  --n-samples N           Max new samples to process per run (default: 5)
  --tools TOOLS           Comma-separated tools (default: hlahd,spechla,arcashla,optitype,xhla)
  --work-dir DIR          Working directory (default: /mnt/d/hla_calibration)
  --run-config FILE       Nextflow config file (default: run.config)
  --parallel-extract N    Parallel FASTQ extraction jobs (default: 2)
  --target N              Stop when N samples typed per core tool (default: 50)
  --gt-file PATH          Use this CSV as ground truth (overrides downloading 1KGP Phase 1 GT).
                          Columns: sample,A1,A2,B1,B2[,...] (comma-separated, header required).
                          SAMPLES list is built from this file; genes auto-detected from headers.
  --fetch-accessions      Fetch ENA run report and build accession table, then exit
  --status                Show typing progress per tool, then exit
  --prepare-batches       Pre-generate batch_NNN.csv files in WORK_DIR/batches/, then exit
  --batch-size N          Samples per batch for --prepare-batches (default: 5)
  --batch-file FILE       Run a specific pre-prepared batch CSV (overrides --n-samples)
  --skip-extract          Skip FASTQ extraction (use existing FASTQs)
  --skip-typing           Skip Nextflow typing step
  --help                  Show this help
EOF
    exit 0
}

SKIP_EXTRACT=false
SKIP_TYPING=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --source)           SOURCE="$2";           shift 2 ;;
        --n-samples)        N_SAMPLES="$2";        shift 2 ;;
        --tools)            TOOLS="$2";            shift 2 ;;
        --work-dir)         WORK_DIR="$2";         shift 2 ;;
        --run-config)       RUN_CONFIG="$2";       shift 2 ;;
        --parallel-extract) PARALLEL_EXTRACT="$2"; shift 2 ;;
        --target)           TARGET_N="$2";         shift 2 ;;
        --gt-file)          USER_GT_FILE="$2";         shift 2 ;;
        --fetch-accessions) FETCH_ACCESSIONS=true;   shift ;;
        --status)           STATUS_ONLY=true;        shift ;;
        --prepare-batches)  PREPARE_BATCHES=true;    shift ;;
        --batch-size)       BATCH_SIZE="$2";         shift 2 ;;
        --batch-file)       BATCH_FILE="$2";         shift 2 ;;
        --skip-extract)     SKIP_EXTRACT=true;       shift ;;
        --skip-typing)      SKIP_TYPING=true;        shift ;;
        --help|-h)          usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

#-----------------------------------------------------------------------------
# Paths — separate result dirs for 30x vs phase3 to avoid mixing
#-----------------------------------------------------------------------------
HLA_BAM_DIR="${WORK_DIR}/hla_bams"
if [[ "$SOURCE" == "30x" ]]; then
    RESULTS_DIR="${WORK_DIR}/typing_results_30x"
    HLA_REGION="chr6:28000000-34000000"   # GRCh38 (UCSC chr prefix)
else
    RESULTS_DIR="${WORK_DIR}/typing_results"
    HLA_REGION="6:28000000-34000000"      # hg19 (ENSEMBL, no prefix)
fi
GT_FILE="${WORK_DIR}/1kgp_hla_gt.tsv"
LOGS_DIR="${WORK_DIR}/logs"

# User-provided GT file overrides the downloaded Gourraud 2014 GT
if [[ -n "$USER_GT_FILE" ]]; then
    if [[ ! -f "$USER_GT_FILE" ]]; then
        echo "[ERROR] --gt-file not found: $USER_GT_FILE"
        exit 1
    fi
    GT_FILE="$USER_GT_FILE"
    echo "[INFO] Using provided GT file: $GT_FILE"
fi

if [[ "$SOURCE" == "30x" ]]; then
    WEIGHTS_OUT="${PIPELINE_DIR}/conf/tool_weights_wgs_v2.json"
    TABLE_OUT="${PIPELINE_DIR}/conf/tool_accuracy_wgs_v2.tsv"
    GENES="A,B,C,DRB1,DQB1"
else
    WEIGHTS_OUT="${PIPELINE_DIR}/conf/tool_weights_wgs.json"
    TABLE_OUT="${PIPELINE_DIR}/conf/tool_accuracy_wgs.tsv"
    GENES="A,B,C"
fi

# Auto-detect genes from user-provided GT CSV header (overrides SOURCE-based default)
if [[ -n "$USER_GT_FILE" ]]; then
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
fi

mkdir -p "$HLA_BAM_DIR" "$RESULTS_DIR" "$LOGS_DIR" "${WORK_DIR}/conf"

#-----------------------------------------------------------------------------
# Ordered sample list — priority: diverse populations, all with 1KGP GT
# Populated from known 1KGP 30x NYGC CRAM accessions (PRJEB31736)
#-----------------------------------------------------------------------------
ALL_SAMPLES_ORDERED=(
    # YRI (West African)
    NA19238 NA19239 NA19240 NA19209 NA19210 NA19129 NA19130 NA19131 NA19152
    NA19153 NA19159 NA19171 NA19172 NA19190 NA19200 NA19201 NA19204 NA19207
    # CHB (East Asian)
    NA18526 NA18542 NA18524 NA18529 NA18532 NA18537 NA18561 NA18562 NA18563
    NA18564 NA18566 NA18567 NA18570 NA18571 NA18573 NA18574 NA18576 NA18577
    # GBR (British)
    HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00105 HG00106
    HG00107 HG00108 HG00109 HG00110 HG00111 HG00112 HG00113 HG00114 HG00116
    # TSI (Italian)
    NA20502 NA20503 NA20504 NA20505 NA20506 NA20507 NA20508 NA20509 NA20510
    NA20511 NA20512 NA20513 NA20514 NA20515 NA20516 NA20517 NA20518 NA20519
    # CEU (European-American)
    NA12878 NA12889 NA12891 NA12892 NA06984 NA06985 NA06986 NA06989 NA06994
    # PEL (Peruvian)
    HG01565 HG01566 HG01567 HG01568 HG01570 HG01571 HG01572 HG01573 HG01574
    # MXL (Mexican-American)
    NA19648 NA19649 NA19651 NA19652 NA19654 NA19655 NA19657 NA19658 NA19659
    # PUR (Puerto Rican)
    HG00731 HG00732 HG00733 HG00734 HG00736 HG00737 HG00738 HG00739 HG00740
    # JPT (Japanese)
    NA18939 NA18940 NA18941 NA18942 NA18943 NA18944 NA18945 NA18946 NA18947
)

declare -A SAMPLE_POP=(
    [NA19238]=YRI [NA19239]=YRI [NA19240]=YRI [NA19209]=YRI [NA19210]=YRI
    [NA19129]=YRI [NA19130]=YRI [NA19131]=YRI [NA19152]=YRI [NA19153]=YRI
    [NA19159]=YRI [NA19171]=YRI [NA19172]=YRI [NA19190]=YRI [NA19200]=YRI
    [NA19201]=YRI [NA19204]=YRI [NA19207]=YRI
    [NA18526]=CHB [NA18542]=CHB [NA18524]=CHB [NA18529]=CHB [NA18532]=CHB
    [NA18537]=CHB [NA18561]=CHB [NA18562]=CHB [NA18563]=CHB [NA18564]=CHB
    [NA18566]=CHB [NA18567]=CHB [NA18570]=CHB [NA18571]=CHB [NA18573]=CHB
    [NA18574]=CHB [NA18576]=CHB [NA18577]=CHB
    [HG00096]=GBR [HG00097]=GBR [HG00099]=GBR [HG00100]=GBR [HG00101]=GBR
    [HG00102]=GBR [HG00103]=GBR [HG00105]=GBR [HG00106]=GBR [HG00107]=GBR
    [HG00108]=GBR [HG00109]=GBR [HG00110]=GBR [HG00111]=GBR [HG00112]=GBR
    [HG00113]=GBR [HG00114]=GBR [HG00116]=GBR
    [NA20502]=TSI [NA20503]=TSI [NA20504]=TSI [NA20505]=TSI [NA20506]=TSI
    [NA20507]=TSI [NA20508]=TSI [NA20509]=TSI [NA20510]=TSI [NA20511]=TSI
    [NA20512]=TSI [NA20513]=TSI [NA20514]=TSI [NA20515]=TSI [NA20516]=TSI
    [NA20517]=TSI [NA20518]=TSI [NA20519]=TSI
    [NA12878]=CEU [NA12889]=CEU [NA12891]=CEU [NA12892]=CEU [NA06984]=CEU
    [NA06985]=CEU [NA06986]=CEU [NA06989]=CEU [NA06994]=CEU
    [HG01565]=PEL [HG01566]=PEL [HG01567]=PEL [HG01568]=PEL [HG01570]=PEL
    [HG01571]=PEL [HG01572]=PEL [HG01573]=PEL [HG01574]=PEL
    [NA19648]=MXL [NA19649]=MXL [NA19651]=MXL [NA19652]=MXL [NA19654]=MXL
    [NA19655]=MXL [NA19657]=MXL [NA19658]=MXL [NA19659]=MXL
    [HG00731]=PUR [HG00732]=PUR [HG00733]=PUR [HG00734]=PUR [HG00736]=PUR
    [HG00737]=PUR [HG00738]=PUR [HG00739]=PUR [HG00740]=PUR
    [NA18939]=JPT [NA18940]=JPT [NA18941]=JPT [NA18942]=JPT [NA18943]=JPT
    [NA18944]=JPT [NA18945]=JPT [NA18946]=JPT [NA18947]=JPT
)

#-----------------------------------------------------------------------------
# fetch_30x_accessions: download ENA run report, cross-ref with GT, write TSV
#-----------------------------------------------------------------------------
fetch_30x_accessions() {
    local acc_file="${WORK_DIR}/conf/1kgp_30x_accessions.tsv"
    local ena_tsv="${WORK_DIR}/conf/ena_run_report.tsv"
    local gt_file="${GT_FILE}"

    echo "[INFO] Fetching ENA run report for PRJEB31736 (1KGP 30x NYGC) ..."
    # Use submitted_ftp: the CRAM filename encodes the 1KGP sample ID (e.g. NA19238.final.cram)
    curl -sL \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB31736&result=read_run&fields=run_accession,submitted_ftp&format=tsv" \
        > "$ena_tsv"

    local nrows
    nrows=$(tail -n +2 "$ena_tsv" | wc -l || echo 0)
    echo "[INFO] Downloaded $nrows run entries from ENA"

    if [[ ! -f "$gt_file" ]]; then
        echo "[WARN] GT file not found at $gt_file — will mark all as has_gt=unknown"
    fi

    python3 - << PYEOF
import csv, sys, os, re

gt_file = "${gt_file}"
ena_tsv = "${ena_tsv}"
out_tsv = "${acc_file}"

# Load GT samples
gt = set()
if os.path.exists(gt_file):
    with open(gt_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # skip header
        for row in reader:
            if row:
                gt.add(row[0].strip())
print(f"[INFO] GT samples loaded: {len(gt)}")

# Parse ENA report: extract sample ID from CRAM filename in submitted_ftp
# e.g. "ftp.sra.ebi.ac.uk/.../NA19238.final.cram;..." → "NA19238"
written = 0
with open(ena_tsv) as f, open(out_tsv, 'w') as out:
    out.write("sample_id\terr_accession\thas_gt\n")
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ftp = row.get('submitted_ftp', '').strip()
        err = row.get('run_accession', '').strip()
        if not ftp or not err:
            continue
        # Extract sample ID from CRAM filename: pick .final.cram entry
        m = re.search(r'/([A-Z0-9]+)\.final\.cram(?:;|$)', ftp)
        if not m:
            continue
        s = m.group(1)
        has_gt = 'yes' if s in gt else 'no'
        out.write(f"{s}\t{err}\t{has_gt}\n")
        written += 1

print(f"[INFO] Wrote {written} accessions to ${acc_file}")
gt_count = sum(1 for line in open(out_tsv).readlines()[1:] if line.split('\t')[2].strip() == 'yes')
print(f"[INFO] Samples with ground truth: {gt_count}")
PYEOF

    echo "[INFO] Accession file written: $acc_file"
}

#-----------------------------------------------------------------------------
# load_accessions: populate CRAM_ERR from downloaded TSV (keeps 7-sample fallback)
#-----------------------------------------------------------------------------
load_accessions() {
    local acc_file="${WORK_DIR}/conf/1kgp_30x_accessions.tsv"
    if [[ ! -f "$acc_file" ]]; then
        echo "[INFO] No accession file found at $acc_file — using hardcoded 7-sample fallback"
        echo "[INFO] Run with --fetch-accessions to download full accession table"
        return
    fi
    local count=0
    while IFS=$'\t' read -r sample err has_gt; do
        [[ "$has_gt" == "yes" ]] && CRAM_ERR["$sample"]="$err" && count=$((count+1))
    done < <(tail -n +2 "$acc_file")
    echo "[INFO] Loaded ${count} samples with 30x accessions + ground truth"
}

#-----------------------------------------------------------------------------
# needs_typing: returns 0 (true) if sample needs at least one tool result
#-----------------------------------------------------------------------------
needs_typing() {
    local sample="$1"
    local tool
    for tool in $(echo "$TOOLS" | tr ',' ' '); do
        local dst="${RESULTS_DIR}/by_tool/${tool}/${sample}_${tool}.txt"
        [[ ! -f "$dst" ]] && return 0
    done
    return 1
}

#-----------------------------------------------------------------------------
# all_tools_at_target: returns 0 (true) if core tools have >= TARGET_N results
#-----------------------------------------------------------------------------
all_tools_at_target() {
    for tool in hlahd spechla optitype; do
        local n
        n=$(ls "${RESULTS_DIR}/by_tool/${tool}/"*_${tool}.txt 2>/dev/null | wc -l || echo 0)
        [[ $n -lt $TARGET_N ]] && return 1
    done
    return 0
}

#-----------------------------------------------------------------------------
# prepare_batches: pre-generate numbered batch CSV files (dry run, no downloads)
#-----------------------------------------------------------------------------
prepare_batches() {
    local batches_dir="${WORK_DIR}/batches"
    mkdir -p "$batches_dir"

    # Collect all samples in priority order that still need typing
    local -a pending=()
    for S in "${ALL_SAMPLES_ORDERED[@]}"; do
        [[ -n "${CRAM_ERR[$S]:-}" ]] || continue
        needs_typing "$S" || continue
        pending+=("$S")
    done

    if [[ ${#pending[@]} -eq 0 ]]; then
        echo "[INFO] No samples pending — all tools at target or no accessions loaded."
        echo "       Run --fetch-accessions first, then --prepare-batches."
        return
    fi

    echo "=== Batch plan: ${#pending[@]} untyped samples → batches of ${BATCH_SIZE} ==="
    echo ""
    printf "%-12s  %-6s  %-25s  %s\n" "Batch file" "N" "Populations" "Samples"
    printf -- "%-12s  %-6s  %-25s  %s\n" "----------" "--" "-----------" "-------"

    local batch_num=1
    local -a batch_samples=()
    for S in "${pending[@]}"; do
        batch_samples+=("$S")
        if [[ ${#batch_samples[@]} -ge $BATCH_SIZE ]]; then
            _write_batch "$batch_num" "${batch_samples[@]}"
            batch_num=$((batch_num + 1))
            batch_samples=()
        fi
    done
    # Write remaining partial batch
    if [[ ${#batch_samples[@]} -gt 0 ]]; then
        _write_batch "$batch_num" "${batch_samples[@]}"
    fi

    echo ""
    echo "[INFO] Batch files written to: $batches_dir"
    echo ""
    echo "[INFO] To run a batch:"
    echo "  bash scripts/run_calibration_local.sh --source $SOURCE \\"
    echo "      --batch-file ${batches_dir}/batch_001.csv \\"
    echo "      --tools \"$TOOLS\" --work-dir $WORK_DIR"
}

_write_batch() {
    local batch_num="$1"; shift
    local samples=("$@")
    local batch_csv="${WORK_DIR}/batches/batch_$(printf '%03d' "$batch_num").csv"
    echo "sample_id,fastq_1,fastq_2" > "$batch_csv"
    local pops=""
    for ss in "${samples[@]}"; do
        echo "${ss},${HLA_BAM_DIR}/${ss}_R1.fastq.gz,${HLA_BAM_DIR}/${ss}_R2.fastq.gz" >> "$batch_csv"
        pops="${pops}${SAMPLE_POP[$ss]:-UNK} "
    done
    printf "batch_%03d     %-6d  %-25s  %s\n" \
        "$batch_num" "${#samples[@]}" "$pops" "${samples[*]}"
}

#-----------------------------------------------------------------------------
# --fetch-accessions mode: download and exit
#-----------------------------------------------------------------------------
if [[ "$FETCH_ACCESSIONS" == "true" ]]; then
    if [[ ! -f "$GT_FILE" ]]; then
        echo "[INFO] GT file not found — downloading first ..."
        python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" download-gt \
            --output "$GT_FILE"
    fi
    fetch_30x_accessions
    exit 0
fi

#-----------------------------------------------------------------------------
# Load dynamic accessions (30x mode only)
#-----------------------------------------------------------------------------
if [[ "$SOURCE" == "30x" ]]; then
    load_accessions
fi

#-----------------------------------------------------------------------------
# --status mode: show per-tool progress and exit
#-----------------------------------------------------------------------------
if [[ "$STATUS_ONLY" == "true" ]]; then
    echo "=== HLA Calibration Status ==="
    echo ""
    if [[ -f "$GT_FILE" ]]; then
        GT_COUNT=$(tail -n +2 "$GT_FILE" | wc -l || echo '?')
        echo "Ground truth samples: $GT_COUNT"
    else
        echo "Ground truth: NOT FOUND ($GT_FILE)"
    fi
    echo ""
    echo "FASTQs extracted: $(ls "${HLA_BAM_DIR}/"*_R1.fastq.gz 2>/dev/null | wc -l)"
    echo ""
    echo "Typing results by tool (target: ${TARGET_N}):"
    printf "  %-12s %5s / %s\n" "Tool" "Done" "$TARGET_N"
    printf "  %-12s %5s\n" "----" "----"
    for tool in hlahd spechla optitype arcashla xhla hlala; do
        N=$(ls "${RESULTS_DIR}/by_tool/${tool}/"*_${tool}.txt 2>/dev/null | wc -l || echo 0)
        printf "  %-12s %5d / %d\n" "$tool" "$N" "$TARGET_N"
    done
    echo ""
    ACC_FILE="${WORK_DIR}/conf/1kgp_30x_accessions.tsv"
    if [[ -f "$ACC_FILE" ]]; then
        ACC_COUNT=$(tail -n +2 "$ACC_FILE" | grep -c $'\tyes$' || echo 0)
        echo "ENA accessions with GT: $ACC_COUNT"
        echo "  (from ${ACC_FILE})"
    else
        echo "ENA accession file: NOT FOUND (run --fetch-accessions)"
    fi
    exit 0
fi

#-----------------------------------------------------------------------------
# --prepare-batches mode: pre-generate batch CSVs and exit
#-----------------------------------------------------------------------------
if [[ "$PREPARE_BATCHES" == "true" ]]; then
    prepare_batches
    exit 0
fi

#-----------------------------------------------------------------------------
# Validate run.config (required for typing phases)
#-----------------------------------------------------------------------------
if [[ ! -f "$RUN_CONFIG" ]]; then
    echo "ERROR: run.config not found: $RUN_CONFIG"
    echo "Set with --run-config or run from the pipeline directory"
    exit 1
fi

# Check required tools
for cmd in samtools nextflow singularity python3; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found in PATH"
        exit 1
    fi
done

#-----------------------------------------------------------------------------
# Build effective sample list
#   --batch-file FILE  → read exactly those sample IDs from the CSV
#   otherwise          → pick next N from ordered list (filtered by accession)
#-----------------------------------------------------------------------------
if [[ -n "$BATCH_FILE" ]]; then
    if [[ ! -f "$BATCH_FILE" ]]; then
        echo "ERROR: batch file not found: $BATCH_FILE"
        exit 1
    fi
    mapfile -t SAMPLES < <(tail -n +2 "$BATCH_FILE" | cut -d',' -f1)
    echo "[INFO] Loaded ${#SAMPLES[@]} samples from batch file: $BATCH_FILE"
elif [[ -n "$USER_GT_FILE" ]]; then
    # Build sample list from the provided GT CSV (first column, skip header)
    # Then filter to those with an ENA accession in CRAM_ERR map
    mapfile -t _GT_SAMPLES < <(python3 -c "
import csv
with open('${GT_FILE}') as f:
    for i, row in enumerate(csv.reader(f)):
        if i == 0: continue
        s = row[0].strip()
        if s: print(s)
")
    SAMPLES=()
    for S in "${_GT_SAMPLES[@]}"; do
        [[ -n "${CRAM_ERR[$S]:-}" ]] && SAMPLES+=("$S")
    done
    echo "[INFO] GT CSV: ${#_GT_SAMPLES[@]} samples → ${#SAMPLES[@]} have ENA accessions"
elif [[ "$SOURCE" == "30x" ]]; then
    SAMPLES=()
    for S in "${ALL_SAMPLES_ORDERED[@]}"; do
        [[ -n "${CRAM_ERR[$S]:-}" ]] && SAMPLES+=("$S")
    done
else
    SAMPLES=("${ALL_SAMPLES_ORDERED[@]:0:$N_SAMPLES}")
fi

echo "=============================================="
echo " HLA Calibration — Local"
echo "=============================================="
echo " Source:           $SOURCE"
echo " Work dir:         $WORK_DIR"
echo " Results dir:      $RESULTS_DIR"
echo " Run config:       $RUN_CONFIG"
echo " Tools:            $TOOLS"
echo " Genes:            $GENES"
[[ -n "$USER_GT_FILE" ]] && echo " GT file:          $USER_GT_FILE"
echo " Target per tool:  $TARGET_N"
echo " Parallel extract: $PARALLEL_EXTRACT"
if [[ -n "$BATCH_FILE" ]]; then
echo " Batch file:       $BATCH_FILE"
echo " Samples in batch: ${#SAMPLES[@]}"
else
echo " Candidate pool:   ${#SAMPLES[@]} samples"
fi
echo "=============================================="
echo ""

# Early exit if already at target
if all_tools_at_target; then
    echo "[INFO] Target of $TARGET_N samples already reached for core tools."
    echo "       Run with --status to see details."
    exit 0
fi

#-----------------------------------------------------------------------------
# Phase 1: Extract HLA region reads (parallel, skip if already done)
#-----------------------------------------------------------------------------
if [[ "$SKIP_EXTRACT" != "true" ]]; then
    if [[ "$SOURCE" == "30x" ]]; then
        echo "=== Phase 1: Extracting HLA region reads from 1KGP 30x CRAMs (GRCh38) ==="
        export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"
        export REF_CACHE="/tmp/ref_cache/%2s/%2s/%s"
        mkdir -p /tmp/ref_cache
    else
        echo "=== Phase 1: Extracting HLA region reads from 1KGP Phase 3 BAMs (hg19) ==="
    fi
    echo ""

    _extract_bam() {
        local SAMPLE="$1"
        local OUT_BAM="${HLA_BAM_DIR}/${SAMPLE}.bam"

        if [[ -f "${HLA_BAM_DIR}/${SAMPLE}_R1.fastq.gz" ]] && \
           [[ -f "${HLA_BAM_DIR}/${SAMPLE}_R2.fastq.gz" ]]; then
            echo "[SKIP] $SAMPLE — FASTQ already exists"
            return 0
        fi

        if [[ -f "$OUT_BAM" ]] && [[ -f "${OUT_BAM}.bai" ]]; then
            echo "[SKIP] $SAMPLE — HLA BAM already exists ($(du -sh "$OUT_BAM" | cut -f1))"
        else
            local POP="${SAMPLE_POP[$SAMPLE]:-UNK}"
            echo "[EXTRACT] $SAMPLE ($POP) ..."

            local REMOTE_URL=""
            if [[ "$SOURCE" == "30x" ]]; then
                local ERR="${CRAM_ERR[$SAMPLE]:-}"
                if [[ -z "$ERR" ]]; then
                    echo "[WARN] $SAMPLE: no 30x CRAM ERR accession — skipping"
                    return 1
                fi
                local PREFIX3="${ERR:0:6}"
                REMOTE_URL="ftp://ftp.sra.ebi.ac.uk/vol1/run/${PREFIX3}/${ERR}/${SAMPLE}.final.cram"
            else
                local POP_LC="${SAMPLE_POP[$SAMPLE]:-UNK}"
                REMOTE_URL=$(curl -s --list-only \
                    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/${SAMPLE}/alignment/" \
                    2>/dev/null \
                    | grep "${SAMPLE}.mapped.ILLUMINA.bwa.${POP_LC}.low_coverage.*\.bam$" \
                    | grep -v '\.bai$' \
                    | head -1 \
                    | awk -v base="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/${SAMPLE}/alignment" \
                      '{print base "/" $NF}' \
                    || true)
                if [[ -z "$REMOTE_URL" ]]; then
                    echo "[WARN] $SAMPLE: could not find phase3 BAM URL — skipping"
                    return 1
                fi
            fi

            echo "       URL: $REMOTE_URL"

            if samtools view -b \
                -o "${OUT_BAM}.tmp" \
                "${REMOTE_URL}" \
                "$HLA_REGION" 2>>"${LOGS_DIR}/extract_${SAMPLE}.log"; then

                mv "${OUT_BAM}.tmp" "$OUT_BAM"
                samtools sort -o "${OUT_BAM%.bam}.sorted.bam" "$OUT_BAM" \
                    2>>"${LOGS_DIR}/extract_${SAMPLE}.log"
                mv "${OUT_BAM%.bam}.sorted.bam" "$OUT_BAM"
                samtools index "$OUT_BAM" 2>>"${LOGS_DIR}/extract_${SAMPLE}.log"

                READS=$(samtools view -c "$OUT_BAM" 2>/dev/null || echo "?")
                SIZE=$(du -sh "$OUT_BAM" | cut -f1)
                echo "       [OK] $SAMPLE: $SIZE, $READS reads"
            else
                echo "[WARN] $SAMPLE: samtools extraction failed (see ${LOGS_DIR}/extract_${SAMPLE}.log)"
                rm -f "${OUT_BAM}.tmp"
                return 1
            fi
        fi

        # Convert to FASTQ (avoids chr naming issues inside Singularity containers)
        echo "       [FASTQ] $SAMPLE — converting to paired FASTQ..."
        samtools sort -n -@ 4 "$OUT_BAM" -o "${OUT_BAM%.bam}_namesorted.bam" \
            2>>"${LOGS_DIR}/extract_${SAMPLE}.log"
        samtools fastq -@ 4 \
            -1 "${HLA_BAM_DIR}/${SAMPLE}_R1.fastq.gz" \
            -2 "${HLA_BAM_DIR}/${SAMPLE}_R2.fastq.gz" \
            -0 /dev/null -s /dev/null \
            "${OUT_BAM%.bam}_namesorted.bam" \
            2>>"${LOGS_DIR}/extract_${SAMPLE}.log"
        rm -f "${OUT_BAM%.bam}_namesorted.bam"
        echo "       [OK] $SAMPLE: FASTQ conversion done"
    }

    # Parallel extraction: up to PARALLEL_EXTRACT concurrent jobs
    EXTRACT_PIDS=()
    EXTRACT_COUNT=0
    for SAMPLE in "${SAMPLES[@]}"; do
        # Only extract up to N_SAMPLES new samples that still need typing
        if needs_typing "$SAMPLE" || \
           [[ ! -f "${HLA_BAM_DIR}/${SAMPLE}_R1.fastq.gz" ]]; then
            _extract_bam "$SAMPLE" &
            EXTRACT_PIDS+=($!)
            EXTRACT_COUNT=$((EXTRACT_COUNT + 1))

            # Throttle parallel jobs
            if [[ ${#EXTRACT_PIDS[@]} -ge $PARALLEL_EXTRACT ]]; then
                wait "${EXTRACT_PIDS[0]}" || true
                EXTRACT_PIDS=("${EXTRACT_PIDS[@]:1}")
            fi

            [[ $EXTRACT_COUNT -ge $N_SAMPLES ]] && break
        fi
    done
    # Drain remaining background jobs
    for pid in "${EXTRACT_PIDS[@]}"; do
        wait "$pid" || true
    done
    echo ""
fi

#-----------------------------------------------------------------------------
# Phase 2: Build batch samplesheet of all FASTQs that need typing
#-----------------------------------------------------------------------------
if [[ "$SKIP_TYPING" != "true" ]]; then
    echo "=== Phase 2: HLA typing with tools: $TOOLS ==="
    echo ""

    BATCH_CSV="${WORK_DIR}/conf/batch_samplesheet_$(date +%Y%m%d_%H%M%S).csv"
    echo "sample_id,fastq_1,fastq_2" > "$BATCH_CSV"
    BATCH_COUNT=0

    for SAMPLE in "${SAMPLES[@]}"; do
        R1="${HLA_BAM_DIR}/${SAMPLE}_R1.fastq.gz"
        R2="${HLA_BAM_DIR}/${SAMPLE}_R2.fastq.gz"
        if [[ -f "$R1" ]] && [[ -f "$R2" ]] && needs_typing "$SAMPLE"; then
            echo "${SAMPLE},${R1},${R2}" >> "$BATCH_CSV"
            BATCH_COUNT=$((BATCH_COUNT + 1))
            echo "  [QUEUED] $SAMPLE"
        fi
        [[ $BATCH_COUNT -ge $N_SAMPLES ]] && break
    done

    if [[ $BATCH_COUNT -eq 0 ]]; then
        echo "[INFO] No samples need typing — all results present or no FASTQs available."
    else
        echo ""
        echo "[BATCH] Running Nextflow for $BATCH_COUNT samples ..."
        WORK_NF="${WORK_DIR}/work/batch_$(date +%Y%m%d_%H%M%S)"
        mkdir -p "$WORK_NF"
        LOG_BATCH="${LOGS_DIR}/batch_$(date +%Y%m%d_%H%M%S).log"

        nextflow run "${PIPELINE_DIR}/main.nf" \
            -c "$RUN_CONFIG" \
            --input_samplesheet "$BATCH_CSV" \
            --outdir "$RESULTS_DIR" \
            --tools "$TOOLS" \
            --seq_type dna \
            --resolution 2-field \
            --max_cpus 4 \
            --max_memory "12GB" \
            -profile singularity \
            -work-dir "$WORK_NF" \
            -resume \
            2>&1 | tee "$LOG_BATCH"

        echo "[DONE] Nextflow batch complete."
    fi

    echo ""
fi

#-----------------------------------------------------------------------------
# Phase 3: Collect results into by_tool layout
#-----------------------------------------------------------------------------
echo "=== Phase 3: Collecting results ==="
for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
    mkdir -p "${RESULTS_DIR}/by_tool/${TOOL}"
done

for SAMPLE in "${SAMPLES[@]}"; do
    for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
        SRC="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
        DST="${RESULTS_DIR}/by_tool/${TOOL}/${SAMPLE}_${TOOL}.txt"
        if [[ -f "$SRC" ]] && [[ ! -e "$DST" ]]; then
            ln -sf "$SRC" "$DST"
        fi
    done
done

echo "Results by tool:"
for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
    N=$(ls "${RESULTS_DIR}/by_tool/${TOOL}/"*_${TOOL}.txt 2>/dev/null | wc -l || echo 0)
    printf "  %-12s %3d / %d\n" "$TOOL" "$N" "$TARGET_N"
done
echo ""

#-----------------------------------------------------------------------------
# Phase 4: Metrics from Nextflow trace
#-----------------------------------------------------------------------------
echo "=== Phase 4: Collecting pipeline metrics ==="
TRACE_FILE=$(ls -t "${RESULTS_DIR}/pipeline_info/trace_"*.txt 2>/dev/null | head -1 || true)
if [[ -n "$TRACE_FILE" ]]; then
    METRICS_DIR="${WORK_DIR}/metrics"
    mkdir -p "$METRICS_DIR"
    echo "[INFO] Parsing trace: $TRACE_FILE"
    python3 "${PIPELINE_DIR}/bin/hla_pipeline_metrics.py" \
        --trace "$TRACE_FILE" \
        --outdir "$METRICS_DIR" \
        2>/dev/null || true
    echo "[INFO] Metrics written to $METRICS_DIR"
else
    echo "[INFO] No Nextflow trace file found — skipping metrics"
fi
echo ""

#-----------------------------------------------------------------------------
# Phase 5: Download ground truth + population file + run calibration
#-----------------------------------------------------------------------------
echo "=== Phase 5: Calibration ==="

if [[ ! -f "$GT_FILE" ]]; then
    echo "Downloading 1KGP ground truth ..."
    python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" download-gt \
        --output "$GT_FILE"
else
    echo "Ground truth exists: $GT_FILE"
fi

# Generate / refresh population file (1KGP panel → sample→superpopulation TSV)
POP_FILE="${WORK_DIR}/conf/1kgp_populations.tsv"
if [[ ! -f "$POP_FILE" ]]; then
    echo ""
    echo "Downloading 1KGP population panel ..."
    python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" generate-population-file \
        --output "$POP_FILE" \
        && echo "[OK] Population file: $POP_FILE" \
        || echo "[WARN] Could not download population file — continuing without stratification"
else
    echo "Population file exists: $POP_FILE"
fi

echo ""
echo "Running calibration (genes: $GENES) ..."

POP_FILE_ARG=""
[[ -f "$POP_FILE" ]] && POP_FILE_ARG="--population-file $POP_FILE"

# shellcheck disable=SC2086
python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" calibrate \
    --ground-truth "$GT_FILE" \
    --results-dir "${RESULTS_DIR}/by_tool" \
    --data-type wgs \
    --genes "$GENES" \
    --resolution 2-field \
    --output-weights "$WEIGHTS_OUT" \
    --output-table "$TABLE_OUT" \
    --target-n "$TARGET_N" \
    $POP_FILE_ARG

echo ""
echo "[INFO] Calibration output:"
echo "         Weights: $WEIGHTS_OUT"
echo "         Table:   $TABLE_OUT"
if [[ -f "$POP_FILE" ]]; then
    POP_TABLE_DIR="$(dirname "$TABLE_OUT")"
    echo "         Per-pop:  ${POP_TABLE_DIR}/tool_accuracy_wgs_v2_*.tsv"
fi

#-----------------------------------------------------------------------------
# Phase 6: Compare voting strategies (Wilcoxon signed-rank test)
#-----------------------------------------------------------------------------
echo ""
echo "=== Phase 6: Voting strategy comparison ==="

COMPARE_OUT="${WORK_DIR}/conf/strategy_comparison_$(date +%Y%m%d).tsv"

if [[ -f "$WEIGHTS_OUT" ]]; then
    echo "Comparing equal vs calibrated voting (Wilcoxon signed-rank test) ..."
    python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" compare-strategies \
        --ground-truth "$GT_FILE" \
        --results-dir "${RESULTS_DIR}/by_tool" \
        --genes "$GENES" \
        --resolution 2-field \
        --weights-file "$WEIGHTS_OUT" \
        --ref-mode equal \
        --test-mode calibrated \
        --output "$COMPARE_OUT" \
        && echo "[OK] Strategy comparison: $COMPARE_OUT" \
        || echo "[WARN] Strategy comparison failed — check logs"

    echo ""
    echo "Comparing equal vs read_confidence voting ..."
    COMPARE_RC_OUT="${WORK_DIR}/conf/strategy_comparison_rc_$(date +%Y%m%d).tsv"
    python3 "${PIPELINE_DIR}/bin/calibrate_tool_weights.py" compare-strategies \
        --ground-truth "$GT_FILE" \
        --results-dir "${RESULTS_DIR}/by_tool" \
        --genes "$GENES" \
        --resolution 2-field \
        --ref-mode equal \
        --test-mode read_confidence \
        --output "$COMPARE_RC_OUT" \
        && echo "[OK] RC comparison: $COMPARE_RC_OUT" \
        || echo "[WARN] RC comparison failed"
else
    echo "[SKIP] Calibrated weights not found — skipping strategy comparison."
    echo "       (Run this step again after calibration produces $WEIGHTS_OUT)"
fi

echo ""
echo "=============================================="
echo " Calibration + comparison complete"
echo "=============================================="
echo " Weights:         $WEIGHTS_OUT"
echo " Accuracy table:  $TABLE_OUT"
if [[ -f "$COMPARE_OUT" ]]; then
echo " Strategy comp:   $COMPARE_OUT"
fi
echo ""
echo " Next steps:"
echo "   # Check progress"
echo "   bash scripts/run_calibration_local.sh --source $SOURCE --status \\"
echo "       --work-dir $WORK_DIR"
echo ""
echo "   # Run next batch (add more samples)"
echo "   bash scripts/run_calibration_local.sh --source $SOURCE \\"
echo "       --n-samples $N_SAMPLES --tools \"$TOOLS\" \\"
echo "       --work-dir $WORK_DIR --target $TARGET_N"
echo "=============================================="
