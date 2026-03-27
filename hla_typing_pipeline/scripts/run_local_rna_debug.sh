#!/usr/bin/env bash
# Local reduced RNA debug runner for the RNA-compatible tool set.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
TOOLS="arcashla,optitype,t1k,seq2hla"
STAMP="$(date +%Y%m%d_%H%M%S)"

FASTQ1="${1:-}"
FASTQ2="${2:-}"
RUN_NAME="${3:-local_rna_debug_${STAMP}}"
OUTDIR="${PIPELINE_DIR}/${RUN_NAME}"
WORKDIR="${PIPELINE_DIR}/work_${RUN_NAME}"
USER_CONFIG="${PIPELINE_DIR}/conf/user.config"
DEBUG_CONFIG="${PIPELINE_DIR}/conf/local_debug_rna.config"
RUN_CONFIG="${PIPELINE_DIR}/run.config"

require_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: required file not found: $path" >&2
        exit 1
    fi
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || -z "${1:-}" || -z "${2:-}" ]]; then
    cat <<EOF
Usage: $0 <fastq_r1> <fastq_r2> [run_name]

Runs one local reduced RNA debug pass with:
  ${TOOLS}

Output:
  ${PIPELINE_DIR}/<run_name>
  ${PIPELINE_DIR}/work_<run_name>
EOF
    [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && exit 0 || exit 1
fi

require_file "$FASTQ1"
require_file "$FASTQ2"
require_file "$USER_CONFIG"
require_file "$DEBUG_CONFIG"
require_file "$RUN_CONFIG"

command -v nextflow >/dev/null || { echo "ERROR: nextflow not found in PATH" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity not found in PATH" >&2; exit 1; }

require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/arcashla.sif"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/optitype.sif"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/t1k.sif"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/seq2hla.sif"

mkdir -p "$OUTDIR"

echo "==================================================================="
echo " Local RNA Debug"
echo "==================================================================="
echo "FASTQ R1  : $FASTQ1"
echo "FASTQ R2  : $FASTQ2"
echo "Run name  : $RUN_NAME"
echo "Output    : $OUTDIR"
echo "Work dir  : $WORKDIR"
echo "Tools     : $TOOLS"
echo "Configs   : conf/user.config + run.config + conf/local_debug_rna.config"
echo "==================================================================="

cd "$PIPELINE_DIR"

set +e
nextflow run main.nf \
    -c "$USER_CONFIG" \
    -c "$RUN_CONFIG" \
    -c "$DEBUG_CONFIG" \
    -profile singularity \
    -work-dir "$WORKDIR" \
    --input_fastq_1 "$FASTQ1" \
    --input_fastq_2 "$FASTQ2" \
    --outdir "$OUTDIR" \
    --tools "$TOOLS"
nf_exit=$?
set -e

echo
echo "==================================================================="
echo " Run Summary"
echo "==================================================================="
echo "Nextflow exit code: $nf_exit"
for tool in arcashla optitype t1k seq2hla; do
    result_file="$(find "$OUTDIR" -path "*/${tool}/*_${tool}.txt" -type f | head -n 1 || true)"
    if [[ -n "$result_file" ]]; then
        echo "[OK] ${tool}: ${result_file}"
    else
        echo "[WARN] ${tool}: no published result file"
    fi
done

TRACE_FILE="$(find "$OUTDIR/pipeline_info" -maxdepth 1 -name 'trace_*.txt' -type f | sort | tail -n 1 || true)"
if [[ -n "$TRACE_FILE" ]]; then
    echo
    echo "Latest trace: $TRACE_FILE"
    awk -F'\t' 'NR==1 || $4 ~ /ARCASHLA_FASTQ|OPTITYPE_FASTQ|T1K_FASTQ|SEQ2HLA/' "$TRACE_FILE"
fi

exit "$nf_exit"
