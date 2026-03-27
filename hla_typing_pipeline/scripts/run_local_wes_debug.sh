#!/usr/bin/env bash
# Local reduced WES debug runner using the shared BAM preprocessing path.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
TOOLS="hlahd,spechla,arcashla,optitype"
STAMP="$(date +%Y%m%d_%H%M%S)"

INPUT_BAM="${1:-}"
RUN_NAME="${2:-local_wes_debug_${STAMP}}"
OUTDIR="${PIPELINE_DIR}/${RUN_NAME}"
WORKDIR="${PIPELINE_DIR}/work_${RUN_NAME}"
USER_CONFIG="${PIPELINE_DIR}/conf/user.config"
DEBUG_CONFIG="${PIPELINE_DIR}/conf/local_debug_wes.config"
RUN_CONFIG="${PIPELINE_DIR}/run.config"
INPUT_BAI=""

require_file() {
    local path="$1"
    if [[ ! -f "$path" ]]; then
        echo "ERROR: required file not found: $path" >&2
        exit 1
    fi
}

require_dir() {
    local path="$1"
    if [[ ! -d "$path" ]]; then
        echo "ERROR: required directory not found: $path" >&2
        exit 1
    fi
}

resolve_bam_index() {
    local bam="$1"
    if [[ -f "${bam}.bai" ]]; then
        INPUT_BAI="${bam}.bai"
        return 0
    fi
    if [[ "$bam" == *.bam ]]; then
        local alt_bai="${bam%.bam}.bai"
        if [[ -f "$alt_bai" ]]; then
            INPUT_BAI="$alt_bai"
            return 0
        fi
    fi
    echo "ERROR: BAM index not found for $bam (checked ${bam}.bai and ${bam%.bam}.bai)" >&2
    exit 1
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || -z "${1:-}" ]]; then
    cat <<EOF
Usage: $0 <input_bam> [run_name]

Runs one local reduced WES debug pass with:
  ${TOOLS}

Output:
  ${PIPELINE_DIR}/<run_name>
  ${PIPELINE_DIR}/work_<run_name>
EOF
    [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && exit 0 || exit 1
fi

require_file "$INPUT_BAM"
resolve_bam_index "$INPUT_BAM"
require_file "$USER_CONFIG"
require_file "$DEBUG_CONFIG"
require_file "$RUN_CONFIG"

command -v nextflow >/dev/null || { echo "ERROR: nextflow not found in PATH" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity not found in PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

require_dir "/home/umut/SpecHLA"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/hlahd.sif"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/arcashla.sif"
require_file "/home/umut/projects/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_references/containers/optitype.sif"

mkdir -p "$OUTDIR"

echo "==================================================================="
echo " Local WES Debug"
echo "==================================================================="
echo "Input BAM : $INPUT_BAM"
echo "BAM index : $INPUT_BAI"
echo "Run name  : $RUN_NAME"
echo "Output    : $OUTDIR"
echo "Work dir  : $WORKDIR"
echo "Tools     : $TOOLS"
echo "Configs   : conf/user.config + run.config + conf/local_debug_wes.config"
echo "==================================================================="

cd "$PIPELINE_DIR"

set +e
nextflow run main.nf \
    -c "$USER_CONFIG" \
    -c "$RUN_CONFIG" \
    -c "$DEBUG_CONFIG" \
    -profile singularity \
    -work-dir "$WORKDIR" \
    --input_bam "$INPUT_BAM" \
    --outdir "$OUTDIR" \
    --tools "$TOOLS"
nf_exit=$?
set -e

echo
echo "==================================================================="
echo " Run Summary"
echo "==================================================================="
echo "Nextflow exit code: $nf_exit"
for tool in hlahd spechla arcashla optitype; do
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
    awk -F'\t' 'NR==1 || $4 ~ /EXTRACT_HLA_READS|SPECHLA|HLAHD|ARCASHLA|OPTITYPE_FASTQ/' "$TRACE_FILE"
fi

exit "$nf_exit"
