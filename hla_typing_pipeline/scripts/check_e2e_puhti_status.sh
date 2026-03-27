#!/bin/bash
# =============================================================================
# check_e2e_puhti_status.sh
# Summarize Puhti E2E test state without modifying scratch contents.
#
# Usage:
#   bash scripts/check_e2e_puhti_status.sh [--project PROJECT_ID] [--tail N]
#
# Default:
#   --project  project_2008084
#   --tail     20
# =============================================================================
set -euo pipefail

PROJECT_ID="${SLURM_JOB_ACCOUNT:-project_2008084}"
TAIL_LINES=20

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project) PROJECT_ID="$2"; shift 2 ;;
        --tail)    TAIL_LINES="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

SCRATCH_ROOT="/scratch/${PROJECT_ID}/hla_calibration"

print_file_mtime() {
    local path="$1"
    if [[ -f "$path" ]]; then
        stat -c '%y' "$path"
    else
        echo "MISSING"
    fi
}

print_log_tail() {
    local title="$1"
    local path="$2"
    echo ""
    echo "${title}: ${path}"
    if [[ -f "$path" ]]; then
        tail -n "${TAIL_LINES}" "$path"
    else
        echo "MISSING"
    fi
}

report_phase0_state() {
    local type="$1"
    local sample="$2"
    local base="$3"

    if [[ "$type" == "rna" ]]; then
        local r1="${base}/fastqs/${sample}_R1.fastq.gz"
        local r2="${base}/fastqs/${sample}_R2.fastq.gz"
        if [[ -f "$r1" && -f "$r2" ]]; then
            echo "Phase 0: FASTQs present"
        else
            echo "Phase 0: FASTQs missing"
        fi
    else
        local bam="${base}/bams/${sample}_hla.bam"
        if [[ -f "$bam" ]]; then
            echo "Phase 0: HLA BAM present"
        else
            echo "Phase 0: HLA BAM missing"
        fi
    fi
}

report_phase1_state() {
    local type="$1"
    local sample="$2"
    local base="$3"
    shift 3
    local tools=("$@")
    local found=0
    local total="${#tools[@]}"

    for tool in "${tools[@]}"; do
        local res="${base}/results/${sample}/${tool}/${sample}_${tool}.txt"
        [[ -f "$res" ]] && found=$((found+1))
    done

    echo "Phase 1: ${found}/${total} tool outputs present"
}

summarize_type() {
    local type="$1"
    local sample="$2"
    shift 2
    local tools=("$@")
    local base="${SCRATCH_ROOT}/test_${type}"
    local logs="${base}/logs"
    local report="${base}/report_${sample}.txt"

    echo "==================================================================="
    echo "E2E ${type^^} / ${sample}"
    echo "Scratch: ${base}"
    report_phase0_state "$type" "$sample" "$base"
    report_phase1_state "$type" "$sample" "$base" "${tools[@]}"
    echo "Report mtime: $(print_file_mtime "$report")"

    if command -v squeue >/dev/null 2>&1; then
        echo ""
        echo "Queue snapshot:"
        squeue -u "$USER" | awk 'NR==1 || /test_wgs|test_wes|test_rna|nf-HLAHD|nf-SPECH|nf-ARCAS|nf-OPTIT|nf-POLYS|nf-T1K|nf-SEQ2H/'
    fi

    print_log_tail "Phase 0 stderr" "${logs}/test_phase0.err"
    print_log_tail "Phase 1 stderr" "${logs}/test_typing.err"
    print_log_tail "Phase 2 stderr" "${logs}/test_report.err"

    if [[ -f "$report" ]]; then
        print_log_tail "Report tail" "$report"
    fi
}

summarize_type "wgs" "NA19238" hlahd spechla arcashla optitype polysolver
summarize_type "wes" "NA18501" hlahd spechla arcashla optitype polysolver
summarize_type "rna" "NA18502" hlahd spechla arcashla optitype t1k seq2hla
