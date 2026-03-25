#!/bin/bash
#=============================================================================
# Re-parse already-completed seq2HLA / Kourami / POLYSOLVER work directories.
#
# HLA typing already ran successfully on all 131 1KGP samples, but the parse
# step failed because PIPELINE_BIN had a doubled path. This script re-runs
# only the parse step without re-running the expensive HLA typing.
#
# Run on the Puhti login node after `git pull` (once PIPELINE_BIN is fixed).
# Usage: bash scripts/reparse_completed_results.sh
#=============================================================================
set -euo pipefail

PROJECT_ID="project_2008084"
BASE="/scratch/${PROJECT_ID}/ozcanumu/hla_calibration"
PIPELINE_BIN="/scratch/${PROJECT_ID}/ozcanumu/new_pipeline_2/hla_typing_pipeline/bin"
RESULTS_DIR="${BASE}/1kgp_typing_results"

module load python-data

echo "[INFO] PIPELINE_BIN = ${PIPELINE_BIN}"
echo "[INFO] RESULTS_DIR  = ${RESULTS_DIR}"
echo ""

# ── seq2HLA ──────────────────────────────────────────────────────────────────
echo "=== seq2HLA ==="
WORK_BASE="${BASE}/seq2hla_work"
mkdir -p "${RESULTS_DIR}/by_tool/seq2hla"
N_OK=0; N_SKIP=0; N_FAIL=0

for WORKDIR in "${WORK_BASE}"/*/; do
    [[ -d "$WORKDIR" ]] || continue
    SAMPLE=$(basename "$WORKDIR")
    OUT_TSV="${RESULTS_DIR}/by_tool/seq2hla/${SAMPLE}_seq2hla.txt"
    [[ -f "$OUT_TSV" ]] && { (( N_OK++ )) || true; continue; }

    PREFIX="${WORKDIR}${SAMPLE}."
    if [[ ! -f "${PREFIX}-ClassI-class.HLAgenotype4digits" ]]; then
        echo "  [SKIP] ${SAMPLE} — no seq2hla output in ${WORKDIR}"
        (( N_SKIP++ )) || true
        continue
    fi

    echo "  [PARSE] ${SAMPLE}"
    if python3 "${PIPELINE_BIN}/parse_seq2hla_results.py" \
            --sample  "$SAMPLE" \
            --prefix  "$PREFIX" \
            --output  "$OUT_TSV"; then
        mkdir -p "${RESULTS_DIR}/${SAMPLE}/seq2hla"
        cp "$OUT_TSV" "${RESULTS_DIR}/${SAMPLE}/seq2hla/${SAMPLE}_seq2hla.txt"
        (( N_OK++ )) || true
    else
        echo "  [WARN] parse failed for ${SAMPLE}"
        (( N_FAIL++ )) || true
    fi
done
echo "  seq2HLA: ${N_OK} parsed, ${N_SKIP} skipped (no output), ${N_FAIL} failed"

# ── Kourami ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Kourami ==="
WORK_BASE="${BASE}/kourami_work"
mkdir -p "${RESULTS_DIR}/by_tool/kourami"
N_OK=0; N_SKIP=0; N_FAIL=0

for WORKDIR in "${WORK_BASE}"/*/; do
    [[ -d "$WORKDIR" ]] || continue
    SAMPLE=$(basename "$WORKDIR")
    OUT_TSV="${RESULTS_DIR}/by_tool/kourami/${SAMPLE}_kourami.txt"
    [[ -f "$OUT_TSV" ]] && { (( N_OK++ )) || true; continue; }

    RESULT="${WORKDIR}${SAMPLE}.kourami.result"
    if [[ ! -f "$RESULT" ]]; then
        echo "  [SKIP] ${SAMPLE} — no kourami result in ${WORKDIR}"
        (( N_SKIP++ )) || true
        continue
    fi

    echo "  [PARSE] ${SAMPLE}"
    if python3 "${PIPELINE_BIN}/parse_kourami_results.py" \
            --input   "$RESULT" \
            --sample  "$SAMPLE" \
            --output  "$OUT_TSV"; then
        mkdir -p "${RESULTS_DIR}/${SAMPLE}/kourami"
        cp "$OUT_TSV" "${RESULTS_DIR}/${SAMPLE}/kourami/${SAMPLE}_kourami.txt"
        (( N_OK++ )) || true
    else
        echo "  [WARN] parse failed for ${SAMPLE}"
        (( N_FAIL++ )) || true
    fi
done
echo "  Kourami: ${N_OK} parsed, ${N_SKIP} skipped (no output), ${N_FAIL} failed"

# ── POLYSOLVER ────────────────────────────────────────────────────────────────
echo ""
echo "=== POLYSOLVER ==="
WORK_BASE="${BASE}/polysolver_work"
mkdir -p "${RESULTS_DIR}/by_tool/polysolver"
N_OK=0; N_SKIP=0; N_FAIL=0

for WORKDIR in "${WORK_BASE}"/*/; do
    [[ -d "$WORKDIR" ]] || continue
    SAMPLE=$(basename "$WORKDIR")
    OUT_TSV="${RESULTS_DIR}/by_tool/polysolver/${SAMPLE}_polysolver.txt"
    [[ -f "$OUT_TSV" ]] && { (( N_OK++ )) || true; continue; }

    # POLYSOLVER produces winners.hla.nofreq.txt (include_freq=0)
    # or winners.hla.txt (include_freq=1); support both
    WINNERS="${WORKDIR}/polysolver_out/winners.hla.nofreq.txt"
    [[ -f "$WINNERS" ]] || WINNERS="${WORKDIR}/polysolver_out/winners.hla.txt"
    if [[ ! -f "$WINNERS" ]]; then
        echo "  [SKIP] ${SAMPLE} — no polysolver winners file in ${WORKDIR}/polysolver_out/"
        (( N_SKIP++ )) || true
        continue
    fi

    echo "  [PARSE] ${SAMPLE} (winners: $(basename "$WINNERS"))"
    if python3 "${PIPELINE_BIN}/parse_polysolver_results.py" \
            --input   "$WINNERS" \
            --sample  "$SAMPLE" \
            --output  "$OUT_TSV"; then
        mkdir -p "${RESULTS_DIR}/${SAMPLE}/polysolver"
        cp "$OUT_TSV" "${RESULTS_DIR}/${SAMPLE}/polysolver/${SAMPLE}_polysolver.txt"
        (( N_OK++ )) || true
    else
        echo "  [WARN] parse failed for ${SAMPLE}"
        (( N_FAIL++ )) || true
    fi
done
echo "  POLYSOLVER: ${N_OK} parsed, ${N_SKIP} skipped (no output), ${N_FAIL} failed"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== by_tool/ counts ==="
for TOOL in seq2hla kourami polysolver; do
    N=$(ls "${RESULTS_DIR}/by_tool/${TOOL}/"*_${TOOL}.txt 2>/dev/null | wc -l || echo 0)
    echo "  ${TOOL}: ${N} samples"
done
echo ""
echo "[DONE] Run scripts/run_calibration_puhti.sh next to update weights."
