#!/bin/bash
#=============================================================================
# WGS calibration runner for Puhti HPC
#=============================================================================
# Runs calibrate_tool_weights.py against all tools present in by_tool/.
# Auto-discovers tool directories — just ensure seq2hla/, kourami/, polysolver/
# are populated before running (via run_*_array.sh --refresh after git pull).
#
# Prerequisites:
#   - All tool array jobs completed + parsed (winners/results → {sample}_{tool}.txt)
#   - ground-truth file at ${BASE}/hla_calibration/conf/1kgp_hla_gt.tsv
#   - python-data module available
#
# Usage (from Puhti login node):
#   bash scripts/run_calibration_puhti.sh          # produces tool_weights_wgs_v3.json
#   bash scripts/run_calibration_puhti.sh v4       # produces tool_weights_wgs_v4.json
#=============================================================================
set -euo pipefail

PROJECT_ID="project_2008084"
BASE="/scratch/${PROJECT_ID}/ozcanumu"
PIPELINE="${BASE}/new_pipeline_2/hla_typing_pipeline"
RESULTS_DIR="${BASE}/hla_calibration/1kgp_typing_results/by_tool"
GT="${BASE}/hla_calibration/conf/1kgp_hla_gt.tsv"
VERSION="${1:-v3}"

echo "=== WGS Calibration — version ${VERSION} ==="
echo "Date: $(date)"
echo "Results dir: ${RESULTS_DIR}"
echo ""

#-----------------------------------------------------------------------------
# Load modules
#-----------------------------------------------------------------------------
module purge
module load python-data
# Note: python-data on Puhti includes numpy/scipy/pandas needed by calibrate script

#-----------------------------------------------------------------------------
# Sanity checks
#-----------------------------------------------------------------------------
if [[ ! -f "$GT" ]]; then
    echo "ERROR: Ground-truth file not found: ${GT}"
    echo "       Run: python3 ${PIPELINE}/bin/calibrate_tool_weights.py download-gt --output ${GT}"
    exit 1
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: Results directory not found: ${RESULTS_DIR}"
    exit 1
fi

echo "[INFO] Tools found in by_tool/:"
ls "${RESULTS_DIR}/"
echo ""

# Require at least one new tool to be present
for TOOL in seq2hla kourami polysolver; do
    if [[ ! -d "${RESULTS_DIR}/${TOOL}" ]]; then
        echo "WARNING: ${TOOL}/ not found in by_tool/ — will be skipped in calibration"
    else
        N=$(ls "${RESULTS_DIR}/${TOOL}/" 2>/dev/null | wc -l)
        echo "[INFO] ${TOOL}: ${N} result files"
    fi
done
echo ""

#-----------------------------------------------------------------------------
# Run calibration
#-----------------------------------------------------------------------------
OUT_WEIGHTS="${PIPELINE}/conf/tool_weights_wgs_${VERSION}.json"
OUT_TABLE="${PIPELINE}/conf/tool_accuracy_wgs_${VERSION}.tsv"

echo "[Step 1] Running calibrate_tool_weights.py calibrate ..."
python3 "${PIPELINE}/bin/calibrate_tool_weights.py" calibrate \
    --ground-truth   "$GT" \
    --results-dir    "$RESULTS_DIR" \
    --data-type      wgs \
    --output-weights "$OUT_WEIGHTS" \
    --output-table   "$OUT_TABLE"

echo ""
echo "[OK] Weights written to: ${OUT_WEIGHTS}"
echo "[OK] Accuracy table:     ${OUT_TABLE}"
echo ""

#-----------------------------------------------------------------------------
# Print weight summary
#-----------------------------------------------------------------------------
echo "=== Weight summary (Class I) ==="
python3 - <<'PYEOF'
import json, sys, os

path = os.environ.get('OUT_WEIGHTS', '')
if not path:
    # fallback: read from first arg or stdin
    sys.exit(0)
try:
    w = json.load(open(path))
    genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
    tools = sorted(w['genes']['A'].keys())
    header = f"{'Tool':12s}" + "".join(f"  {g:>7s}" for g in genes)
    print(header)
    print("-" * len(header))
    for t in tools:
        row = f"{t:12s}" + "".join(
            f"  {w['genes'][g].get(t, 0.0):7.3f}" for g in genes
        )
        print(row)
    print(f"\nn_samples = {w.get('n_samples', '?')}")
except Exception as e:
    print(f"(could not pretty-print: {e})")
PYEOF

OUT_WEIGHTS="$OUT_WEIGHTS" python3 - <<'PYEOF'
import json, sys, os
path = os.environ['OUT_WEIGHTS']
try:
    w = json.load(open(path))
    genes = ['A', 'B', 'C', 'DRB1', 'DQB1']
    tools = sorted(w['genes']['A'].keys())
    header = f"{'Tool':12s}" + "".join(f"  {g:>7s}" for g in genes)
    print(header)
    print("-" * len(header))
    for t in tools:
        row = f"{t:12s}" + "".join(
            f"  {w['genes'][g].get(t, 0.0):7.3f}" for g in genes
        )
        print(row)
    print(f"\nn_samples = {w.get('n_samples', '?')}")
except Exception as e:
    print(f"(could not pretty-print: {e})")
PYEOF

echo ""
echo "=== Next steps ==="
echo "1. git add conf/tool_weights_wgs_${VERSION}.json conf/tool_accuracy_wgs_${VERSION}.tsv"
echo "2. git commit -m 'feat: WGS calibration ${VERSION} with seq2hla/kourami/polysolver (N=131)'"
echo "3. git push"
echo "4. Locally: git pull && cp conf/tool_weights_wgs_${VERSION}.json conf/tool_weights_wgs.json"
