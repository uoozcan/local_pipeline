#!/bin/bash
#=============================================================================
# HLA Typing Pipeline — In-Place Update Script for CSC Puhti
#=============================================================================
# Updates an EXISTING Puhti installation from v1.3.x to v1.4.0.
# Does NOT touch your data (scratch directory is read-only to this script).
#
# Usage
# -----
#   # Auto-detect project from $CSC_PROJECT or current directory
#   bash scripts/update_puhti.sh
#
#   # Specify project explicitly
#   bash scripts/update_puhti.sh project_2008084
#
#   # Preview what would happen without making changes
#   bash scripts/update_puhti.sh --dry-run
#
#   # Skip container rebuild/pull (just update code)
#   bash scripts/update_puhti.sh --skip-containers
#
# What this script does
# ---------------------
#   1. git pull origin main (in the install directory)
#   2. Rebuild hla_postprocess.sif  (adds plotly, kaleido for pipeline metrics)
#   3. Pull t1k.sif                  (new in v1.4.0; needed for calibration)
#   4. Pull hifihla.sif              (optional; only if --include-longreads)
#   5. Validate new v1.4.0 files are present
#   6. Check conf/user.config for missing v1.4.0 params
#
# Requirements
# ------------
#   - Pipeline was previously installed with scripts/install_puhti.sh
#   - Run on a Puhti login node (or inside an sinteractive session)
#=============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

#-----------------------------------------------------------------------------
# Colours
#-----------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
ok()   { echo -e "${GREEN}[OK]${NC}   $1"; }
err()  { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
step() { echo -e "\n${CYAN}=== $1 ===${NC}"; }

#-----------------------------------------------------------------------------
# Parse arguments
#-----------------------------------------------------------------------------
PROJECT_ID=""
INSTALL_DIR=""
DRY_RUN=false
SKIP_CONTAINERS=false
INCLUDE_LONGREADS=false

usage() {
    cat << EOF
Usage: $(basename "$0") [project_id] [options]

Arguments:
  project_id          CSC project ID (auto-detected from \$CSC_PROJECT or PWD if omitted)

Options:
  --install-dir DIR   Override install directory (default: /projappl/<project>/hla_typing)
  --dry-run           Print what would be done; make no changes
  --skip-containers   Skip container rebuild and pull steps
  --include-longreads Also pull hifihla.sif (PacBio HiFi-HLA; large download ~4 GB)
  --help, -h          Show this help
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir)      INSTALL_DIR="$2";      shift 2 ;;
        --dry-run)          DRY_RUN=true;           shift ;;
        --skip-containers)  SKIP_CONTAINERS=true;   shift ;;
        --include-longreads) INCLUDE_LONGREADS=true; shift ;;
        --help|-h)          usage ;;
        project_*)          PROJECT_ID="$1";        shift ;;
        *)  echo "Unknown argument: $1"; usage ;;
    esac
done

#-----------------------------------------------------------------------------
# Auto-detect project ID
#-----------------------------------------------------------------------------
if [[ -z "$PROJECT_ID" ]]; then
    if [[ -n "${CSC_PROJECT:-}" ]]; then
        PROJECT_ID="$CSC_PROJECT"
    elif [[ "$PWD" =~ (project_[0-9]+) ]]; then
        PROJECT_ID="${BASH_REMATCH[1]}"
    elif [[ -n "${SLURM_JOB_ACCOUNT:-}" ]]; then
        PROJECT_ID="$SLURM_JOB_ACCOUNT"
    else
        err "Could not auto-detect project ID. Pass it as: bash $(basename "$0") project_2008084"
    fi
fi

#-----------------------------------------------------------------------------
# Validate we are on Puhti (or at least an HPC with /appl and sbatch)
#-----------------------------------------------------------------------------
if [[ ! -d "/appl" ]]; then
    warn "'/appl' not found — this script is designed for CSC Puhti."
    warn "Continuing anyway, but some module commands may fail."
fi

#-----------------------------------------------------------------------------
# Resolve install directory
#-----------------------------------------------------------------------------
if [[ -z "$INSTALL_DIR" ]]; then
    # Try standard locations — projappl first, then scratch
    for candidate in \
        "/projappl/${PROJECT_ID}/hla_typing" \
        "/projappl/${PROJECT_ID}/hla_typing_pipeline" \
        "/projappl/${PROJECT_ID}/${USER}/hla_typing" \
        "/scratch/${PROJECT_ID}/${USER}/hla_typing" \
        "/scratch/${PROJECT_ID}/${USER}/hla_typing_pipeline" \
        "/scratch/${PROJECT_ID}/${USER}/local_pipeline/hla_typing_pipeline" \
        "/scratch/${PROJECT_ID}/ozcanumu/local_pipeline/hla_typing_pipeline"; do
        if [[ -d "$candidate" ]] && [[ -f "${candidate}/main.nf" ]]; then
            INSTALL_DIR="$candidate"
            break
        fi
    done
fi

if [[ -z "$INSTALL_DIR" ]]; then
    err "Could not find existing install directory. Pass --install-dir /projappl/${PROJECT_ID}/hla_typing"
fi

if [[ ! -f "${INSTALL_DIR}/main.nf" ]]; then
    err "main.nf not found in ${INSTALL_DIR} — does not look like an HLA pipeline install."
fi

CONTAINERS_DIR="${INSTALL_DIR}/containers"
USER_CONFIG="${INSTALL_DIR}/conf/user.config"

#-----------------------------------------------------------------------------
# Print banner
#-----------------------------------------------------------------------------
echo ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  HLA Typing Pipeline — Puhti In-Place Update${NC}"
echo -e "${CYAN}============================================================${NC}"
echo "  Project:       $PROJECT_ID"
echo "  Install dir:   $INSTALL_DIR"
echo "  Containers:    $CONTAINERS_DIR"
[[ "$DRY_RUN" == "true" ]] && echo -e "  ${YELLOW}Mode: DRY RUN — no changes will be made${NC}"
echo ""

#-----------------------------------------------------------------------------
# Helper: run command or print it in dry-run mode
#-----------------------------------------------------------------------------
run() {
    if [[ "$DRY_RUN" == "true" ]]; then
        echo -e "  ${YELLOW}[dry-run]${NC} $*"
    else
        "$@"
    fi
}

#-----------------------------------------------------------------------------
# Step 1: Read current version before update
#-----------------------------------------------------------------------------
step "Step 1: Current version"
OLD_VERSION=$(grep -oP "version\s*=\s*'\K[^']+" "${INSTALL_DIR}/nextflow.config" 2>/dev/null || echo "unknown")
OLD_COMMIT=$(git -C "$INSTALL_DIR" rev-parse --short HEAD 2>/dev/null || echo "unknown")
info "Current version: v${OLD_VERSION} (${OLD_COMMIT})"

#-----------------------------------------------------------------------------
# Step 2: Git pull
#-----------------------------------------------------------------------------
step "Step 2: Updating pipeline code (git pull)"

# Check for local modifications
if git -C "$INSTALL_DIR" diff --quiet HEAD 2>/dev/null; then
    : # Clean working tree
else
    warn "Local modifications detected in ${INSTALL_DIR}"
    warn "Run 'git -C ${INSTALL_DIR} diff HEAD' to inspect them."
    warn "If you want to discard them: git -C ${INSTALL_DIR} checkout -- ."
    warn "Continuing with git pull (may fail if there are conflicts)..."
fi

# Check remote URL
REMOTE_URL=$(git -C "$INSTALL_DIR" remote get-url origin 2>/dev/null || echo "")
if [[ -z "$REMOTE_URL" ]]; then
    warn "No git remote 'origin' found. Skipping git pull."
    warn "You may need to manually copy the updated files from your local machine:"
    warn "  rsync -avz --exclude='work/' --exclude='.git/' \\"
    warn "      /path/to/local/hla_typing_pipeline/ \\"
    warn "      ${INSTALL_DIR}/"
else
    info "Remote: $REMOTE_URL"
    run git -C "$INSTALL_DIR" fetch origin
    run git -C "$INSTALL_DIR" pull --ff-only origin main \
        || run git -C "$INSTALL_DIR" pull --ff-only origin master \
        || warn "git pull failed. Try: git -C ${INSTALL_DIR} pull --rebase origin main"
fi

NEW_VERSION=$(grep -oP "version\s*=\s*'\K[^']+" "${INSTALL_DIR}/nextflow.config" 2>/dev/null || echo "unknown")
NEW_COMMIT=$(git -C "$INSTALL_DIR" rev-parse --short HEAD 2>/dev/null || echo "unknown")

if [[ "$OLD_VERSION" != "$NEW_VERSION" ]]; then
    ok "Pipeline updated: v${OLD_VERSION} → v${NEW_VERSION} (${NEW_COMMIT})"
else
    info "Version unchanged (${NEW_VERSION}) — pipeline may already be up to date or update failed"
fi

#-----------------------------------------------------------------------------
# Step 3: Validate v1.4.0 files
#-----------------------------------------------------------------------------
step "Step 3: Validating v1.4.0 new files"

NEW_FILES=(
    "modules/hifihla.nf"
    "bin/generate_hml.py"
    "bin/parse_hifihla_results.py"
    "bin/hla_pipeline_metrics.py"
)
MISSING=()
for f in "${NEW_FILES[@]}"; do
    if [[ -f "${INSTALL_DIR}/${f}" ]]; then
        ok "  ${f}"
    else
        warn "  MISSING: ${f}"
        MISSING+=("$f")
    fi
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    warn "${#MISSING[@]} v1.4.0 file(s) missing. Git pull may not have worked."
    warn "Try: rsync from local machine (see above)."
fi

#-----------------------------------------------------------------------------
# Step 4: Container updates
#-----------------------------------------------------------------------------
if [[ "$SKIP_CONTAINERS" == "true" ]]; then
    info "Skipping container updates (--skip-containers)"
else
    step "Step 4: Container updates"
    mkdir -p "$CONTAINERS_DIR"

    # Load Singularity module
    module load singularity 2>/dev/null || true
    if ! command -v singularity &>/dev/null; then
        warn "singularity not found — try 'module load singularity' first"
    fi

    export SINGULARITY_CACHEDIR="${INSTALL_DIR}/singularity_cache"
    mkdir -p "$SINGULARITY_CACHEDIR"

    #------------------------------------------------------------------
    # 4a: Rebuild hla_postprocess.sif (plotly + kaleido added in v1.4.0)
    #------------------------------------------------------------------
    POSTPROCESS_DEF="${INSTALL_DIR}/containers/hla_postprocess.def"
    POSTPROCESS_SIF="${CONTAINERS_DIR}/hla_postprocess.sif"
    POSTPROCESS_STAMP="${CONTAINERS_DIR}/.hla_postprocess_md5"

    if [[ ! -f "$POSTPROCESS_DEF" ]]; then
        warn "hla_postprocess.def not found — skipping rebuild"
    else
        CURRENT_MD5=$(md5sum "$POSTPROCESS_DEF" | cut -d' ' -f1)
        STORED_MD5=$(cat "$POSTPROCESS_STAMP" 2>/dev/null || echo "none")

        if [[ "$CURRENT_MD5" == "$STORED_MD5" ]] && [[ -f "$POSTPROCESS_SIF" ]]; then
            ok "hla_postprocess.sif is up to date (MD5 matches)"
        else
            info "hla_postprocess.sif needs rebuild (def file changed or SIF missing)"
            info "Checking if we are in an interactive allocation..."

            if [[ -n "${SLURM_JOB_ID:-}" ]]; then
                # Inside interactive allocation — build directly
                info "In SLURM job ${SLURM_JOB_ID} — building directly..."
                run singularity build --fakeroot \
                    "$POSTPROCESS_SIF" \
                    "$POSTPROCESS_DEF" \
                    && echo "$CURRENT_MD5" > "$POSTPROCESS_STAMP" \
                    && ok "hla_postprocess.sif rebuilt successfully"
            else
                # Login node — submit as sbatch
                info "On login node — submitting rebuild as sbatch job..."
                REBUILD_SCRIPT="${INSTALL_DIR}/conf/rebuild_postprocess.sh"
                cat > "$REBUILD_SCRIPT" << RBEOF
#!/bin/bash
#SBATCH --job-name=rebuild_sif
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=${INSTALL_DIR}/conf/rebuild_sif_%j.out
#SBATCH --error=${INSTALL_DIR}/conf/rebuild_sif_%j.err

module load singularity 2>/dev/null || true
export SINGULARITY_CACHEDIR="${INSTALL_DIR}/singularity_cache"
singularity build --fakeroot \
    "${POSTPROCESS_SIF}" \
    "${POSTPROCESS_DEF}" \
    && echo "${CURRENT_MD5}" > "${POSTPROCESS_STAMP}" \
    && echo "[OK] hla_postprocess.sif rebuilt"
RBEOF
                if [[ "$DRY_RUN" == "true" ]]; then
                    echo -e "  ${YELLOW}[dry-run]${NC} sbatch $REBUILD_SCRIPT"
                else
                    REBUILD_JOB=$(sbatch --parsable "$REBUILD_SCRIPT")
                    ok "Rebuild job submitted: ${REBUILD_JOB}"
                    info "Monitor: tail -f ${INSTALL_DIR}/conf/rebuild_sif_${REBUILD_JOB}.out"
                    info "NOTE: Other container steps below proceed in parallel."
                fi
            fi
        fi
    fi

    #------------------------------------------------------------------
    # 4b: Pull T1K container (new in v1.4.0)
    #------------------------------------------------------------------
    T1K_SIF="${CONTAINERS_DIR}/t1k.sif"
    if [[ -f "$T1K_SIF" ]]; then
        ok "t1k.sif already present ($(du -sh "$T1K_SIF" | cut -f1))"
    else
        info "Pulling T1K container (~1.5 GB)..."
        run singularity pull "$T1K_SIF" \
            docker://quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
            && ok "t1k.sif pulled" \
            || warn "T1K pull failed — only needed for T1K typing and long-read mode"
    fi

    #------------------------------------------------------------------
    # 4c: Pull HiFi-HLA container (optional, PacBio long-read only)
    #------------------------------------------------------------------
    if [[ "$INCLUDE_LONGREADS" == "true" ]]; then
        HIFIHLA_SIF="${CONTAINERS_DIR}/hifihla.sif"
        if [[ -f "$HIFIHLA_SIF" ]]; then
            ok "hifihla.sif already present ($(du -sh "$HIFIHLA_SIF" | cut -f1))"
        else
            info "Pulling HiFi-HLA container (~4 GB)..."
            run singularity pull "$HIFIHLA_SIF" \
                docker://quay.io/pacbio/hifihla:latest \
                && ok "hifihla.sif pulled" \
                || warn "HiFi-HLA pull failed — only needed for PacBio HiFi mode"
        fi
    else
        info "Skipping HiFi-HLA container (pass --include-longreads to pull)"
    fi

    #------------------------------------------------------------------
    # 4d: Report container inventory
    #------------------------------------------------------------------
    echo ""
    echo "Container inventory (${CONTAINERS_DIR}):"
    printf "  %-35s %s\n" "Container" "Size"
    printf "  %-35s %s\n" "---------" "----"
    for sif in "${CONTAINERS_DIR}"/*.sif; do
        [[ -f "$sif" ]] || continue
        printf "  %-35s %s\n" "$(basename "$sif")" "$(du -sh "$sif" | cut -f1)"
    done
fi

#-----------------------------------------------------------------------------
# Step 5: Check conf/user.config for new v1.4.0 params
#-----------------------------------------------------------------------------
step "Step 5: Checking conf/user.config"

if [[ ! -f "$USER_CONFIG" ]]; then
    warn "conf/user.config not found at $USER_CONFIG"
    info "You may need to re-run scripts/install_puhti.sh to regenerate it, or create manually."
else
    ok "conf/user.config found"

    # New params added in v1.4.0
    declare -A NEW_PARAMS=(
        [output_format]="'text'     // text|gl_string|hml|all"
        [allele_db_version]="'3.57.0'  // IMGT/HLA release for HML output"
        [hml_center_id]="'HLA-PIPELINE' // Typing centre ID for HML XML"
        [weights_file]="null       // Path to calibrated weights JSON"
    )

    MISSING_PARAMS=()
    for param in "${!NEW_PARAMS[@]}"; do
        if grep -q "^\s*${param}\s*=" "$USER_CONFIG" 2>/dev/null; then
            ok "  ${param} already in user.config"
        else
            warn "  ${param} missing from user.config"
            MISSING_PARAMS+=("$param")
        fi
    done

    if [[ ${#MISSING_PARAMS[@]} -gt 0 ]]; then
        echo ""
        info "Add these lines to ${USER_CONFIG} (inside the params { } block):"
        for param in "${MISSING_PARAMS[@]}"; do
            echo "    ${param} = ${NEW_PARAMS[$param]}"
        done
    fi
fi

#-----------------------------------------------------------------------------
# Step 6: Verify pipeline can be imported by Nextflow (syntax check)
#-----------------------------------------------------------------------------
step "Step 6: Nextflow syntax check"

if command -v nextflow &>/dev/null || module load nextflow 2>/dev/null; then
    if [[ "$DRY_RUN" == "false" ]]; then
        cd "$INSTALL_DIR"
        if nextflow run main.nf --help > /dev/null 2>&1; then
            ok "Nextflow syntax check passed"
        else
            warn "Nextflow syntax check had warnings (run 'nextflow run main.nf --help' to inspect)"
        fi
    else
        echo -e "  ${YELLOW}[dry-run]${NC} nextflow run main.nf --help"
    fi
else
    info "Nextflow not loaded — skipping syntax check"
    info "To check: module load nextflow && nextflow run main.nf --help"
fi

#-----------------------------------------------------------------------------
# Summary
#-----------------------------------------------------------------------------
echo ""
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  Update Summary${NC}"
echo -e "${CYAN}============================================================${NC}"
echo "  Pipeline:  v${OLD_VERSION} → v${NEW_VERSION}"
echo "  Install:   ${INSTALL_DIR}"
echo ""
echo "  Next steps:"
echo ""
echo "  1. If hla_postprocess.sif rebuild was submitted as a job:"
echo "     squeue -u \$USER   # wait for rebuild job to finish"
echo ""
echo "  2. Run the 1KGP calibration analysis (50 samples, 5 superpopulations):"
echo "     bash ${INSTALL_DIR}/scripts/submit_calibration_puhti.sh \\"
echo "         --project ${PROJECT_ID} \\"
echo "         --source 30x \\"
echo "         --n-samples 50 \\"
echo "         --tools hlahd,spechla,optitype,arcashla"
echo ""
echo "  3. Or run a single clinical sample:"
echo "     bash ${INSTALL_DIR}/scripts/run_hla_puhti.sh \\"
echo "         --bam /scratch/${PROJECT_ID}/sample.bam \\"
echo "         --tools hlahd,spechla,optitype,arcashla"
echo -e "${CYAN}============================================================${NC}"
echo ""
