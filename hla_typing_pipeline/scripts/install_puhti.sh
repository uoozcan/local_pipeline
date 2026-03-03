#!/bin/bash
#=============================================================================
# HLA Typing Pipeline - Installation Script for CSC Puhti
#=============================================================================
# Following CSC best practices:
#   - https://docs.csc.fi/computing/installing/
#   - https://docs.csc.fi/computing/containers/tykky/
#
# Usage:
#   bash install_puhti.sh <project_id>
#
# Example:
#   bash install_puhti.sh project_2001234
#
#=============================================================================

set -e

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
log() { echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }
info() { echo -e "${BLUE}[INFO]${NC} $1"; }

show_help() {
    cat << EOF
HLA Typing Pipeline - Installation for CSC Puhti

Usage:
    $(basename "$0") <project_id> [options]

Arguments:
    project_id          Your CSC project ID (e.g., project_2001234)

Options:
    --install-dir DIR   Installation directory (default: /projappl/<project>/hla_typing)
    --scratch-dir DIR   Scratch directory for outputs (default: /scratch/<project>/hla_results)
    --skip-containers   Skip container building (use existing)
    --skip-postprocess  Skip post-processing container
    --help, -h          Show this help

Examples:
    # Standard installation
    $(basename "$0") project_2001234

    # Custom directories
    $(basename "$0") project_2001234 --install-dir /projappl/project_2001234/mytools

Notes:
    - Run this script on Puhti login node
    - Container building may require an interactive session for large builds
    - First run takes ~30-60 minutes depending on downloads

EOF
    exit 0
}

#-----------------------------------------------------------------------------
# Parse Arguments
#-----------------------------------------------------------------------------
PROJECT_ID=""
INSTALL_DIR=""
SCRATCH_DIR=""
SKIP_CONTAINERS=false
SKIP_POSTPROCESS=false

# First argument should be project ID
if [[ $# -lt 1 ]] || [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
    show_help
fi

PROJECT_ID="$1"
shift

# Validate project ID format
if [[ ! "$PROJECT_ID" =~ ^project_[0-9]+$ ]]; then
    error "Invalid project ID format. Expected: project_XXXXXXX (e.g., project_2001234)"
fi

# Parse remaining options
while [[ $# -gt 0 ]]; do
    case $1 in
        --install-dir) INSTALL_DIR="$2"; shift 2 ;;
        --scratch-dir) SCRATCH_DIR="$2"; shift 2 ;;
        --skip-containers) SKIP_CONTAINERS=true; shift ;;
        --skip-postprocess) SKIP_POSTPROCESS=true; shift ;;
        --help|-h) show_help ;;
        *) error "Unknown option: $1" ;;
    esac
done

# Set default directories
INSTALL_DIR="${INSTALL_DIR:-/projappl/${PROJECT_ID}/hla_typing}"
SCRATCH_DIR="${SCRATCH_DIR:-/scratch/${PROJECT_ID}/hla_results}"
CONTAINER_DIR="${INSTALL_DIR}/containers"
LOG_FILE="${INSTALL_DIR}/install_${TIMESTAMP}.log"

#-----------------------------------------------------------------------------
# Pre-flight Checks
#-----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "HLA Typing Pipeline Installer for CSC Puhti"
echo "=============================================="
echo ""
log "Starting installation..."
info "Project:      $PROJECT_ID"
info "Install dir:  $INSTALL_DIR"
info "Scratch dir:  $SCRATCH_DIR"
info "Pipeline src: $PIPELINE_DIR"
echo ""

# Check if on Puhti
if [[ ! -d "/appl" ]]; then
    warn "This doesn't appear to be CSC Puhti"
    warn "Some features may not work correctly"
fi

# Check project directory access
if [[ ! -d "/projappl/${PROJECT_ID}" ]]; then
    error "Cannot access /projappl/${PROJECT_ID}. Check your project ID and permissions."
fi

#-----------------------------------------------------------------------------
# Create Directory Structure
#-----------------------------------------------------------------------------
log "Creating directory structure..."

mkdir -p "$INSTALL_DIR"
mkdir -p "$CONTAINER_DIR"
mkdir -p "$SCRATCH_DIR"
mkdir -p "${INSTALL_DIR}/singularity_cache"
mkdir -p "${INSTALL_DIR}/logs"

# Start logging
exec > >(tee -a "$LOG_FILE") 2>&1
log "Log file: $LOG_FILE"

#-----------------------------------------------------------------------------
# Load Required Modules
#-----------------------------------------------------------------------------
log "Loading Puhti modules..."

module purge 2>/dev/null || true

# Load modules (with fallbacks for different Puhti configurations)
module load gcc/11.3.0 2>/dev/null || module load gcc 2>/dev/null || warn "gcc module not loaded"
module load cmake/3.23.1 2>/dev/null || module load cmake 2>/dev/null || warn "cmake module not loaded"
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || warn "singularity/apptainer not loaded"
module load nextflow/23.04.1 2>/dev/null || module load nextflow 2>/dev/null || warn "nextflow module not loaded"
module load tykky 2>/dev/null || info "Tykky module not available (optional)"

info "Modules loaded"

#-----------------------------------------------------------------------------
# Copy Pipeline Files
#-----------------------------------------------------------------------------
log "Copying pipeline files..."

# Copy pipeline to install directory
if [[ "$PIPELINE_DIR" != "$INSTALL_DIR/pipeline" ]]; then
    mkdir -p "${INSTALL_DIR}/pipeline"
    rsync -av --exclude='.git' --exclude='work' --exclude='results' --exclude='.nextflow*' \
        "${PIPELINE_DIR}/" "${INSTALL_DIR}/pipeline/"
    info "Pipeline copied to ${INSTALL_DIR}/pipeline"
fi

#-----------------------------------------------------------------------------
# Build/Download Containers
#-----------------------------------------------------------------------------
if [[ "$SKIP_CONTAINERS" != true ]]; then
    log "Setting up containers..."

    # Set Singularity cache directory
    export SINGULARITY_CACHEDIR="${INSTALL_DIR}/singularity_cache"
    export APPTAINER_CACHEDIR="${INSTALL_DIR}/singularity_cache"

    # 1. ArcasHLA container (from biocontainers)
    log "Pulling arcasHLA container..."
    if [[ ! -f "${CONTAINER_DIR}/arcashla.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/arcashla.sif" \
            docker://quay.io/biocontainers/arcas-hla:0.5.0--hdfd78af_1 2>&1 || \
            warn "Failed to pull arcasHLA container"
    else
        info "arcasHLA container already exists"
    fi

    # 2. OptiType container
    log "Pulling OptiType container..."
    if [[ ! -f "${CONTAINER_DIR}/optitype.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/optitype.sif" \
            docker://fred2/optitype:1.3.5 2>&1 || \
            warn "Failed to pull OptiType container"
    else
        info "OptiType container already exists"
    fi

    # 3. FastQC container
    log "Pulling FastQC container..."
    if [[ ! -f "${CONTAINER_DIR}/fastqc.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/fastqc.sif" \
            docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0 2>&1 || \
            warn "Failed to pull FastQC container"
    else
        info "FastQC container already exists"
    fi

    # 4. xHLA container
    log "Pulling xHLA container..."
    if [[ ! -f "${CONTAINER_DIR}/xhla.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/xhla.sif" \
            docker://quay.io/biocontainers/xhla:latest 2>&1 || \
            warn "Failed to pull xHLA container"
    else
        info "xHLA container already exists"
    fi

    # 5. HLA*LA container
    log "Pulling HLA*LA container..."
    if [[ ! -f "${CONTAINER_DIR}/hlala.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/hlala.sif" \
            docker://quay.io/biocontainers/hla-la:latest 2>&1 || \
            warn "Failed to pull HLA*LA container"
    else
        info "HLA*LA container already exists"
    fi

    # 6. BAMQC container
    log "Pulling BAMQC container..."
    if [[ ! -f "${CONTAINER_DIR}/bamqc.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/bamqc.sif" \
            docker://kennethlim206/bamqc:latest 2>&1 || \
            warn "Failed to pull BAMQC container"
    else
        info "BAMQC container already exists"
    fi

    # 7. flow-OptiType container
    log "Pulling flow-OptiType container..."
    if [[ ! -f "${CONTAINER_DIR}/flow_optitype.sif" ]]; then
        singularity pull "${CONTAINER_DIR}/flow_optitype.sif" \
            docker://nmdpbioinformatics/flow-optitype:latest 2>&1 || \
            warn "Failed to pull flow-OptiType container"
    else
        info "flow-OptiType container already exists"
    fi

    # 8. SpecHLA with SpecHap container
    log "Building SpecHLA container with SpecHap..."
    if [[ ! -f "${CONTAINER_DIR}/spechla_with_spechap.sif" ]]; then
        SPECHLA_DEF="${INSTALL_DIR}/pipeline/containers/spechla_with_spechap.def"
        if [[ -f "$SPECHLA_DEF" ]]; then
            info "Building from definition file (this may take 15-30 minutes)..."
            info "If this fails on login node, run in an interactive session:"
            info "  sinteractive --account=${PROJECT_ID} --mem=16G --time=02:00:00"

            singularity build --fakeroot "${CONTAINER_DIR}/spechla_with_spechap.sif" "$SPECHLA_DEF" 2>&1 || {
                warn "Container build failed with --fakeroot"
                info "Trying alternative: pulling from registry..."
                singularity pull "${CONTAINER_DIR}/spechla_with_spechap.sif" \
                    docker://quay.io/biocontainers/spechla:1.0.7--hdfd78af_3 2>&1 || {
                    warn "Could not pull pre-built container"
                    info ""
                    info "To build SpecHLA container manually:"
                    info "  1. Start interactive session: sinteractive --account=${PROJECT_ID} --mem=16G"
                    info "  2. Load modules: module load singularity"
                    info "  3. Build: singularity build --fakeroot ${CONTAINER_DIR}/spechla_with_spechap.sif $SPECHLA_DEF"
                    info ""
                    info "Or copy from local machine:"
                    info "  scp spechla_with_spechap.sif ${USER}@puhti.csc.fi:${CONTAINER_DIR}/"
                }
            }
        else
            warn "SpecHLA definition file not found at $SPECHLA_DEF"
        fi
    else
        info "SpecHLA container already exists"
    fi

    # 9. Post-processing container with matplotlib
    if [[ "$SKIP_POSTPROCESS" != true ]]; then
        log "Building post-processing container..."
        if [[ ! -f "${CONTAINER_DIR}/hla_postprocess.sif" ]]; then
            POSTPROCESS_DEF="${INSTALL_DIR}/pipeline/containers/hla_postprocess.def"
            if [[ -f "$POSTPROCESS_DEF" ]]; then
                singularity build --fakeroot "${CONTAINER_DIR}/hla_postprocess.sif" "$POSTPROCESS_DEF" 2>&1 || {
                    warn "Post-processing container build failed"
                    info "Trying to create with Tykky..."

                    # Alternative: Use Tykky to create Python environment
                    if command -v pip-containerize &>/dev/null; then
                        log "Creating visualization environment with Tykky..."

                        cat > /tmp/hla_viz_requirements.txt << 'REQEOF'
numpy>=1.21.0
pandas>=1.4.0
matplotlib>=3.5.0
seaborn>=0.11.0
jinja2>=3.0.0
multiqc>=1.14
REQEOF

                        pip-containerize new --prefix "${INSTALL_DIR}/viz_env" /tmp/hla_viz_requirements.txt || \
                            warn "Tykky installation also failed"
                        rm -f /tmp/hla_viz_requirements.txt
                    fi
                }
            else
                warn "Post-processing definition file not found"
            fi
        else
            info "Post-processing container already exists"
        fi
    fi

    info "Container setup complete"
fi

#-----------------------------------------------------------------------------
# Create Configuration Files
#-----------------------------------------------------------------------------
log "Creating configuration files..."

# Create user configuration for Puhti
cat > "${INSTALL_DIR}/pipeline/conf/user.config" << EOF
/*
 * HLA Typing Pipeline - User Configuration for Puhti
 * Generated: $(date)
 * Project: ${PROJECT_ID}
 */

params {
    // CSC Project
    project = '${PROJECT_ID}'

    // Paths
    install_dir   = '${INSTALL_DIR}'
    container_dir = '${CONTAINER_DIR}'
    outdir        = '${SCRATCH_DIR}'

    // SpecHLA settings
    spechla_path     = '/opt/SpecHLA'
    use_local_spechla = false

    // Default tools
    tools = 'spechla,arcashla'

    // Reference genome
    reference = 'hg38'
}

// Include Puhti-specific settings
includeConfig 'puhti.config'

// Container overrides for user installation
process {
    withName: 'SPECHLA|SPECHLA_FASTQ' {
        container = '${CONTAINER_DIR}/spechla_with_spechap.sif'
    }
    withName: 'ARCASHLA|ARCASHLA_FASTQ' {
        container = '${CONTAINER_DIR}/arcashla.sif'
    }
    withName: 'OPTITYPE|OPTITYPE_FASTQ' {
        container = '${CONTAINER_DIR}/optitype.sif'
    }
    withName: 'XHLA|XHLA_FASTQ' {
        container = '${CONTAINER_DIR}/xhla.sif'
    }
    withName: 'HLALA' {
        container = '${CONTAINER_DIR}/hlala.sif'
        containerOptions = "--bind ${CONTAINER_DIR}/../hlala_graphs:/usr/local/bin/HLA-LA/graphs"
    }
    withName: 'BAMQC|BAMQC_FASTQ' {
        container = '${CONTAINER_DIR}/bamqc.sif'
    }
    withName: 'FLOW_OPTITYPE' {
        container = '${CONTAINER_DIR}/flow_optitype.sif'
    }
    withName: 'FASTQC_BAM|FASTQC_FASTQ' {
        container = '${CONTAINER_DIR}/fastqc.sif'
    }
    withName: 'QC_BAM|QC_FASTQ' {
        container = '${CONTAINER_DIR}/arcashla.sif'
    }
    withName: 'CONSENSUS|HLA_VISUALIZE|HLA_SUMMARY_REPORT|MULTIQC' {
        container = '${CONTAINER_DIR}/hla_postprocess.sif'
    }
}
EOF

info "Configuration created: ${INSTALL_DIR}/pipeline/conf/user.config"

#-----------------------------------------------------------------------------
# Create Convenience Scripts
#-----------------------------------------------------------------------------
log "Creating convenience scripts..."

# Environment loader script
cat > "${INSTALL_DIR}/load_env.sh" << EOF
#!/bin/bash
# Load HLA Typing Pipeline environment
# Usage: source ${INSTALL_DIR}/load_env.sh

# Load required modules
module purge
module load singularity
module load nextflow

# Set environment variables
export HLA_PIPELINE="${INSTALL_DIR}/pipeline"
export HLA_CONTAINERS="${CONTAINER_DIR}"
export HLA_PROJECT="${PROJECT_ID}"
export SINGULARITY_CACHEDIR="${INSTALL_DIR}/singularity_cache"

# Add pipeline bin to PATH
export PATH="\${HLA_PIPELINE}/bin:\${PATH}"

# Convenience function
run_hla() {
    nextflow run \${HLA_PIPELINE}/main.nf \\
        -c \${HLA_PIPELINE}/conf/user.config \\
        -profile singularity \\
        "\$@"
}

echo "HLA Typing Pipeline environment loaded"
echo "  Pipeline: \${HLA_PIPELINE}"
echo "  Project:  ${PROJECT_ID}"
echo ""
echo "Usage: run_hla --input_bam /path/to/sample.bam"
EOF
chmod +x "${INSTALL_DIR}/load_env.sh"

# SLURM batch script template
cat > "${INSTALL_DIR}/submit_hla.sh" << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_typing
#SBATCH --account=PROJECTID
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=hla_%j.out
#SBATCH --error=hla_%j.err

# Load environment
source INSTALLDIR/load_env.sh

# Parse arguments
INPUT_BAM=""
INPUT_FASTQ1=""
INPUT_FASTQ2=""
OUTDIR="SCRATCHDIR"
TOOLS="spechla,arcashla"

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam) INPUT_BAM="$2"; shift 2 ;;
        --fastq1) INPUT_FASTQ1="$2"; shift 2 ;;
        --fastq2) INPUT_FASTQ2="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --tools) TOOLS="$2"; shift 2 ;;
        *) shift ;;
    esac
done

# Build and run command
CMD="run_hla --outdir ${OUTDIR} --tools ${TOOLS}"

if [[ -n "$INPUT_BAM" ]]; then
    CMD+=" --input_bam ${INPUT_BAM}"
elif [[ -n "$INPUT_FASTQ1" ]]; then
    CMD+=" --input_fastq_1 ${INPUT_FASTQ1} --input_fastq_2 ${INPUT_FASTQ2}"
else
    echo "Error: Provide --bam or --fastq1/--fastq2"
    exit 1
fi

echo "Running: $CMD"
eval $CMD
EOF

# Replace placeholders
sed -i "s|PROJECTID|${PROJECT_ID}|g" "${INSTALL_DIR}/submit_hla.sh"
sed -i "s|INSTALLDIR|${INSTALL_DIR}|g" "${INSTALL_DIR}/submit_hla.sh"
sed -i "s|SCRATCHDIR|${SCRATCH_DIR}|g" "${INSTALL_DIR}/submit_hla.sh"
chmod +x "${INSTALL_DIR}/submit_hla.sh"

info "Scripts created"

#-----------------------------------------------------------------------------
# Verify Installation
#-----------------------------------------------------------------------------
log "Verifying installation..."

echo ""
echo "=============================================="
echo "Installation Summary"
echo "=============================================="

# Check pipeline
if [[ -f "${INSTALL_DIR}/pipeline/main.nf" ]]; then
    echo -e "${GREEN}[OK]${NC} Pipeline files installed"
else
    echo -e "${RED}[FAIL]${NC} Pipeline files missing"
fi

# Check containers
for container in arcashla optitype fastqc xhla hlala bamqc flow_optitype spechla_with_spechap hla_postprocess; do
    if [[ -f "${CONTAINER_DIR}/${container}.sif" ]]; then
        echo -e "${GREEN}[OK]${NC} Container: ${container}.sif"
    else
        echo -e "${YELLOW}[WARN]${NC} Container missing: ${container}.sif"
    fi
done

# Check configuration
if [[ -f "${INSTALL_DIR}/pipeline/conf/user.config" ]]; then
    echo -e "${GREEN}[OK]${NC} User configuration"
else
    echo -e "${RED}[FAIL]${NC} User configuration missing"
fi

echo ""
echo "=============================================="
echo "Usage Instructions"
echo "=============================================="
echo ""
echo "1. Load the environment:"
echo "   source ${INSTALL_DIR}/load_env.sh"
echo ""
echo "2. Run on a single BAM file (interactive):"
echo "   run_hla --input_bam /path/to/sample.bam --outdir ./results"
echo ""
echo "3. Submit as batch job:"
echo "   sbatch ${INSTALL_DIR}/submit_hla.sh --bam /path/to/sample.bam"
echo ""
echo "4. For multiple samples, create a samplesheet (CSV):"
echo "   sample_id,bam_path"
echo "   sample1,/path/to/sample1.bam"
echo "   sample2,/path/to/sample2.bam"
echo ""
echo "   Then run:"
echo "   run_hla --input_samplesheet samples.csv"
echo ""
echo "=============================================="
echo "Log file: $LOG_FILE"
echo "=============================================="

log "Installation complete!"
