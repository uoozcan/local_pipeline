#!/bin/bash
#=============================================================================
# HLA Typing Pipeline - Full Setup Script for CSC Puhti
#=============================================================================
# This script automates the installation of all HLA typing tools and
# dependencies on CSC Puhti supercomputer.
#
# Usage:
#   ./puhti_full_setup.sh [--project PROJECT_ID] [--install-dir DIR]
#
# Example:
#   ./puhti_full_setup.sh --project project_2001234
#   ./puhti_full_setup.sh --project project_2001234 --install-dir /scratch/project_2001234/myuser/hla_analysis
#
#=============================================================================

set -e  # Exit on error

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="setup_log_${TIMESTAMP}.txt"

# Default values
PROJECT_ID=""
INSTALL_DIR=""
SKIP_SPECHLA=false
SKIP_HLAHD=false
SKIP_ARCASHLA=false
SKIP_OPTITYPE=false
SKIP_HLALA=false
SKIP_XHLA=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

#-----------------------------------------------------------------------------
# Helper Functions
#-----------------------------------------------------------------------------
log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
    exit 1
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

check_command() {
    if command -v "$1" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

#-----------------------------------------------------------------------------
# Parse Arguments
#-----------------------------------------------------------------------------
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --project)
                PROJECT_ID="$2"
                shift 2
                ;;
            --install-dir)
                INSTALL_DIR="$2"
                shift 2
                ;;
            --skip-spechla)
                SKIP_SPECHLA=true
                shift
                ;;
            --skip-hlahd)
                SKIP_HLAHD=true
                shift
                ;;
            --skip-arcashla)
                SKIP_ARCASHLA=true
                shift
                ;;
            --skip-optitype)
                SKIP_OPTITYPE=true
                shift
                ;;
            --skip-hlala)
                SKIP_HLALA=true
                shift
                ;;
            --skip-xhla)
                SKIP_XHLA=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                error "Unknown option: $1"
                ;;
        esac
    done
}

show_help() {
    cat << EOF
HLA Typing Pipeline - Full Setup Script for CSC Puhti

Usage:
    $0 [OPTIONS]

Options:
    --project PROJECT_ID    CSC project ID (e.g., project_2001234)
    --install-dir DIR       Installation directory (default: /scratch/PROJECT/USER/hla_analysis)
    --skip-spechla          Skip SpecHLA installation
    --skip-hlahd            Skip HLA-HD installation
    --skip-arcashla         Skip arcasHLA installation
    --skip-optitype         Skip OptiType installation
    --skip-hlala            Skip HLA*LA installation
    --skip-xhla             Skip xHLA installation
    -h, --help              Show this help message

Examples:
    # Full installation
    $0 --project project_2001234

    # Install only SpecHLA and HLA-HD
    $0 --project project_2001234 --skip-arcashla --skip-optitype --skip-hlala

    # Custom installation directory
    $0 --project project_2001234 --install-dir /projappl/project_2001234/hla_tools
EOF
}

#-----------------------------------------------------------------------------
# Detect Environment
#-----------------------------------------------------------------------------
detect_environment() {
    log "Detecting environment..."

    # Check if on Puhti
    if [[ $(hostname) == *"puhti"* ]] || [[ -d "/appl/soft" ]]; then
        info "Running on CSC Puhti"
        ON_PUHTI=true
    else
        warn "Not running on Puhti - some features may not work"
        ON_PUHTI=false
    fi

    # Auto-detect project ID if not provided
    if [[ -z "$PROJECT_ID" ]]; then
        if [[ "$PWD" =~ project_[0-9]+ ]]; then
            PROJECT_ID=$(echo "$PWD" | grep -o 'project_[0-9]*')
            info "Auto-detected project ID: $PROJECT_ID"
        else
            error "Could not detect project ID. Please provide with --project"
        fi
    fi

    # Set default installation directory
    if [[ -z "$INSTALL_DIR" ]]; then
        INSTALL_DIR="/scratch/${PROJECT_ID}/${USER}/hla_analysis"
    fi

    info "Installation directory: $INSTALL_DIR"
}

#-----------------------------------------------------------------------------
# Load Puhti Modules
#-----------------------------------------------------------------------------
load_modules() {
    log "Loading required modules..."

    if [[ "$ON_PUHTI" == true ]]; then
        module purge 2>/dev/null || true

        # Core modules
        module load gcc/11.3.0 2>/dev/null || module load gcc
        module load cmake/3.23.1 2>/dev/null || module load cmake
        module load python-data/3.10-22.09 2>/dev/null || module load python

        # Bioinformatics tools
        module load samtools/1.15 2>/dev/null || module load samtools
        module load bwa/0.7.17 2>/dev/null || module load bwa
        module load nextflow/22.10.6 2>/dev/null || module load nextflow

        # Container runtime
        module load singularity 2>/dev/null || true

        info "Modules loaded successfully"
    else
        warn "Not on Puhti - skipping module loading"
        warn "Make sure gcc, cmake, samtools, bwa, python3 are available"
    fi
}

#-----------------------------------------------------------------------------
# Create Directory Structure
#-----------------------------------------------------------------------------
create_directories() {
    log "Creating directory structure..."

    mkdir -p "$INSTALL_DIR"
    cd "$INSTALL_DIR"

    # Create subdirectories
    mkdir -p hla_typing_pipeline
    mkdir -p hla_references/containers
    mkdir -p hla_references/databases/hlahd_db
    mkdir -p hla_references/databases/hlala_graphs
    mkdir -p logs

    info "Directory structure created at $INSTALL_DIR"
}

#-----------------------------------------------------------------------------
# Clone Pipeline Repository
#-----------------------------------------------------------------------------
clone_pipeline() {
    log "Cloning HLA typing pipeline..."

    cd "$INSTALL_DIR"

    if [[ -d "hla_typing_pipeline/.git" ]]; then
        info "Pipeline already exists, pulling latest changes..."
        cd hla_typing_pipeline
        git pull origin main || warn "Could not pull latest changes"
        cd ..
    else
        rm -rf hla_typing_pipeline
        git clone https://github.com/uoozcan/local_pipeline.git hla_typing_pipeline
        info "Pipeline cloned successfully"
    fi
}

#-----------------------------------------------------------------------------
# Install SpecHLA (Container-based)
#-----------------------------------------------------------------------------
install_spechla() {
    if [[ "$SKIP_SPECHLA" == true ]]; then
        info "Skipping SpecHLA installation"
        return
    fi

    log "Setting up SpecHLA (container-based)..."

    cd "$INSTALL_DIR"

    # SpecHLA via Singularity container (recommended for Puhti)
    # This avoids complex local builds of SpecHap with ARPACK dependencies
    if check_command singularity; then
        log "Pulling SpecHLA container..."

        # First try the pre-built container from a registry
        singularity pull "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" \
            docker://quay.io/biocontainers/spechla:1.0.7--hdfd78af_3 2>&1 | tee -a "$LOG_FILE" || {

            warn "Could not pull SpecHLA container from biocontainers"
            info "Trying alternative source..."

            # Try alternative sources or provide manual instructions
            singularity pull "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" \
                docker://deepomicslab/spechla:latest 2>&1 | tee -a "$LOG_FILE" || {

                warn "Could not pull SpecHLA container from dockerhub"
                info ""
                info "SpecHLA container needs to be copied manually."
                info "If you have the container file locally, copy it to:"
                info "  $INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif"
                info ""
                info "Alternative: Build container from definition file:"
                info "  singularity build spechla_1.0.7-3.sif spechla.def"
            }
        }
    else
        warn "Singularity not available - cannot install SpecHLA container"
        info "Load singularity module and re-run setup"
    fi

    # Verify container
    if [[ -f "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" ]]; then
        info "SpecHLA container available"

        # Test container
        log "Testing SpecHLA container..."
        singularity exec "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" \
            ls /opt/SpecHLA/script/whole/SpecHLA.sh 2>/dev/null && \
            info "SpecHLA container verified - SpecHLA.sh found" || \
            warn "SpecHLA container may have issues - SpecHLA.sh not found at expected path"
    else
        warn "SpecHLA container not found"
        info "The pipeline will not be able to run SpecHLA without the container."
        info ""
        info "To obtain the container, either:"
        info "  1. Copy from local machine: scp spechla_1.0.7-3.sif user@puhti.csc.fi:$INSTALL_DIR/hla_references/containers/"
        info "  2. Build from definition file (requires fakeroot): singularity build --fakeroot spechla_1.0.7-3.sif spechla.def"
    fi

    cd "$INSTALL_DIR"
}

#-----------------------------------------------------------------------------
# Install HLA-HD
#-----------------------------------------------------------------------------
install_hlahd() {
    if [[ "$SKIP_HLAHD" == true ]]; then
        info "Skipping HLA-HD installation"
        return
    fi

    log "Setting up HLA-HD..."

    cd "$INSTALL_DIR"

    # HLA-HD requires a license and manual download
    # We'll set up the container approach

    info "HLA-HD requires manual download due to licensing"
    info "Options:"
    info "  1. Download from: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/"
    info "  2. Use Singularity container if available"

    # Create placeholder script
    cat > "$INSTALL_DIR/hla_references/setup_hlahd.sh" << 'HLAHD_SETUP'
#!/bin/bash
# HLA-HD Setup Instructions
#
# 1. Register and download HLA-HD from:
#    https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
#
# 2. Extract to this directory:
#    tar -xzf hlahd.*.tar.gz -C /path/to/hla_references/databases/hlahd_db/
#
# 3. Update the IMGT database:
#    cd hlahd.*/
#    sh update.dictionary.sh
#
# 4. If using Singularity container:
#    singularity pull hla_references/containers/hlahd.sif docker://quay.io/biocontainers/hla-hd:1.4.0--hdfd78af_0

echo "Please follow the instructions above to set up HLA-HD"
HLAHD_SETUP
    chmod +x "$INSTALL_DIR/hla_references/setup_hlahd.sh"

    # Try to pull container if available
    if check_command singularity; then
        log "Attempting to pull HLA-HD container..."
        singularity pull "$INSTALL_DIR/hla_references/containers/hlahd.sif" \
            docker://quay.io/biocontainers/hla-hd:1.4.0--hdfd78af_0 2>/dev/null || \
            warn "Could not pull HLA-HD container - may need manual setup"
    fi
}

#-----------------------------------------------------------------------------
# Install arcasHLA
#-----------------------------------------------------------------------------
install_arcashla() {
    if [[ "$SKIP_ARCASHLA" == true ]]; then
        info "Skipping arcasHLA installation"
        return
    fi

    log "Installing arcasHLA..."

    cd "$INSTALL_DIR"

    # Use Singularity container (recommended for Puhti)
    if check_command singularity; then
        log "Pulling arcasHLA container..."
        singularity pull "$INSTALL_DIR/hla_references/containers/arcashla.sif" \
            docker://quay.io/biocontainers/arcas-hla:0.5.0--hdfd78af_1 2>&1 | tee -a "$LOG_FILE" || \
            warn "Could not pull arcasHLA container"
    else
        warn "Singularity not available - cannot install arcasHLA"
    fi

    # Verify
    if [[ -f "$INSTALL_DIR/hla_references/containers/arcashla.sif" ]]; then
        info "arcasHLA container available"
    else
        warn "arcasHLA installation incomplete"
        info "Try manually: singularity pull arcashla.sif docker://quay.io/biocontainers/arcas-hla:0.5.0--hdfd78af_1"
    fi
}

#-----------------------------------------------------------------------------
# Install OptiType
#-----------------------------------------------------------------------------
install_optitype() {
    if [[ "$SKIP_OPTITYPE" == true ]]; then
        info "Skipping OptiType installation"
        return
    fi

    log "Installing OptiType..."

    cd "$INSTALL_DIR"

    # OptiType via Singularity container
    if check_command singularity; then
        log "Pulling OptiType container..."
        singularity pull "$INSTALL_DIR/hla_references/containers/optitype.sif" \
            docker://fred2/optitype:1.3.5 2>&1 | tee -a "$LOG_FILE" || \
            warn "Could not pull OptiType container"
    fi

    # Verify
    if [[ -f "$INSTALL_DIR/hla_references/containers/optitype.sif" ]]; then
        info "OptiType container available"
    else
        warn "OptiType installation incomplete"
        info "Try: singularity pull optitype.sif docker://fred2/optitype:1.3.5"
    fi
}

#-----------------------------------------------------------------------------
# Install HLA*LA
#-----------------------------------------------------------------------------
install_hlala() {
    if [[ "$SKIP_HLALA" == true ]]; then
        info "Skipping HLA*LA installation"
        return
    fi

    log "Installing HLA*LA..."

    cd "$INSTALL_DIR"

    # HLA*LA via Singularity container
    if check_command singularity; then
        log "Pulling HLA*LA container..."
        singularity pull "$INSTALL_DIR/hla_references/containers/hlala.sif" \
            docker://quay.io/biocontainers/hla-la:1.0.3--hd03093a_0 2>&1 | tee -a "$LOG_FILE" || \
            warn "Could not pull HLA*LA container"
    fi

    # Download HLA*LA graphs
    log "Downloading HLA*LA reference graphs..."
    mkdir -p "$INSTALL_DIR/hla_references/databases/hlala_graphs"
    cd "$INSTALL_DIR/hla_references/databases/hlala_graphs"

    if [[ ! -d "PRG_MHC_GRCh38_withIMGT" ]]; then
        # Download from HLA*LA repository
        log "Downloading PRG_MHC_GRCh38_withIMGT (this may take a while - ~30GB)..."
        wget "http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz" 2>&1 | tee -a "$INSTALL_DIR/$LOG_FILE" || \
            warn "Could not download HLA*LA graphs - may need manual download"

        if [[ -f "PRG_MHC_GRCh38_withIMGT.tar.gz" ]]; then
            log "Extracting HLA*LA graphs..."
            tar -xzf PRG_MHC_GRCh38_withIMGT.tar.gz
            rm -f PRG_MHC_GRCh38_withIMGT.tar.gz
            info "HLA*LA graphs downloaded and extracted"
        fi
    else
        info "HLA*LA graphs already exist"
    fi

    # Index the graphs if not already done (required for first run)
    if [[ -d "PRG_MHC_GRCh38_withIMGT" ]] && [[ ! -f "PRG_MHC_GRCh38_withIMGT/serializedGRAPH" ]]; then
        log "HLA*LA graphs need indexing on first run - this will happen automatically"
        info "Note: First HLA*LA run will take longer due to graph indexing"
    fi

    cd "$INSTALL_DIR"

    # Create HLA*LA working directory structure
    mkdir -p "$INSTALL_DIR/hla_references/hlala_workdir"

    # Verify
    if [[ -f "$INSTALL_DIR/hla_references/containers/hlala.sif" ]]; then
        info "HLA*LA container available"
    else
        warn "HLA*LA container not available"
    fi

    if [[ -d "$INSTALL_DIR/hla_references/databases/hlala_graphs/PRG_MHC_GRCh38_withIMGT" ]]; then
        info "HLA*LA reference graphs available"
    else
        warn "HLA*LA reference graphs not found"
        info "Download manually: wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz"
    fi
}

#-----------------------------------------------------------------------------
# Install xHLA
#-----------------------------------------------------------------------------
install_xhla() {
    if [[ "$SKIP_XHLA" == true ]]; then
        info "Skipping xHLA installation"
        return
    fi

    log "Installing xHLA..."

    cd "$INSTALL_DIR"

    # xHLA via Singularity container
    if check_command singularity; then
        log "Pulling xHLA container..."
        singularity pull "$INSTALL_DIR/hla_references/containers/xhla.sif" \
            docker://humanlongevity/hla:latest 2>&1 | tee -a "$LOG_FILE" || \
            warn "Could not pull xHLA container from docker hub"

        # Alternative: try biocontainers
        if [[ ! -f "$INSTALL_DIR/hla_references/containers/xhla.sif" ]]; then
            singularity pull "$INSTALL_DIR/hla_references/containers/xhla.sif" \
                docker://quay.io/biocontainers/xhla:1.0--hdfd78af_0 2>&1 | tee -a "$LOG_FILE" || \
                warn "Could not pull xHLA container"
        fi
    fi

    # Verify
    if [[ -f "$INSTALL_DIR/hla_references/containers/xhla.sif" ]]; then
        info "xHLA container available"
    else
        warn "xHLA installation incomplete"
        info "xHLA container may need to be copied manually if not available from docker hub"
    fi
}

#-----------------------------------------------------------------------------
# Install Python Dependencies
#-----------------------------------------------------------------------------
install_python_deps() {
    log "Installing Python dependencies..."

    cd "$INSTALL_DIR"

    # Create requirements file for pipeline visualization only
    # Note: SpecHLA dependencies are inside the container
    cat > requirements.txt << 'EOF'
# Pipeline visualization dependencies
matplotlib>=3.5.0
pandas>=1.3.0
seaborn>=0.11.0
jinja2>=3.0.0
numpy>=1.21.0
EOF

    # Install with pip
    pip install --user -r requirements.txt 2>&1 | tee -a "$LOG_FILE" || \
        warn "Some Python packages may not have installed correctly"

    # Verify visualization packages
    python3 -c "import matplotlib; import pandas; import numpy" 2>/dev/null && \
        info "Python visualization dependencies verified" || \
        warn "Some Python visualization dependencies may be missing"

    info "Python dependencies installed"
}

#-----------------------------------------------------------------------------
# Configure Pipeline
#-----------------------------------------------------------------------------
configure_pipeline() {
    log "Configuring pipeline..."

    cd "$INSTALL_DIR/hla_typing_pipeline"

    # Create user configuration file
    cat > conf/user.config << EOF
/*
 * User-specific configuration for Puhti
 * Generated by setup script on $(date)
 */

params {
    // SpecHLA settings - use container (recommended for Puhti)
    use_local_spechla = false
    spechla_path      = "/opt/SpecHLA"  // Path inside container

    // Container paths
    container_dir     = "${INSTALL_DIR}/hla_references/containers"

    // Database paths
    hlahd_db          = "${INSTALL_DIR}/hla_references/databases/hlahd_db"
    hlala_graphs      = "${INSTALL_DIR}/hla_references/databases/hlala_graphs"

    // Default output directory
    outdir            = "${INSTALL_DIR}/results"
}

// Puhti-specific process settings
process {
    executor = 'slurm'
    clusterOptions = '--account=${PROJECT_ID}'

    // Queue selection based on time
    queue = { task.time <= 2.h ? 'small' : task.time <= 3.d ? 'small' : 'longrun' }
}

// Singularity settings
singularity {
    enabled = true
    autoMounts = true
    runOptions = '--bind /scratch --bind /projappl'
}
EOF

    info "Pipeline configured with user settings"
    info "Configuration file: $INSTALL_DIR/hla_typing_pipeline/conf/user.config"
}

#-----------------------------------------------------------------------------
# Create Convenience Scripts
#-----------------------------------------------------------------------------
create_convenience_scripts() {
    log "Creating convenience scripts..."

    cd "$INSTALL_DIR"

    # Script to load environment
    cat > load_hla_env.sh << EOF
#!/bin/bash
# Load HLA typing pipeline environment on Puhti
# Usage: source load_hla_env.sh

# Load modules
module purge
module load gcc/11.3.0
module load samtools/1.15
module load bwa/0.7.17
module load nextflow/22.10.6
module load singularity

# Set paths
export HLA_PIPELINE_DIR="${INSTALL_DIR}/hla_typing_pipeline"
export HLA_CONTAINERS="${INSTALL_DIR}/hla_references/containers"
export PATH="\${HLA_PIPELINE_DIR}/bin:\${PATH}"

# Convenience aliases
alias run_hla="nextflow run \${HLA_PIPELINE_DIR}/main.nf -c \${HLA_PIPELINE_DIR}/conf/user.config"

echo "HLA typing pipeline environment loaded"
echo "Run 'run_hla --help' for usage information"
EOF
    chmod +x load_hla_env.sh

    # Quick run script
    cat > run_hla_typing.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_typing
#SBATCH --output=logs/hla_%j.out
#SBATCH --error=logs/hla_%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Load environment
source "$(dirname "$0")/load_hla_env.sh"

# Parse arguments
INPUT_BAM=""
INPUT_FASTQ_1=""
INPUT_FASTQ_2=""
INPUT_SAMPLESHEET=""
OUTDIR="./results"
TOOLS="spechla,hlahd"

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam) INPUT_BAM="$2"; shift 2 ;;
        --fastq1) INPUT_FASTQ_1="$2"; shift 2 ;;
        --fastq2) INPUT_FASTQ_2="$2"; shift 2 ;;
        --samplesheet) INPUT_SAMPLESHEET="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --tools) TOOLS="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Build command
CMD="nextflow run ${HLA_PIPELINE_DIR}/main.nf -c ${HLA_PIPELINE_DIR}/conf/user.config"
CMD+=" --outdir ${OUTDIR} --tools ${TOOLS}"

if [[ -n "$INPUT_BAM" ]]; then
    CMD+=" --input_bam ${INPUT_BAM}"
elif [[ -n "$INPUT_FASTQ_1" ]]; then
    CMD+=" --input_fastq_1 ${INPUT_FASTQ_1} --input_fastq_2 ${INPUT_FASTQ_2}"
elif [[ -n "$INPUT_SAMPLESHEET" ]]; then
    CMD+=" --input_samplesheet ${INPUT_SAMPLESHEET}"
else
    echo "Error: Please provide input (--bam, --fastq1/--fastq2, or --samplesheet)"
    exit 1
fi

echo "Running: $CMD"
eval $CMD
EOF
    chmod +x run_hla_typing.sh

    info "Convenience scripts created"
}

#-----------------------------------------------------------------------------
# Verify Installation
#-----------------------------------------------------------------------------
verify_installation() {
    log "Verifying installation..."

    echo ""
    echo "=============================================="
    echo "Installation Summary"
    echo "=============================================="

    # Check pipeline
    if [[ -f "$INSTALL_DIR/hla_typing_pipeline/main.nf" ]]; then
        echo -e "${GREEN}[OK]${NC} Pipeline installed"
    else
        echo -e "${RED}[FAIL]${NC} Pipeline not found"
    fi

    # Check SpecHLA container
    if [[ -f "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" ]]; then
        echo -e "${GREEN}[OK]${NC} SpecHLA container available"
    else
        echo -e "${RED}[FAIL]${NC} SpecHLA container not found"
        echo "       Copy container to: $INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif"
    fi

    # Check other containers
    for tool in hlahd arcashla optitype hlala xhla; do
        if [[ -f "$INSTALL_DIR/hla_references/containers/${tool}.sif" ]]; then
            echo -e "${GREEN}[OK]${NC} ${tool} container available"
        else
            echo -e "${YELLOW}[WARN]${NC} ${tool} container not available"
        fi
    done

    # Check HLA*LA reference graphs
    if [[ -d "$INSTALL_DIR/hla_references/databases/hlala_graphs/PRG_MHC_GRCh38_withIMGT" ]]; then
        echo -e "${GREEN}[OK]${NC} HLA*LA reference graphs available"
    else
        echo -e "${YELLOW}[WARN]${NC} HLA*LA reference graphs not found (required for HLA*LA)"
    fi

    # Check configuration
    if [[ -f "$INSTALL_DIR/hla_typing_pipeline/conf/user.config" ]]; then
        echo -e "${GREEN}[OK]${NC} User configuration created"
    else
        echo -e "${YELLOW}[WARN]${NC} User configuration not created"
    fi

    echo ""
    echo "=============================================="
    echo "Next Steps"
    echo "=============================================="
    echo "1. Load environment: source $INSTALL_DIR/load_hla_env.sh"
    echo "2. Test with a sample:"
    echo "   run_hla --input_bam /path/to/sample.bam --outdir ./test_results"
    echo ""
    echo "3. For batch processing, use SLURM submission:"
    echo "   sbatch run_hla_typing.sh --bam /path/to/sample.bam"
    echo ""
    echo "4. If some tools are missing, run their setup scripts:"
    echo "   - HLA-HD: $INSTALL_DIR/hla_references/setup_hlahd.sh"
    echo ""
    echo "Log file: $INSTALL_DIR/$LOG_FILE"
    echo "=============================================="
}

#-----------------------------------------------------------------------------
# Main
#-----------------------------------------------------------------------------
main() {
    echo "=============================================="
    echo "HLA Typing Pipeline - Full Setup for Puhti"
    echo "=============================================="
    echo ""

    parse_args "$@"
    detect_environment

    # Move log file to install directory once created
    touch "$LOG_FILE"

    load_modules
    create_directories

    # Move log to install dir
    mv "$LOG_FILE" "$INSTALL_DIR/$LOG_FILE"
    LOG_FILE="$INSTALL_DIR/$LOG_FILE"

    clone_pipeline
    install_spechla
    install_hlahd
    install_arcashla
    install_optitype
    install_hlala
    install_xhla
    install_python_deps
    configure_pipeline
    create_convenience_scripts
    verify_installation

    log "Setup complete!"
}

# Run main
main "$@"
