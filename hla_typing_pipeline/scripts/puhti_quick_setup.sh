#!/bin/bash
#=============================================================================
# HLA Typing Pipeline - Quick Setup Script for CSC Puhti
#=============================================================================
# Minimal setup using SpecHLA container (recommended)
#
# Usage:
#   ./puhti_quick_setup.sh [PROJECT_ID]
#
# Example:
#   ./puhti_quick_setup.sh project_2001234
#
#=============================================================================

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

#-----------------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------------

# Get project ID from argument or auto-detect
if [[ -n "$1" ]]; then
    PROJECT_ID="$1"
elif [[ "$PWD" =~ project_[0-9]+ ]]; then
    PROJECT_ID=$(echo "$PWD" | grep -o 'project_[0-9]*')
    log "Auto-detected project: $PROJECT_ID"
else
    error "Please provide project ID: $0 project_XXXXXXX"
fi

INSTALL_DIR="/scratch/${PROJECT_ID}/${USER}/hla_analysis"

#-----------------------------------------------------------------------------
# Main Setup
#-----------------------------------------------------------------------------

log "Setting up HLA typing pipeline in: $INSTALL_DIR"

# Load modules
log "Loading modules..."
module purge 2>/dev/null || true
module load gcc/11.3.0 2>/dev/null || module load gcc
module load samtools 2>/dev/null || true
module load bwa 2>/dev/null || true
module load nextflow 2>/dev/null || true
module load singularity 2>/dev/null || true

# Create directories
log "Creating directories..."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"
mkdir -p hla_typing_pipeline hla_references/containers logs results

# Clone pipeline
log "Cloning pipeline repository..."
if [[ -d "hla_typing_pipeline/.git" ]]; then
    cd hla_typing_pipeline && git pull && cd ..
else
    rm -rf hla_typing_pipeline
    git clone https://github.com/uoozcan/local_pipeline.git hla_typing_pipeline
fi

# Pull SpecHLA container
log "Setting up SpecHLA container..."
if [[ -f "hla_references/containers/spechla_1.0.7-3.sif" ]]; then
    log "SpecHLA container already exists"
else
    log "Attempting to pull SpecHLA container..."
    singularity pull "hla_references/containers/spechla_1.0.7-3.sif" \
        docker://quay.io/biocontainers/spechla:1.0.7--hdfd78af_3 2>/dev/null || \
    singularity pull "hla_references/containers/spechla_1.0.7-3.sif" \
        docker://deepomicslab/spechla:latest 2>/dev/null || {
        warn "Could not pull SpecHLA container automatically"
        warn "Please copy the container manually to: $INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif"
    }
fi

cd "$INSTALL_DIR"

# Create user config
log "Creating configuration..."
cat > hla_typing_pipeline/conf/user.config << EOF
/*
 * User configuration for Puhti - Generated $(date)
 * Uses container-based SpecHLA (recommended)
 */
params {
    // SpecHLA via container
    use_local_spechla = false
    spechla_path      = "/opt/SpecHLA"  // Path inside container
    container_dir     = "${INSTALL_DIR}/hla_references/containers"
    outdir            = "${INSTALL_DIR}/results"
}

process {
    executor = 'slurm'
    clusterOptions = '--account=${PROJECT_ID}'
}

singularity {
    enabled = true
    autoMounts = true
    runOptions = '--bind /scratch --bind /projappl'
}
EOF

# Create environment loader
cat > load_env.sh << EOF
#!/bin/bash
module purge
module load gcc/11.3.0 samtools bwa nextflow singularity
export HLA_DIR="${INSTALL_DIR}"
export HLA_CONTAINERS="${INSTALL_DIR}/hla_references/containers"
alias run_hla="nextflow run \${HLA_DIR}/hla_typing_pipeline/main.nf -c \${HLA_DIR}/hla_typing_pipeline/conf/user.config -profile singularity"
echo "Environment loaded. Use 'run_hla --help' for options."
EOF
chmod +x load_env.sh

# Create simple run script
cat > run_sample.sh << 'RUNSCRIPT'
#!/bin/bash
#SBATCH --job-name=hla_typing
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/hla_%j.out
#SBATCH --error=logs/hla_%j.err

source "$(dirname "$0")/load_env.sh"

if [[ -z "$1" ]]; then
    echo "Usage: sbatch run_sample.sh <input.bam> [output_dir]"
    exit 1
fi

INPUT_BAM="$1"
OUTDIR="${2:-./results}"
SAMPLE=$(basename "$INPUT_BAM" .bam)

nextflow run ${HLA_DIR}/hla_typing_pipeline/main.nf \
    -c ${HLA_DIR}/hla_typing_pipeline/conf/user.config \
    --input_bam "$INPUT_BAM" \
    --outdir "$OUTDIR" \
    --tools spechla \
    -profile singularity \
    -resume
RUNSCRIPT
chmod +x run_sample.sh

#-----------------------------------------------------------------------------
# Verification
#-----------------------------------------------------------------------------

echo ""
echo "=============================================="
echo "Setup Complete!"
echo "=============================================="

# Check installation
if [[ -f "$INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif" ]]; then
    echo -e "${GREEN}[OK]${NC} SpecHLA container available"
else
    echo -e "${RED}[FAIL]${NC} SpecHLA container not found"
    echo "       Please copy container to: $INSTALL_DIR/hla_references/containers/spechla_1.0.7-3.sif"
fi

if [[ -f "$INSTALL_DIR/hla_typing_pipeline/main.nf" ]]; then
    echo -e "${GREEN}[OK]${NC} Pipeline installed"
else
    echo -e "${RED}[FAIL]${NC} Pipeline installation incomplete"
fi

echo ""
echo "Quick Start:"
echo "  1. cd $INSTALL_DIR"
echo "  2. source load_env.sh"
echo "  3. run_hla --input_bam /path/to/sample.bam --outdir ./results"
echo ""
echo "For SLURM submission:"
echo "  sbatch run_sample.sh /path/to/sample.bam"
echo ""
echo "=============================================="
