#!/bin/bash
#
# HLA Typing Pipeline - Puhti Setup Script
# This script sets up the environment for running the HLA typing pipeline on CSC Puhti
#

set -e

echo "=============================================="
echo "HLA Typing Pipeline - Puhti Setup"
echo "=============================================="

# Detect project ID from current path or environment
if [[ "$PWD" =~ /scratch/(project_[0-9]+)/ ]]; then
    PROJECT_ID="${BASH_REMATCH[1]}"
elif [[ "$PWD" =~ /projappl/(project_[0-9]+)/ ]]; then
    PROJECT_ID="${BASH_REMATCH[1]}"
elif [[ -n "${CSC_PROJECT:-}" ]]; then
    PROJECT_ID="$CSC_PROJECT"
else
    echo "ERROR: Could not detect CSC project ID."
    echo ""
    echo "Please either:"
    echo "  1. Run this script from /scratch/project_XXXXXXX/ directory"
    echo "  2. Set CSC_PROJECT environment variable:"
    echo "     export CSC_PROJECT=project_XXXXXXX"
    echo "     ./puhti_setup.sh"
    exit 1
fi

echo "Detected project: $PROJECT_ID"

# Get username
USERNAME="${USER}"
echo "Username: $USERNAME"

# Set base directories
BASE_DIR="/scratch/${PROJECT_ID}/${USERNAME}/hla_analysis"
PIPELINE_DIR="${BASE_DIR}/hla_typing_pipeline"
SINGULARITY_CACHE="${BASE_DIR}/singularity_cache"
WORK_DIR="${BASE_DIR}/work"
RESULTS_DIR="${BASE_DIR}/results"
INPUT_DIR="${BASE_DIR}/input"
LOGS_DIR="${BASE_DIR}/logs"

echo ""
echo "Setting up directories..."
echo "Base directory: $BASE_DIR"

# Create directory structure
mkdir -p "$SINGULARITY_CACHE"
mkdir -p "$WORK_DIR"
mkdir -p "$RESULTS_DIR"
mkdir -p "$INPUT_DIR/bam"
mkdir -p "$INPUT_DIR/fastq"
mkdir -p "$LOGS_DIR"

echo "  - Created $SINGULARITY_CACHE"
echo "  - Created $WORK_DIR"
echo "  - Created $RESULTS_DIR"
echo "  - Created $INPUT_DIR"
echo "  - Created $LOGS_DIR"

# Create Puhti-specific Nextflow config
CONFIG_FILE="${PIPELINE_DIR}/conf/puhti.config"
mkdir -p "$(dirname "$CONFIG_FILE")"

cat > "$CONFIG_FILE" << EOF
/*
 * Puhti-specific configuration for HLA Typing Pipeline
 * CSC - IT Center for Science, Finland
 */

params {
    // Puhti-optimized resource settings
    max_cpus   = 40
    max_memory = '180.GB'
    max_time   = '72.h'

    // Output directory
    outdir = '${RESULTS_DIR}'

    // Container directory
    container_dir = '${SINGULARITY_CACHE}'
}

// Singularity settings for Puhti
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '${SINGULARITY_CACHE}'
    runOptions = '--bind /scratch,/projappl,/users'
}

// Process settings optimized for Puhti
process {
    executor = 'slurm'
    queue = 'small'
    clusterOptions = '--account=${PROJECT_ID}'

    // Default resources
    cpus = 4
    memory = '16.GB'
    time = '4.h'

    // Process-specific settings
    withName: 'SPECHLA|SPECHLA_FASTQ' {
        cpus = 20
        memory = '64.GB'
        time = '8.h'
    }

    withName: 'HLAHD|HLAHD_FASTQ' {
        cpus = 20
        memory = '64.GB'
        time = '8.h'
    }

    withName: 'HLALA' {
        cpus = 20
        memory = '80.GB'
        time = '12.h'
    }

    withName: 'ARCASHLA|ARCASHLA_FASTQ' {
        cpus = 10
        memory = '32.GB'
        time = '4.h'
    }

    withName: 'OPTITYPE|OPTITYPE_FASTQ' {
        cpus = 8
        memory = '32.GB'
        time = '4.h'
    }

    withName: 'QC_BAM|QC_FASTQ|FASTQC_BAM|FASTQC_FASTQ' {
        cpus = 4
        memory = '16.GB'
        time = '2.h'
    }

    withName: 'CONSENSUS|HLA_VISUALIZE|HLA_SUMMARY_REPORT|MULTIQC' {
        cpus = 2
        memory = '8.GB'
        time = '1.h'
    }
}

// Work directory
workDir = '${WORK_DIR}'

// Timeline and reports
def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
timeline {
    enabled = true
    file = '${RESULTS_DIR}/pipeline_info/timeline_\${timestamp}.html'
}
report {
    enabled = true
    file = '${RESULTS_DIR}/pipeline_info/report_\${timestamp}.html'
}
trace {
    enabled = true
    file = '${RESULTS_DIR}/pipeline_info/trace_\${timestamp}.txt'
}
EOF

echo "  - Created Puhti configuration: $CONFIG_FILE"

# Create environment setup script
ENV_FILE="${PIPELINE_DIR}/scripts/load_modules.sh"
cat > "$ENV_FILE" << 'EOF'
#!/bin/bash
# Load required modules for HLA typing pipeline on Puhti

# Clean module environment
module purge

# Load Nextflow
module load nextflow/23.10.0 2>/dev/null || module load nextflow

# Load Singularity (usually available by default)
module load singularity 2>/dev/null || true

# Load samtools for BAM processing
module load samtools 2>/dev/null || module load biokit

# Set Singularity cache directory
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$PWD/singularity_cache}"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"

# Disable Nextflow telemetry
export NXF_ANSI_LOG=false

echo "Modules loaded successfully"
echo "  - Nextflow: $(nextflow -version 2>&1 | head -1)"
echo "  - Singularity cache: $SINGULARITY_CACHEDIR"
EOF

chmod +x "$ENV_FILE"
echo "  - Created module loader: $ENV_FILE"

# Create a quick reference card
QUICKREF="${PIPELINE_DIR}/docs/QUICK_REFERENCE.txt"
mkdir -p "$(dirname "$QUICKREF")"
cat > "$QUICKREF" << EOF
=====================================================
HLA Typing Pipeline - Quick Reference for Puhti
=====================================================

PROJECT: ${PROJECT_ID}
USER: ${USERNAME}
BASE: ${BASE_DIR}

-----------------------------------------------------
QUICK COMMANDS
-----------------------------------------------------

# Load modules
source ${PIPELINE_DIR}/scripts/load_modules.sh

# Single BAM analysis
sbatch ${PIPELINE_DIR}/scripts/submit_hla_single.sh sample.bam

# Batch analysis
sbatch ${PIPELINE_DIR}/scripts/submit_hla_batch.sh samples.csv

# Interactive run (testing)
sinteractive -A ${PROJECT_ID} -t 04:00:00 -m 32G -c 8
source ${PIPELINE_DIR}/scripts/load_modules.sh
nextflow run ${PIPELINE_DIR}/main.nf --input_bam sample.bam -profile singularity

-----------------------------------------------------
DIRECTORY STRUCTURE
-----------------------------------------------------

${BASE_DIR}/
├── input/
│   ├── bam/        <- Place BAM files here
│   └── fastq/      <- Place FASTQ files here
├── results/        <- Output will be here
├── work/           <- Nextflow work directory
├── logs/           <- SLURM logs
└── singularity_cache/

-----------------------------------------------------
USEFUL COMMANDS
-----------------------------------------------------

# Check job status
squeue -u \$USER

# Cancel job
scancel <jobid>

# Check billing units
csc-projects

# View results
ls results/*/
cat results/summary/hla_summary_statistics.tsv

# Clean up after completion
rm -rf work/

=====================================================
EOF

echo "  - Created quick reference: $QUICKREF"

# Verify modules are available
echo ""
echo "Verifying module availability..."
if module avail nextflow 2>&1 | grep -q nextflow; then
    echo "  - Nextflow: Available"
else
    echo "  - Nextflow: Not found (may need to load manually)"
fi

if command -v singularity &> /dev/null; then
    echo "  - Singularity: Available ($(singularity --version 2>/dev/null || echo 'version unknown'))"
else
    echo "  - Singularity: Not in PATH (load module if needed)"
fi

# Final summary
echo ""
echo "=============================================="
echo "Setup Complete!"
echo "=============================================="
echo ""
echo "Next steps:"
echo ""
echo "1. Place your input files:"
echo "   - BAM files: ${INPUT_DIR}/bam/"
echo "   - FASTQ files: ${INPUT_DIR}/fastq/"
echo ""
echo "2. Load modules before running:"
echo "   source ${PIPELINE_DIR}/scripts/load_modules.sh"
echo ""
echo "3. Run the pipeline:"
echo "   # Single sample"
echo "   sbatch ${PIPELINE_DIR}/scripts/submit_hla_single.sh input/bam/sample.bam"
echo ""
echo "   # Multiple samples"
echo "   sbatch ${PIPELINE_DIR}/scripts/submit_hla_batch.sh samples.csv"
echo ""
echo "4. Check results:"
echo "   ls ${RESULTS_DIR}/"
echo ""
echo "For detailed instructions, see:"
echo "   ${PIPELINE_DIR}/docs/PUHTI_GUIDE.md"
echo ""
