#!/bin/bash
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=reorganize_pipeline

# Reorganize HLA typing pipeline for GitHub repository
# Usage: sbatch reorganize_for_github.sh

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Configuration
CURRENT_DIR="/scratch/project_2008084/hla_rnaseq_analysis"
NEW_REPO_DIR="/scratch/project_2008084/hla-typing-pipeline"
ARCHIVE_DIR="/scratch/project_2008084/hla_archive_$(date +%Y%m%d)"

log_info "Starting pipeline reorganization for GitHub..."
log_info "Current directory: $CURRENT_DIR"
log_info "New repository: $NEW_REPO_DIR"
log_info "Archive directory: $ARCHIVE_DIR"

# Create directories
log_info "Creating directory structure..."
mkdir -p "$NEW_REPO_DIR"/{modules,bin,conf,docs,examples,scripts/{setup,utils}}
mkdir -p "$ARCHIVE_DIR"/{backups,old_results,raw_data,logs}

# ============================================================================
# STEP 1: Copy core pipeline files
# ============================================================================
log_info "Copying core pipeline files..."

# Main pipeline files
cp "$CURRENT_DIR/main.nf" "$NEW_REPO_DIR/" 2>/dev/null || log_warn "main.nf not found in root"
cp "$CURRENT_DIR/nextflow.config" "$NEW_REPO_DIR/" 2>/dev/null || log_warn "nextflow.config not found"
cp "$CURRENT_DIR/LICENSE" "$NEW_REPO_DIR/" 2>/dev/null || log_warn "LICENSE not found"

# Check pipeline subdirectory
if [ -d "$CURRENT_DIR/pipeline" ]; then
    log_info "Found pipeline subdirectory, copying from there..."
    cp "$CURRENT_DIR/pipeline/main.nf" "$NEW_REPO_DIR/" 2>/dev/null || true
    cp "$CURRENT_DIR/pipeline/nextflow.config" "$NEW_REPO_DIR/" 2>/dev/null || true
fi

# Modules
log_info "Copying modules..."
if [ -d "$CURRENT_DIR/pipeline/modules" ]; then
    cp "$CURRENT_DIR/pipeline/modules"/*.nf "$NEW_REPO_DIR/modules/" 2>/dev/null || true
elif [ -d "$CURRENT_DIR/modules" ]; then
    cp "$CURRENT_DIR/modules"/*.nf "$NEW_REPO_DIR/modules/" 2>/dev/null || true
fi

# Bin scripts
log_info "Copying bin scripts..."
if [ -d "$CURRENT_DIR/pipeline/bin" ]; then
    cp "$CURRENT_DIR/pipeline/bin"/* "$NEW_REPO_DIR/bin/" 2>/dev/null || true
elif [ -d "$CURRENT_DIR/bin" ]; then
    cp "$CURRENT_DIR/bin"/* "$NEW_REPO_DIR/bin/" 2>/dev/null || true
fi

# Configuration files
log_info "Copying configuration files..."
if [ -d "$CURRENT_DIR/pipeline/conf" ]; then
    cp "$CURRENT_DIR/pipeline/conf"/*.config "$NEW_REPO_DIR/conf/" 2>/dev/null || true
elif [ -d "$CURRENT_DIR/conf" ]; then
    cp "$CURRENT_DIR/conf"/*.config "$NEW_REPO_DIR/conf/" 2>/dev/null || true
fi

# ============================================================================
# STEP 2: Copy and organize documentation
# ============================================================================
log_info "Organizing documentation..."

# Main README
if [ -f "$CURRENT_DIR/README.md" ]; then
    cp "$CURRENT_DIR/README.md" "$NEW_REPO_DIR/"
elif [ -f "$CURRENT_DIR/pipeline/README.md" ]; then
    cp "$CURRENT_DIR/pipeline/README.md" "$NEW_REPO_DIR/"
else
    log_warn "README.md not found, will need to create one"
fi

# CHANGELOG
if [ -f "$CURRENT_DIR/CHANGELOG.md" ]; then
    cp "$CURRENT_DIR/CHANGELOG.md" "$NEW_REPO_DIR/"
elif [ -f "$CURRENT_DIR/pipeline/CHANGELOG.md" ]; then
    cp "$CURRENT_DIR/pipeline/CHANGELOG.md" "$NEW_REPO_DIR/"
fi

# Other documentation
for doc_file in INSTALL_PUHTI.md QUICKSTART.md DEPLOYMENT.md DOWNLOAD_GUIDE.md RELEASE_SUMMARY.md; do
    if [ -f "$CURRENT_DIR/$doc_file" ]; then
        cp "$CURRENT_DIR/$doc_file" "$NEW_REPO_DIR/docs/"
    elif [ -f "$CURRENT_DIR/pipeline/$doc_file" ]; then
        cp "$CURRENT_DIR/pipeline/$doc_file" "$NEW_REPO_DIR/docs/"
    fi
done

# Copy docs directory if exists
if [ -d "$CURRENT_DIR/pipeline/docs" ]; then
    cp "$CURRENT_DIR/pipeline/docs"/* "$NEW_REPO_DIR/docs/" 2>/dev/null || true
elif [ -d "$CURRENT_DIR/docs" ]; then
    cp "$CURRENT_DIR/docs"/* "$NEW_REPO_DIR/docs/" 2>/dev/null || true
fi

# Rename INSTALL_PUHTI.md to INSTALLATION.md
if [ -f "$NEW_REPO_DIR/docs/INSTALL_PUHTI.md" ]; then
    mv "$NEW_REPO_DIR/docs/INSTALL_PUHTI.md" "$NEW_REPO_DIR/docs/INSTALLATION.md"
fi

# ============================================================================
# STEP 3: Create example files
# ============================================================================
log_info "Creating example files..."

# Copy example parameter file
if [ -f "$CURRENT_DIR/params.yaml.example" ]; then
    cp "$CURRENT_DIR/params.yaml.example" "$NEW_REPO_DIR/examples/"
elif [ -f "$CURRENT_DIR/pipeline/params.yaml" ]; then
    cp "$CURRENT_DIR/pipeline/params.yaml" "$NEW_REPO_DIR/examples/params.yaml.example"
fi

# Create example samplesheet
if [ -f "$CURRENT_DIR/samplesheet_batch1_VenEx_DNA.csv" ]; then
    head -n 3 "$CURRENT_DIR/samplesheet_batch1_VenEx_DNA.csv" > "$NEW_REPO_DIR/examples/samplesheet_example.csv"
fi

# ============================================================================
# STEP 4: Copy utility scripts
# ============================================================================
log_info "Organizing utility scripts..."

# Setup scripts
for script in verify_hla_references.sh; do
    if [ -f "$CURRENT_DIR/$script" ]; then
        cp "$CURRENT_DIR/$script" "$NEW_REPO_DIR/scripts/setup/"
    fi
done

# Utility scripts
for script in generate_samplesheet.sh; do
    if [ -f "$CURRENT_DIR/$script" ]; then
        cp "$CURRENT_DIR/$script" "$NEW_REPO_DIR/scripts/utils/"
    fi
done

# Copy selected scripts from scripts directory
if [ -d "$CURRENT_DIR/scripts" ]; then
    log_info "Copying utility scripts from scripts directory..."
    
    # Setup scripts
    for script in step1_setup_puhti_bam.sh; do
        if [ -f "$CURRENT_DIR/scripts/$script" ]; then
            cp "$CURRENT_DIR/scripts/$script" "$NEW_REPO_DIR/scripts/setup/"
        fi
    done
    
    # Utility scripts
    for script in validate_bam_for_spechla.sh check_bam.sh; do
        if [ -f "$CURRENT_DIR/scripts/$script" ]; then
            cp "$CURRENT_DIR/scripts/$script" "$NEW_REPO_DIR/scripts/utils/"
        fi
    done
fi

# ============================================================================
# STEP 5: Create .gitignore
# ============================================================================
log_info "Creating .gitignore..."

cat > "$NEW_REPO_DIR/.gitignore" << 'EOF'
# Nextflow work directories and logs
work/
.nextflow/
.nextflow.log*
*.log

# Results and outputs
results/
results_*/
output/
reports/
analysis_output_*/
analysis_with_mv_*/

# Raw data
raw_fastq/
raw_bam/
raw_data/
pipeline_input/
pipeline_input_bam/
test_input/
merged_fastq/

# Singularity cache
singularity_cache/
.singularity/

# Reference files (too large for git)
references/
*.fa
*.fa.gz
*.fai
*.fa.alt
*.bwt
*.pac
*.ann
*.amb
*.sa

# Backups
*_backup*/
pipeline_backup_*/
hla_pipeline_backup/

# Logs
logs/
*.out
*.err

# Personal files and keys
*.sec
*.key
*_key.*
.passphrase
.c4gh_passphrase
*.c4gh

# Temporary files
*.tmp
*.temp
*_tmp/
*.swp
*.swo
*~

# OS files
.DS_Store
Thumbs.db

# Python
__pycache__/
*.pyc
*.pyo
*.egg-info/
.pytest_cache/

# Archives
*.tar.gz
*.tar.bz2
*.zip
*.tgz

# CSC Allas specific
allas_file_list.txt

# Sample sheets with real data
samplesheet_batch*.csv
batch*_file_list*.txt

# Work-in-progress files
TODO.md
NOTES.md
scratch/

# Editor files
.vscode/
.idea/
*.sublime-*

# Test results
quick_test_*/
*_test_results*/
pipeline_diagnosis.txt
EOF

# ============================================================================
# STEP 6: Create README template (if needed)
# ============================================================================
if [ ! -f "$NEW_REPO_DIR/README.md" ]; then
    log_info "Creating README template..."
    
    cat > "$NEW_REPO_DIR/README.md" << 'EOF'
# HLA Typing Pipeline

A comprehensive Nextflow pipeline for HLA typing from NGS data using multiple specialized tools with confidence-weighted majority voting.

## Features

- **Multiple HLA typing tools**: OptiType, ArcasHLA, SpecHLA, xHLA, HLA-HD, HLA\*LA, Seq2HLA, HLAProfiler, Kourami
- **Consensus calling**: Confidence-weighted majority voting algorithm
- **Flexible input**: Supports BAM and FASTQ files from RNA-seq, WES, and WGS
- **HPC optimized**: SLURM integration for CSC Puhti
- **Containerized**: Singularity containers for reproducibility

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/hla-typing-pipeline.git
cd hla-typing-pipeline

# 2. See installation guide
cat docs/INSTALLATION.md

# 3. Run the pipeline
nextflow run main.nf -profile puhti --input samplesheet.csv
```

## Documentation

- [Installation Guide](docs/INSTALLATION.md)
- [Usage Guide](docs/USAGE.md)
- [Quick Start](docs/QUICKSTART.md)
- [Deployment Guide](docs/DEPLOYMENT.md)

## Citation

If you use this pipeline, please cite:
- The individual HLA typing tools used
- Nextflow: Di Tommaso, P., et al. (2017). Nat Biotechnol. 35, 316-319

## License

See [LICENSE](LICENSE) file for details.

## Contact

For questions and support, please open an issue on GitHub.
EOF
fi

# ============================================================================
# STEP 7: Archive old files
# ============================================================================
log_info "Archiving old files..."

# Backup directories
log_info "Moving backup directories..."
mv "$CURRENT_DIR"/pipeline_backup_* "$ARCHIVE_DIR/backups/" 2>/dev/null || true

# Old results
log_info "Moving old results..."
mv "$CURRENT_DIR"/results "$ARCHIVE_DIR/old_results/" 2>/dev/null || true
mv "$CURRENT_DIR"/results_bam "$ARCHIVE_DIR/old_results/" 2>/dev/null || true
mv "$CURRENT_DIR"/quick_test_results* "$ARCHIVE_DIR/old_results/" 2>/dev/null || true
mv "$CURRENT_DIR"/analysis_output_* "$ARCHIVE_DIR/old_results/" 2>/dev/null || true
mv "$CURRENT_DIR"/analysis_with_mv_* "$ARCHIVE_DIR/old_results/" 2>/dev/null || true

# Raw data (keep references to paths, move actual data)
log_info "Moving raw data..."
mv "$CURRENT_DIR"/raw_fastq "$ARCHIVE_DIR/raw_data/" 2>/dev/null || true
mv "$CURRENT_DIR"/raw_bam "$ARCHIVE_DIR/raw_data/" 2>/dev/null || true

# Logs
log_info "Moving logs..."
mv "$CURRENT_DIR"/logs "$ARCHIVE_DIR/" 2>/dev/null || true
mv "$CURRENT_DIR"/*.log "$ARCHIVE_DIR/logs/" 2>/dev/null || true
mv "$CURRENT_DIR"/*.out "$ARCHIVE_DIR/logs/" 2>/dev/null || true
mv "$CURRENT_DIR"/*.err "$ARCHIVE_DIR/logs/" 2>/dev/null || true

# Work directories
log_info "Moving work directories..."
mv "$CURRENT_DIR"/work "$ARCHIVE_DIR/" 2>/dev/null || true
mv "$CURRENT_DIR"/work_batch_bam "$ARCHIVE_DIR/" 2>/dev/null || true

# Nextflow cache
mv "$CURRENT_DIR"/.nextflow* "$ARCHIVE_DIR/" 2>/dev/null || true

# ============================================================================
# STEP 8: Create summary document
# ============================================================================
log_info "Creating reorganization summary..."

cat > "$NEW_REPO_DIR/REORGANIZATION_NOTES.md" << EOF
# Pipeline Reorganization Notes

**Date**: $(date)
**Original location**: $CURRENT_DIR
**New repository**: $NEW_REPO_DIR
**Archive location**: $ARCHIVE_DIR

## Structure Changes

### New Repository Structure
\`\`\`
hla-typing-pipeline/
├── main.nf
├── nextflow.config
├── modules/          # All Nextflow modules
├── bin/              # Pipeline scripts
├── conf/             # Configuration files
├── docs/             # Documentation
├── examples/         # Example files
└── scripts/          # Utility scripts
    ├── setup/
    └── utils/
\`\`\`

### Archived Content
- Backup directories → $ARCHIVE_DIR/backups/
- Old results → $ARCHIVE_DIR/old_results/
- Raw data → $ARCHIVE_DIR/raw_data/
- Logs → $ARCHIVE_DIR/logs/
- Work directories → $ARCHIVE_DIR/

### Important Paths

Reference data remains at:
- /scratch/project_2008084/hla_references/

Container images remain at:
- /scratch/project_2008084/hla_references/singularity_cache/containers/

## Next Steps

1. **Review the new repository structure**
   \`\`\`bash
   cd $NEW_REPO_DIR
   tree -L 2
   \`\`\`

2. **Initialize Git repository**
   \`\`\`bash
   cd $NEW_REPO_DIR
   git init
   git add .
   git commit -m "Initial commit: Organized HLA typing pipeline"
   \`\`\`

3. **Test the pipeline**
   \`\`\`bash
   # Run a test with example data
   nextflow run main.nf -profile puhti --input examples/samplesheet_example.csv
   \`\`\`

4. **Update documentation**
   - Review and update README.md
   - Update paths in documentation files
   - Add usage examples

5. **Create GitHub repository**
   \`\`\`bash
   # On GitHub, create a new repository
   # Then push:
   git remote add origin https://github.com/YOUR_USERNAME/hla-typing-pipeline.git
   git branch -M main
   git push -u origin main
   \`\`\`

## Files Not Copied (Excluded)

- Personal keys and secrets (.sec, .key files)
- Large data files (raw FASTQ/BAM files)
- Temporary work directories
- Old backup directories
- Test results
- Logs

These files are archived in: $ARCHIVE_DIR

## Configuration Updates Needed

Update these configuration files with correct paths:
- [ ] conf/puhti.config - Singularity cache path
- [ ] conf/params.config - Reference genome paths
- [ ] docs/INSTALLATION.md - Update installation paths

EOF

# ============================================================================
# STEP 9: Set permissions
# ============================================================================
log_info "Setting permissions..."
chmod -R u+rwX,go+rX "$NEW_REPO_DIR"
chmod +x "$NEW_REPO_DIR"/scripts/setup/*.sh 2>/dev/null || true
chmod +x "$NEW_REPO_DIR"/scripts/utils/*.sh 2>/dev/null || true

# ============================================================================
# Summary
# ============================================================================
log_info "Reorganization complete!"
echo ""
echo "=================================================="
echo "Summary:"
echo "=================================================="
echo "New repository: $NEW_REPO_DIR"
echo "Archive: $ARCHIVE_DIR"
echo ""
echo "Next steps:"
echo "1. Review: cd $NEW_REPO_DIR && tree -L 2"
echo "2. Read: cat $NEW_REPO_DIR/REORGANIZATION_NOTES.md"
echo "3. Initialize Git: cd $NEW_REPO_DIR && git init"
echo "4. Test pipeline: nextflow run main.nf -profile puhti --help"
echo ""
echo "See REORGANIZATION_NOTES.md for detailed information"
echo "=================================================="
