#!/bin/bash
# Auto-Fix Script: Finds and configures HLA pipeline automatically

set -e

echo "=========================================="
echo "HLA Pipeline Auto-Fix Utility"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
CONFIG_FILE="${WORK_DIR}/scripts/config_bam_batch.sh"

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "✗ ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

source "${CONFIG_FILE}"

echo "Current configuration:"
echo "  PIPELINE_DIR: ${PIPELINE_DIR}"
echo ""

# Check if current config works
if [ -f "${PIPELINE_DIR}/main.nf" ]; then
    echo "✓ Pipeline is already correctly configured!"
    echo "  Location: ${PIPELINE_DIR}"
    echo ""
    echo "Pipeline appears to be working. If you're still having issues,"
    echo "they may be related to something else (modules, permissions, etc.)"
    exit 0
fi

echo "✗ Pipeline not found at configured location"
echo ""

# Search for pipeline
echo "Searching for HLA typing pipeline..."
echo ""

SEARCH_PATHS=(
    "/scratch/${PROJECT_ID}/hla_typing_pipeline"
    "/scratch/${PROJECT_ID}/ozcanumu/hla_typing_pipeline"
    "/scratch/${PROJECT_ID}/ozcanumu/HLA_typing_pipeline"
    "/scratch/${PROJECT_ID}/pipelines/hla_typing_pipeline"
    "${HOME}/hla_typing_pipeline"
    "${HOME}/HLA_typing_pipeline"
    "/projappl/${PROJECT_ID}/hla_typing_pipeline"
)

FOUND_PATH=""

for search_path in "${SEARCH_PATHS[@]}"; do
    if [ -f "${search_path}/main.nf" ]; then
        echo "✓ Found pipeline at: ${search_path}"
        FOUND_PATH="${search_path}"
        break
    fi
done

# If not found in common locations, do a broader search
if [ -z "${FOUND_PATH}" ]; then
    echo "Not found in common locations. Performing broader search..."
    echo "(This may take a minute...)"
    echo ""
    
    # Search in project space
    MAIN_NF_FILES=$(find /scratch/${PROJECT_ID} -name "main.nf" -type f 2>/dev/null | head -10)
    
    if [ ! -z "${MAIN_NF_FILES}" ]; then
        echo "Found main.nf files at:"
        echo "${MAIN_NF_FILES}"
        echo ""
        
        # Use the first one found
        FOUND_PATH=$(echo "${MAIN_NF_FILES}" | head -1 | xargs dirname)
        echo "Using: ${FOUND_PATH}"
        echo ""
        
        # Verify it looks like HLA pipeline
        if [ -f "${FOUND_PATH}/nextflow.config" ]; then
            echo "✓ Appears to be a valid Nextflow pipeline"
        else
            echo "⚠ WARNING: No nextflow.config found. This may not be the HLA pipeline."
            echo ""
            read -p "Use this path anyway? (y/n) " -n 1 -r
            echo ""
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "Aborting."
                exit 1
            fi
        fi
    fi
fi

if [ -z "${FOUND_PATH}" ]; then
    echo "=========================================="
    echo "✗ Pipeline Not Found"
    echo "=========================================="
    echo ""
    echo "Could not locate the HLA typing pipeline anywhere."
    echo ""
    echo "Options:"
    echo "  1. Clone/download the pipeline to Puhti"
    echo "  2. Check if you have access to the project space"
    echo "  3. Contact your supervisor or CSC support"
    echo ""
    echo "Expected pipeline structure:"
    echo "  pipeline_dir/"
    echo "  ├── main.nf"
    echo "  ├── nextflow.config"
    echo "  ├── modules/"
    echo "  └── workflows/"
    echo ""
    exit 1
fi

# Found the pipeline, now update config
echo "=========================================="
echo "Updating Configuration"
echo "=========================================="
echo ""

# Backup config
BACKUP_FILE="${CONFIG_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
cp "${CONFIG_FILE}" "${BACKUP_FILE}"
echo "✓ Backup created: ${BACKUP_FILE}"

# Update PIPELINE_DIR in config
if sed -i "s|^export PIPELINE_DIR=.*|export PIPELINE_DIR=\"${FOUND_PATH}\"|" "${CONFIG_FILE}"; then
    echo "✓ Configuration updated"
else
    echo "✗ Failed to update configuration"
    echo "  Please manually edit: ${CONFIG_FILE}"
    exit 1
fi

# Verify the change
source "${CONFIG_FILE}"

if [ "${PIPELINE_DIR}" = "${FOUND_PATH}" ] && [ -f "${PIPELINE_DIR}/main.nf" ]; then
    echo ""
    echo "=========================================="
    echo "✓ SUCCESS! Pipeline Configured"
    echo "=========================================="
    echo ""
    echo "Pipeline location: ${PIPELINE_DIR}"
    echo ""
    echo "Verification:"
    ls -lh "${PIPELINE_DIR}/main.nf"
    
    if [ -f "${PIPELINE_DIR}/nextflow.config" ]; then
        echo "✓ nextflow.config found"
    fi
    
    echo ""
    echo "You can now run the pipeline:"
    echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
    echo ""
    
    # Show what was changed
    echo "Changes made:"
    echo "  Old: export PIPELINE_DIR=\"${PIPELINE_DIR}\""
    echo "  New: export PIPELINE_DIR=\"${FOUND_PATH}\""
    echo ""
    echo "Backup saved at: ${BACKUP_FILE}"
    echo ""
else
    echo ""
    echo "✗ Configuration update may have failed"
    echo "  Please check: ${CONFIG_FILE}"
    echo ""
    exit 1
fi
