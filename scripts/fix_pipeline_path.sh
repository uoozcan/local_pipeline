#!/bin/bash
# Fix Pipeline Path in Configuration

if [ $# -eq 0 ]; then
    echo "Usage: bash fix_pipeline_path.sh <pipeline_directory>"
    echo ""
    echo "Example:"
    echo "  bash fix_pipeline_path.sh /scratch/project_2008084/ozcanumu/hla_typing_pipeline"
    echo ""
    echo "Or run diagnostics first:"
    echo "  bash diagnose_pipeline_path.sh"
    exit 1
fi

NEW_PIPELINE_DIR="$1"

echo "=========================================="
echo "Fix Pipeline Path in Configuration"
echo "=========================================="
echo ""

# Verify the new path
if [ ! -f "${NEW_PIPELINE_DIR}/main.nf" ]; then
    echo "✗ ERROR: main.nf not found at: ${NEW_PIPELINE_DIR}"
    echo ""
    echo "Please verify the path and try again."
    exit 1
fi

echo "✓ Verified pipeline exists at: ${NEW_PIPELINE_DIR}"
echo ""

# Load configuration to get current settings
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
CONFIG_FILE="${WORK_DIR}/scripts/config_bam_batch.sh"

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "✗ ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

echo "Configuration file: ${CONFIG_FILE}"
echo ""

# Backup the config file
BACKUP_FILE="${CONFIG_FILE}.backup_$(date +%Y%m%d_%H%M%S)"
cp "${CONFIG_FILE}" "${BACKUP_FILE}"
echo "✓ Backup created: ${BACKUP_FILE}"
echo ""

# Get current pipeline dir
source "${CONFIG_FILE}"
CURRENT_PIPELINE_DIR="${PIPELINE_DIR}"

echo "Current PIPELINE_DIR: ${CURRENT_PIPELINE_DIR}"
echo "New PIPELINE_DIR:     ${NEW_PIPELINE_DIR}"
echo ""

# Update the config file
if grep -q "^export PIPELINE_DIR=" "${CONFIG_FILE}"; then
    # Use proper escaping for sed
    ESCAPED_NEW=$(echo "${NEW_PIPELINE_DIR}" | sed 's/[\/&]/\\&/g')
    sed -i "s|^export PIPELINE_DIR=.*|export PIPELINE_DIR=\"${NEW_PIPELINE_DIR}\"|" "${CONFIG_FILE}"
    echo "✓ Updated PIPELINE_DIR in config file"
else
    echo "✗ Could not find PIPELINE_DIR line in config file"
    echo "  You may need to manually edit: ${CONFIG_FILE}"
    exit 1
fi

echo ""

# Verify the change
source "${CONFIG_FILE}"
if [ "${PIPELINE_DIR}" = "${NEW_PIPELINE_DIR}" ]; then
    echo "=========================================="
    echo "✓ Pipeline Path Successfully Updated!"
    echo "=========================================="
    echo ""
    echo "New configuration:"
    echo "  PIPELINE_DIR=${PIPELINE_DIR}"
    echo ""
    echo "Verification:"
    if [ -f "${PIPELINE_DIR}/main.nf" ]; then
        echo "  ✓ main.nf found"
    fi
    if [ -f "${PIPELINE_DIR}/nextflow.config" ]; then
        echo "  ✓ nextflow.config found"
    fi
    echo ""
    echo "You can now run the pipeline:"
    echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
    echo ""
    echo "If you need to restore the original config:"
    echo "  cp ${BACKUP_FILE} ${CONFIG_FILE}"
else
    echo "✗ Update may have failed. Please check manually."
    exit 1
fi
