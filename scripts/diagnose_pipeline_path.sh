#!/bin/bash
# Diagnostic Script: Find HLA Typing Pipeline Location

echo "=========================================="
echo "HLA Typing Pipeline - Location Diagnostic"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

echo "Current Configuration:"
echo "  PIPELINE_DIR: ${PIPELINE_DIR}"
echo "  Expected main.nf: ${PIPELINE_DIR}/main.nf"
echo ""

# Check if pipeline exists at configured location
echo "Checking configured location..."
if [ -f "${PIPELINE_DIR}/main.nf" ]; then
    echo "  ✓ Pipeline found at configured location"
    echo ""
    echo "Pipeline is correctly configured!"
    exit 0
else
    echo "  ✗ Pipeline NOT found at configured location"
fi

echo ""
echo "=========================================="
echo "Searching for Pipeline..."
echo "=========================================="
echo ""

# Search in common locations
SEARCH_PATHS=(
    "${BASE_DIR}/hla_typing_pipeline"
    "${BASE_DIR}/HLA_typing_pipeline"
    "${BASE_DIR}/pipeline"
    "${WORK_DIR}/hla_typing_pipeline"
    "${WORK_DIR}/HLA_typing_pipeline"
    "/scratch/${PROJECT_ID}/hla_typing_pipeline"
    "/scratch/${PROJECT_ID}/ozcanumu/hla_typing_pipeline"
    "${HOME}/hla_typing_pipeline"
    "${HOME}/HLA_typing_pipeline"
)

FOUND_LOCATIONS=()

for search_path in "${SEARCH_PATHS[@]}"; do
    if [ -f "${search_path}/main.nf" ]; then
        echo "✓ Found at: ${search_path}"
        FOUND_LOCATIONS+=("${search_path}")
    fi
done

echo ""

if [ ${#FOUND_LOCATIONS[@]} -eq 0 ]; then
    echo "=========================================="
    echo "Pipeline Not Found Anywhere!"
    echo "=========================================="
    echo ""
    echo "The HLA typing pipeline could not be located."
    echo ""
    echo "Options:"
    echo "  1. Clone/copy the pipeline to a location on Puhti"
    echo "  2. Check if pipeline is in a different directory"
    echo ""
    echo "Typical pipeline structure should include:"
    echo "  - main.nf"
    echo "  - nextflow.config"
    echo "  - modules/"
    echo "  - workflows/"
    echo ""
    echo "To search manually:"
    echo "  find /scratch/${PROJECT_ID} -name \"main.nf\" -type f 2>/dev/null"
    echo ""
    exit 1
else
    echo "=========================================="
    echo "Found ${#FOUND_LOCATIONS[@]} Pipeline Location(s)"
    echo "=========================================="
    echo ""
    
    # Show first found location
    CORRECT_PATH="${FOUND_LOCATIONS[0]}"
    echo "Recommended pipeline path:"
    echo "  ${CORRECT_PATH}"
    echo ""
    
    # Check pipeline version/contents
    if [ -f "${CORRECT_PATH}/nextflow.config" ]; then
        echo "Pipeline appears complete (has nextflow.config)"
    fi
    
    echo ""
    echo "=========================================="
    echo "How to Fix"
    echo "=========================================="
    echo ""
    echo "Update config_bam_batch.sh:"
    echo ""
    echo "  Change line:"
    echo "    export PIPELINE_DIR=\"${PIPELINE_DIR}\""
    echo ""
    echo "  To:"
    echo "    export PIPELINE_DIR=\"${CORRECT_PATH}\""
    echo ""
    echo "Or run the fix script:"
    echo "  bash ${SCRIPTS_DIR}/fix_pipeline_path.sh \"${CORRECT_PATH}\""
    echo ""
fi

# Additional search - find all main.nf files
echo "=========================================="
echo "All main.nf Files Found"
echo "=========================================="
echo ""
echo "Searching entire project space (this may take a moment)..."
find /scratch/${PROJECT_ID} -name "main.nf" -type f 2>/dev/null | while read nf_file; do
    nf_dir=$(dirname "${nf_file}")
    echo "  ${nf_dir}"
done

echo ""
echo "Diagnostic complete."
