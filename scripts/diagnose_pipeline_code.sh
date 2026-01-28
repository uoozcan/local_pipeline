#!/bin/bash
# Diagnose Pipeline Code Issues

echo "=========================================="
echo "Pipeline Code Diagnostics"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

echo "Pipeline location: ${PIPELINE_DIR}"
echo ""

# Check pipeline structure
echo "=========================================="
echo "Pipeline Structure Check"
echo "=========================================="
echo ""

REQUIRED_FILES=(
    "main.nf"
    "nextflow.config"
    "workflows/hla_typing.nf"
    "modules/local/samplesheet_check.nf"
)

MISSING_FILES=0

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "${PIPELINE_DIR}/${file}" ]; then
        echo "✓ Found: ${file}"
    else
        echo "✗ Missing: ${file}"
        ((MISSING_FILES++))
    fi
done

echo ""

if [ ${MISSING_FILES} -gt 0 ]; then
    echo "⚠ Warning: ${MISSING_FILES} required file(s) missing"
    echo "  Pipeline may be incomplete"
fi

# Check for workflow files
echo ""
echo "=========================================="
echo "Workflow Files"
echo "=========================================="
echo ""

if [ -d "${PIPELINE_DIR}/workflows" ]; then
    echo "Workflow files:"
    ls -lh "${PIPELINE_DIR}/workflows/"*.nf 2>/dev/null || echo "  No .nf files found"
else
    echo "✗ workflows/ directory not found"
fi

# Check for module files
echo ""
echo "=========================================="
echo "Module Files"
echo "=========================================="
echo ""

if [ -d "${PIPELINE_DIR}/modules" ]; then
    echo "Module directories:"
    find "${PIPELINE_DIR}/modules" -name "*.nf" | head -10
    MODULE_COUNT=$(find "${PIPELINE_DIR}/modules" -name "*.nf" | wc -l)
    echo ""
    echo "Total modules: ${MODULE_COUNT}"
else
    echo "✗ modules/ directory not found"
fi

# Check Git status
echo ""
echo "=========================================="
echo "Pipeline Version Info"
echo "=========================================="
echo ""

cd "${PIPELINE_DIR}"

if [ -d ".git" ]; then
    echo "Git repository detected"
    echo ""
    
    echo "Current branch:"
    git branch --show-current 2>/dev/null || echo "  Unable to determine"
    
    echo ""
    echo "Latest commit:"
    git log -1 --oneline 2>/dev/null || echo "  Unable to retrieve"
    
    echo ""
    echo "Remote URL:"
    git remote get-url origin 2>/dev/null || echo "  No remote configured"
    
    echo ""
    echo "Available tags:"
    git tag -l 2>/dev/null | tail -5 || echo "  No tags found"
else
    echo "Not a Git repository"
fi

# Check for README
echo ""
echo "=========================================="
echo "Documentation Check"
echo "=========================================="
echo ""

if [ -f "${PIPELINE_DIR}/README.md" ]; then
    echo "✓ README.md found"
    echo ""
    echo "First 30 lines:"
    head -30 "${PIPELINE_DIR}/README.md"
else
    echo "✗ README.md not found"
fi

# Check nextflow.config
echo ""
echo "=========================================="
echo "Pipeline Configuration"
echo "=========================================="
echo ""

if [ -f "${PIPELINE_DIR}/nextflow.config" ]; then
    echo "✓ nextflow.config found"
    echo ""
    echo "Manifest section:"
    grep -A 10 "manifest {" "${PIPELINE_DIR}/nextflow.config" 2>/dev/null || echo "  Not found"
else
    echo "✗ nextflow.config not found"
fi

# Search for ch_input in pipeline code
echo ""
echo "=========================================="
echo "Channel 'ch_input' References"
echo "=========================================="
echo ""

echo "Searching for 'ch_input' in pipeline code..."
echo ""

CH_INPUT_REFS=$(grep -r "ch_input" "${PIPELINE_DIR}"/*.nf "${PIPELINE_DIR}/workflows/"*.nf 2>/dev/null | wc -l)

if [ ${CH_INPUT_REFS} -gt 0 ]; then
    echo "Found ${CH_INPUT_REFS} references to 'ch_input'"
    echo ""
    echo "Files containing 'ch_input':"
    grep -l "ch_input" "${PIPELINE_DIR}"/*.nf "${PIPELINE_DIR}/workflows/"*.nf 2>/dev/null
    echo ""
    echo "Sample occurrences:"
    grep -n "ch_input" "${PIPELINE_DIR}"/*.nf "${PIPELINE_DIR}/workflows/"*.nf 2>/dev/null | head -5
else
    echo "✗ No references to 'ch_input' found"
    echo "  This variable may be undefined or missing"
fi

# Check for similar channel names
echo ""
echo "Other channel definitions in main.nf:"
grep "^Channel" "${PIPELINE_DIR}/main.nf" 2>/dev/null || echo "  None found"
grep "ch_" "${PIPELINE_DIR}/main.nf" 2>/dev/null | head -10

echo ""
echo "=========================================="
echo "Diagnostic Summary"
echo "=========================================="
echo ""

echo "Issue: 'ch_input' channel not found in pipeline code"
echo ""
echo "Possible causes:"
echo "  1. Pipeline code is incomplete or corrupted"
echo "  2. Pipeline version mismatch"
echo "  3. Missing workflow files"
echo "  4. Incompatible pipeline branch/version"
echo ""
echo "Recommended actions:"
echo "  1. Check pipeline version/branch"
echo "  2. Ensure pipeline is complete (run git pull/fetch)"
echo "  3. Try using a stable release tag"
echo "  4. Check pipeline documentation for requirements"
echo ""
echo "For help:"
echo "  cat ${PIPELINE_DIR}/README.md"
echo "  cat ${PIPELINE_DIR}/docs/usage.md"
echo ""
