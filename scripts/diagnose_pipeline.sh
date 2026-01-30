#!/bin/bash
#
# Pipeline Diagnostic Script
# Purpose: Analyze pipeline structure and identify BAM input issues
#

echo "=========================================="
echo "HLA Pipeline Diagnostic Tool"
echo "=========================================="
echo ""

# Configuration
PIPELINE_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline"

if [ ! -d "${PIPELINE_DIR}" ]; then
    echo "ERROR: Pipeline directory not found: ${PIPELINE_DIR}"
    exit 1
fi

cd ${PIPELINE_DIR}

echo "Pipeline directory: ${PIPELINE_DIR}"
echo ""

# 1. Check main.nf line 20 (where error occurred)
echo "=========================================="
echo "1. Checking Line 20 of main.nf"
echo "=========================================="
echo ""
echo "Context around line 20:"
sed -n '15,25p' main.nf
echo ""

# 2. Check include statements
echo "=========================================="
echo "2. Checking ARCASHLA Include Statements"
echo "=========================================="
echo ""
grep -n "include.*ARCASHLA" main.nf
echo ""

# 3. Check available processes in arcashla module
echo "=========================================="
echo "3. Available ArcasHLA Processes"
echo "=========================================="
echo ""
if [ -f "modules/arcashla.nf" ]; then
    echo "Processes found in modules/arcashla.nf:"
    grep "^process ARCASHLA" modules/arcashla.nf
    echo ""
else
    echo "WARNING: modules/arcashla.nf not found"
    echo "Searching for arcashla module..."
    find . -name "*arcashla*" -type f
    echo ""
fi

# 4. Check workflow conditional logic
echo "=========================================="
echo "4. Workflow Input Type Logic"
echo "=========================================="
echo ""
echo "Checking for input_type conditionals:"
grep -n "input_type" main.nf | head -20
echo ""

# 5. Check workflow section
echo "=========================================="
echo "5. Workflow Structure"
echo "=========================================="
echo ""
echo "Main workflow section:"
grep -A 30 "^workflow {" main.nf | head -35
echo ""

# 6. Check for OptiType BAM process
echo "=========================================="
echo "6. OptiType Process Names"
echo "=========================================="
echo ""
if [ -f "modules/optitype.nf" ]; then
    grep "^process OPTITYPE" modules/optitype.nf
else
    echo "OptiType module not found"
fi
echo ""

# 7. Check for SpecHLA BAM process
echo "=========================================="
echo "7. SpecHLA Process Names"
echo "=========================================="
echo ""
if [ -f "modules/spechla.nf" ]; then
    grep "^process SPECHLA" modules/spechla.nf
else
    echo "SpecHLA module not found"
fi
echo ""

# 8. Summary
echo "=========================================="
echo "8. Diagnostic Summary"
echo "=========================================="
echo ""

# Count processes
ARCASHLA_COUNT=$(grep "^process ARCASHLA" modules/arcashla.nf 2>/dev/null | wc -l)
ARCASHLA_BAM_EXISTS=$(grep "^process ARCASHLA_BAM" modules/arcashla.nf 2>/dev/null | wc -l)

echo "ArcasHLA process count: ${ARCASHLA_COUNT}"
echo "ARCASHLA_BAM exists: ${ARCASHLA_BAM_EXISTS}"
echo ""

# Check include statement
INCLUDE_CHECK=$(grep "include.*ARCASHLA_BAM" main.nf 2>/dev/null | wc -l)
echo "ARCASHLA_BAM imported in main.nf: ${INCLUDE_CHECK}"
echo ""

# Check workflow usage
WORKFLOW_CHECK=$(grep "ARCASHLA_BAM" main.nf 2>/dev/null | wc -l)
echo "ARCASHLA_BAM used in workflow: ${WORKFLOW_CHECK}"
echo ""

echo "=========================================="
echo "9. Problem Identification"
echo "=========================================="
echo ""

if [ ${ARCASHLA_BAM_EXISTS} -gt 0 ]; then
    echo "✓ ARCASHLA_BAM process exists in modules/arcashla.nf"
    
    if [ ${INCLUDE_CHECK} -eq 0 ]; then
        echo "✗ PROBLEM: ARCASHLA_BAM is NOT imported in main.nf"
        echo "  Fix: Add include statement for ARCASHLA_BAM"
    else
        echo "✓ ARCASHLA_BAM is imported in main.nf"
    fi
    
    if [ ${WORKFLOW_CHECK} -eq 0 ]; then
        echo "✗ PROBLEM: ARCASHLA_BAM is NOT used in workflow"
        echo "  Fix: Update workflow to call ARCASHLA_BAM for BAM input"
    else
        echo "✓ ARCASHLA_BAM is used in workflow"
    fi
else
    echo "✗ CRITICAL: ARCASHLA_BAM process does not exist"
    echo "  This pipeline may not support BAM input for ArcasHLA"
fi

echo ""
echo "=========================================="
echo "10. Recommended Actions"
echo "=========================================="
echo ""

if [ ${ARCASHLA_BAM_EXISTS} -gt 0 ] && [ ${INCLUDE_CHECK} -eq 0 ]; then
    echo "ISSUE: Include statement missing"
    echo ""
    echo "The ARCASHLA_BAM process exists but is not imported."
    echo ""
    echo "FIX: Add this line to main.nf near other include statements:"
    echo "  include { ARCASHLA_BAM } from './modules/arcashla'"
    echo ""
    echo "Suggested fix location:"
    grep -n "include.*ARCASHLA[^_]" main.nf | head -1
    echo ""
    
elif [ ${ARCASHLA_BAM_EXISTS} -gt 0 ] && [ ${WORKFLOW_CHECK} -eq 0 ]; then
    echo "ISSUE: Workflow not using ARCASHLA_BAM"
    echo ""
    echo "The process is imported but not called in the workflow."
    echo ""
    echo "FIX: Update workflow to conditionally use ARCASHLA_BAM"
    echo "  Modify line ~20 to use ARCASHLA_BAM when input_type is 'bam'"
    echo ""
    
elif [ ${ARCASHLA_BAM_EXISTS} -eq 0 ]; then
    echo "ISSUE: Pipeline does not support BAM input for ArcasHLA"
    echo ""
    echo "Options:"
    echo "  1. Use a different pipeline version that supports BAM"
    echo "  2. Convert BAM to FASTQ and use FASTQ input"
    echo "  3. Modify the pipeline to add BAM support"
    echo ""
else
    echo "Multiple issues detected. Manual inspection needed."
    echo ""
    echo "Next steps:"
    echo "  1. Review the output above"
    echo "  2. Check line 20 of main.nf"
    echo "  3. Verify the include statements"
    echo "  4. Check workflow conditional logic"
fi

echo ""
echo "=========================================="
echo "Diagnostic Complete"
echo "=========================================="
echo ""
echo "Output saved. Share this with pipeline maintainers if needed."
echo ""
