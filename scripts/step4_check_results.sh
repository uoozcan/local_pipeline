#!/bin/bash
#SBATCH --job-name=analyze_results
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/analyze_results_%j.log

# STEP 4: Analyze Batch HLA Typing Results
# Generate comprehensive summary and quality reports

set -e

echo "=========================================="
echo "Batch HLA Typing Results Analysis"
echo "=========================================="
echo "Started: $(date)"
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

# Verify results directory exists
if [ ! -d "${RESULTS_DIR}" ]; then
    echo "ERROR: Results directory not found: ${RESULTS_DIR}"
    exit 1
fi

echo "Results directory: ${RESULTS_DIR}"
echo ""

# Load Python for analysis
module load python-data/3.10-23.03

# Count samples
TOTAL_SAMPLES=$(tail -n +2 ${SAMPLE_SHEET} | wc -l)
echo "Total samples in batch: ${TOTAL_SAMPLES}"
echo ""

# Check consensus results
echo "=========================================="
echo "Consensus Results"
echo "=========================================="
echo ""

CONSENSUS_FILE="${RESULTS_DIR}/majority_voting/all_samples.consensus.tsv"
if [ -f "${CONSENSUS_FILE}" ]; then
    echo "✓ Consensus file found"
    
    # Count samples with results
    SAMPLES_WITH_RESULTS=$(tail -n +2 ${CONSENSUS_FILE} | cut -f1 | sort -u | wc -l)
    TOTAL_CALLS=$(tail -n +2 ${CONSENSUS_FILE} | wc -l)
    
    echo "  Samples with HLA calls: ${SAMPLES_WITH_RESULTS}"
    echo "  Total HLA calls: ${TOTAL_CALLS}"
    echo ""
    
    # Count calls by gene
    echo "  Calls by gene:"
    for gene in A B C DRB1 DQB1 DPB1; do
        gene_count=$(tail -n +2 ${CONSENSUS_FILE} | grep -c "HLA-${gene}" || echo "0")
        echo "    HLA-${gene}: ${gene_count}"
    done
    echo ""
    
    # Show first 10 samples
    echo "  Preview (first 10 samples):"
    head -11 ${CONSENSUS_FILE} | column -t -s $'\t'
    echo ""
else
    echo "✗ Consensus file not found"
    echo ""
fi

# Check individual tool results
echo "=========================================="
echo "Individual Tool Results"
echo "=========================================="
echo ""

TOOL_SUCCESS=0
TOOL_FAILED=0

for sample_dir in ${RESULTS_DIR}/*/; do
    sample_id=$(basename ${sample_dir})
    
    # Skip non-sample directories
    [ "$sample_id" = "majority_voting" ] && continue
    [ "$sample_id" = "multiqc" ] && continue
    
    has_results=false
    
    # Check each tool
    for tool in arcashla xhla hlahd; do
        if [ -d "${sample_dir}/${tool}" ] && [ "$(ls -A ${sample_dir}/${tool} 2>/dev/null)" ]; then
            has_results=true
            break
        fi
    done
    
    if [ "$has_results" = true ]; then
        ((TOOL_SUCCESS++))
    else
        ((TOOL_FAILED++))
        echo "  ⚠ No results for: ${sample_id}"
    fi
done

echo "Samples with tool results: ${TOOL_SUCCESS}"
echo "Samples without results: ${TOOL_FAILED}"
echo ""

# Generate summary report
SUMMARY_REPORT="${RESULTS_DIR}/BATCH_ANALYSIS_SUMMARY.txt"

cat > ${SUMMARY_REPORT} << EOF
========================================
BATCH HLA TYPING ANALYSIS SUMMARY
========================================

Batch: ${BATCH_NAME}
Analysis Date: $(date)
Data Type: ${DATA_TYPE}

SAMPLE STATISTICS
-----------------
Total samples: ${TOTAL_SAMPLES}
Samples with results: ${SAMPLES_WITH_RESULTS}
Samples with tool output: ${TOOL_SUCCESS}
Samples failed: ${TOOL_FAILED}

TOOL CONFIGURATION
------------------
Tools: ${TOOLS}
HLA genes: ${HLA_GENES}
Sequence type: ${SEQ_TYPE}
Consensus: ${ENABLE_CONSENSUS}
Min tools for consensus: ${MIN_TOOLS_CONSENSUS}

RESULTS LOCATION
----------------
Base directory: ${RESULTS_DIR}
Consensus file: ${CONSENSUS_FILE}
Pipeline report: ${RESULTS_DIR}/pipeline_report.html

KEY FILES
---------
- All samples consensus: ${CONSENSUS_FILE}
- Pipeline execution report: ${RESULTS_DIR}/pipeline_report.html
- Timeline visualization: ${RESULTS_DIR}/pipeline_timeline.html
- Concordance analysis: ${RESULTS_DIR}/majority_voting/concordance_analysis/

DOWNLOAD COMMANDS
-----------------
# Download all results
scp -r ${USER}@puhti.csc.fi:${RESULTS_DIR} .

# Download only consensus
scp ${USER}@puhti.csc.fi:${CONSENSUS_FILE} .

# Download summary
scp ${USER}@puhti.csc.fi:${SUMMARY_REPORT} .

SAMPLE LIST
-----------
EOF

# Add sample list
tail -n +2 ${SAMPLE_SHEET} | cut -d',' -f1 >> ${SUMMARY_REPORT}

echo "=========================================="
echo "Summary report created: ${SUMMARY_REPORT}"
echo "=========================================="
echo ""

cat ${SUMMARY_REPORT}

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Completed: $(date)"
echo ""
