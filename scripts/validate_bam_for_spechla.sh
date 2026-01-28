#!/bin/bash
#
# SpecHLA BAM Compatibility Checker
# Purpose: Diagnose potential issues before running SpecHLA analysis
# Usage: bash validate_bam_for_spechla.sh /path/to/alignment.bam
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Symbols
CHECK="${GREEN}✓${NC}"
CROSS="${RED}✗${NC}"
WARN="${YELLOW}⚠${NC}"
INFO="${BLUE}ℹ${NC}"

echo "=============================================="
echo "  SpecHLA BAM Compatibility Checker"
echo "=============================================="
echo ""

# Check if BAM file provided
if [ $# -eq 0 ]; then
    echo "${CROSS} No BAM file provided"
    echo ""
    echo "Usage: bash $0 <bam_file>"
    echo "Example: bash $0 alignment-sorted.bam"
    exit 1
fi

BAM_FILE="$1"
SCORE=0
MAX_SCORE=0
CRITICAL_ISSUES=0
WARNINGS=0

# Function to print test result
print_result() {
    local status=$1
    local message=$2
    local details=$3
    local critical=$4
    
    MAX_SCORE=$((MAX_SCORE + 1))
    
    if [ "$status" = "pass" ]; then
        echo -e "${CHECK} ${message}"
        SCORE=$((SCORE + 1))
    elif [ "$status" = "fail" ]; then
        echo -e "${CROSS} ${message}"
        if [ "$critical" = "yes" ]; then
            CRITICAL_ISSUES=$((CRITICAL_ISSUES + 1))
        fi
    elif [ "$status" = "warn" ]; then
        echo -e "${WARN} ${message}"
        WARNINGS=$((WARNINGS + 1))
        SCORE=$((SCORE + 0))  # Half credit
    else
        echo -e "${INFO} ${message}"
    fi
    
    if [ -n "$details" ]; then
        echo "     $details"
    fi
    echo ""
}

echo "Analyzing: ${BAM_FILE}"
echo ""

# ============================================================================
# TEST 1: File Existence and Accessibility
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1. File Validation"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -f "${BAM_FILE}" ]; then
    print_result "fail" "BAM file not found" "Path: ${BAM_FILE}" "yes"
    echo "CRITICAL: Cannot proceed without valid BAM file"
    exit 1
else
    FILE_SIZE=$(du -h "${BAM_FILE}" | cut -f1)
    print_result "pass" "BAM file exists" "Size: ${FILE_SIZE}"
fi

# Check if file is readable
if [ ! -r "${BAM_FILE}" ]; then
    print_result "fail" "BAM file not readable" "Check file permissions" "yes"
    exit 1
else
    print_result "pass" "BAM file is readable"
fi

# ============================================================================
# TEST 2: BAM Index
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2. BAM Index"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "${BAM_FILE}.bai" ]; then
    print_result "pass" "BAM index (.bai) exists"
elif [ -f "${BAM_FILE%.*}.bai" ]; then
    print_result "pass" "BAM index exists" "Found: ${BAM_FILE%.*}.bai"
else
    print_result "fail" "BAM index not found" "Run: samtools index ${BAM_FILE}" "yes"
fi

# ============================================================================
# TEST 3: BAM Header and Format
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3. BAM Format and Header"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    print_result "warn" "samtools not found" "Load module: module load biokit"
    echo "Skipping remaining tests that require samtools"
    echo ""
    echo "To continue validation, run:"
    echo "  module load biokit"
    echo "  bash $0 ${BAM_FILE}"
    exit 0
fi

# Validate BAM format
if samtools quickcheck -q "${BAM_FILE}"; then
    print_result "pass" "BAM format is valid"
else
    print_result "fail" "BAM format is corrupt" "File may be truncated or damaged" "yes"
fi

# Check sorting
SORT_ORDER=$(samtools view -H "${BAM_FILE}" | grep "^@HD" | grep -oP "SO:\K\w+" || echo "unknown")
if [ "$SORT_ORDER" = "coordinate" ]; then
    print_result "pass" "BAM is coordinate-sorted"
elif [ "$SORT_ORDER" = "queryname" ]; then
    print_result "fail" "BAM is name-sorted" "SpecHLA requires coordinate sorting" "yes"
else
    print_result "warn" "Sorting order unclear" "Verify with: samtools view -H ${BAM_FILE} | grep @HD"
fi

# ============================================================================
# TEST 4: Chromosome Naming Convention
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4. Chromosome Naming"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

CHR_WITH_PREFIX=$(samtools view -H "${BAM_FILE}" | grep "^@SQ" | grep -c "SN:chr" || echo "0")
CHR_WITHOUT_PREFIX=$(samtools view -H "${BAM_FILE}" | grep "^@SQ" | grep -c -E "SN:[0-9XY]" || echo "0")

if [ $CHR_WITH_PREFIX -gt 0 ]; then
    print_result "pass" "Chromosome naming: WITH 'chr' prefix" "Use chr6 in pipeline parameters"
    CHR6_NAME="chr6"
elif [ $CHR_WITHOUT_PREFIX -gt 0 ]; then
    print_result "pass" "Chromosome naming: WITHOUT 'chr' prefix" "Use 6 in pipeline parameters"
    CHR6_NAME="6"
else
    print_result "warn" "Cannot determine chromosome naming" "Check manually"
    CHR6_NAME="chr6"  # default
fi

# List first few chromosomes
echo "     First 5 sequences:"
samtools view -H "${BAM_FILE}" | grep "^@SQ" | head -5 | awk '{print "     "$2}' | sed 's/SN://'

# ============================================================================
# TEST 5: HLA-Aware Reference (CRITICAL)
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "5. HLA-Aware Reference (CRITICAL)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Check for HLA-specific contigs/ALT sequences
HLA_CONTIGS=$(samtools view -H "${BAM_FILE}" | grep -c -E "SN:(HLA-|chr6_.*_alt)" || echo "0")
ALT_CONTIGS=$(samtools view -H "${BAM_FILE}" | grep "^@SQ" | grep -c "_alt" || echo "0")

if [ $HLA_CONTIGS -gt 5 ]; then
    print_result "pass" "HLA ALT contigs found" "Count: ${HLA_CONTIGS} (Excellent!)"
elif [ $ALT_CONTIGS -gt 10 ]; then
    print_result "pass" "ALT contigs present" "Count: ${ALT_CONTIGS} (Good)"
else
    print_result "fail" "No HLA ALT contigs detected" "SpecHLA will have LIMITED accuracy" "yes"
    echo "     ${WARN} This is the MOST CRITICAL issue for SpecHLA"
    echo "     ${INFO} Your BAM appears to use standard reference WITHOUT HLA ALT contigs"
    echo ""
    echo "     Recommendations:"
    echo "       1. Realign to HLA-aware reference (bwakit hs38DH.fa)"
    echo "       2. OR use OptiType/ArcasHLA instead of SpecHLA"
    echo "       3. If continuing, understand results will be less accurate"
fi

# Show some HLA contigs if present
if [ $HLA_CONTIGS -gt 0 ]; then
    echo "     HLA sequences found:"
    samtools view -H "${BAM_FILE}" | grep -E "SN:(HLA-|chr6_.*_alt)" | head -5 | awk '{print "     "$2}' | sed 's/SN://'
fi
echo ""

# ============================================================================
# TEST 6: Coverage in HLA Region
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "6. HLA Region Coverage"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

echo "     Counting reads (this may take a moment)..."

# Total reads
TOTAL_READS=$(samtools view -c "${BAM_FILE}")
echo "     Total reads in BAM: ${TOTAL_READS}"

# HLA region reads
HLA_START=28000000
HLA_END=34000000

HLA_READS=$(samtools view -c "${BAM_FILE}" "${CHR6_NAME}:${HLA_START}-${HLA_END}" 2>/dev/null || echo "0")

if [ $HLA_READS -gt 10000 ]; then
    PERCENT=$(echo "scale=2; ${HLA_READS}*100/${TOTAL_READS}" | bc)
    print_result "pass" "Good HLA coverage" "HLA reads: ${HLA_READS} (${PERCENT}%)"
elif [ $HLA_READS -gt 1000 ]; then
    PERCENT=$(echo "scale=2; ${HLA_READS}*100/${TOTAL_READS}" | bc)
    print_result "warn" "Moderate HLA coverage" "HLA reads: ${HLA_READS} (${PERCENT}%) - Results may vary"
elif [ $HLA_READS -gt 0 ]; then
    PERCENT=$(echo "scale=2; ${HLA_READS}*100/${TOTAL_READS}" | bc)
    print_result "fail" "Low HLA coverage" "HLA reads: ${HLA_READS} (${PERCENT}%) - Insufficient" "no"
else
    print_result "fail" "No HLA reads found" "Check chromosome naming or data type" "yes"
fi

# ============================================================================
# TEST 7: Data Type Assessment
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "7. Data Type Assessment"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ $TOTAL_READS -lt 50000000 ]; then
    DATA_TYPE="Likely EXOME or Targeted Panel"
    SPECHLA_SUITABLE="NO"
    RECOMMENDATION="Use OptiType instead of SpecHLA"
    print_result "fail" "${DATA_TYPE}" "${RECOMMENDATION}" "no"
elif [ $TOTAL_READS -lt 200000000 ]; then
    DATA_TYPE="Likely RNA-seq or Low-coverage WGS"
    SPECHLA_SUITABLE="LIMITED"
    RECOMMENDATION="OptiType or ArcasHLA recommended over SpecHLA"
    print_result "warn" "${DATA_TYPE}" "${RECOMMENDATION}"
else
    DATA_TYPE="Likely Whole Genome Sequencing"
    SPECHLA_SUITABLE="YES"
    RECOMMENDATION="SpecHLA appropriate (if HLA-aware reference used)"
    print_result "pass" "${DATA_TYPE}" "${RECOMMENDATION}"
fi

# ============================================================================
# TEST 8: Read Quality
# ============================================================================
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "8. Read Quality Metrics"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Sample 10000 reads from HLA region
SAMPLE_READS=10000
echo "     Sampling ${SAMPLE_READS} reads from HLA region..."

STATS=$(samtools view "${BAM_FILE}" "${CHR6_NAME}:${HLA_START}-${HLA_END}" 2>/dev/null | \
    head -n ${SAMPLE_READS} | \
    awk '{
        mapped++; 
        if(and($2,4)==0) properly_mapped++;
        mapq_sum+=$5;
    } 
    END {
        if(mapped>0) {
            print properly_mapped, mapped, mapq_sum/mapped
        } else {
            print 0, 0, 0
        }
    }')

PROPERLY_MAPPED=$(echo $STATS | awk '{print $1}')
MAPPED=$(echo $STATS | awk '{print $2}')
AVG_MAPQ=$(echo $STATS | awk '{printf "%.1f", $3}')

if [ "$MAPPED" -gt 0 ]; then
    PROPER_PERCENT=$(echo "scale=1; ${PROPERLY_MAPPED}*100/${MAPPED}" | bc)
    
    if (( $(echo "$AVG_MAPQ >= 30" | bc -l) )); then
        print_result "pass" "Good mapping quality" "Average MAPQ: ${AVG_MAPQ}"
    elif (( $(echo "$AVG_MAPQ >= 20" | bc -l) )); then
        print_result "warn" "Moderate mapping quality" "Average MAPQ: ${AVG_MAPQ}"
    else
        print_result "warn" "Low mapping quality" "Average MAPQ: ${AVG_MAPQ}"
    fi
    
    if (( $(echo "$PROPER_PERCENT >= 80" | bc -l) )); then
        print_result "pass" "Good proper pair rate" "Properly paired: ${PROPER_PERCENT}%"
    else
        print_result "warn" "Some improper pairs" "Properly paired: ${PROPER_PERCENT}%"
    fi
else
    print_result "warn" "Could not sample reads" "Check HLA region accessibility"
fi

# ============================================================================
# Summary and Recommendations
# ============================================================================
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Calculate percentage
PERCENT=$((SCORE * 100 / MAX_SCORE))

echo "Overall Score: ${SCORE}/${MAX_SCORE} (${PERCENT}%)"
echo "Critical Issues: ${CRITICAL_ISSUES}"
echo "Warnings: ${WARNINGS}"
echo ""

if [ $CRITICAL_ISSUES -gt 0 ]; then
    echo -e "${CROSS} CRITICAL ISSUES DETECTED"
    echo "   Your BAM has ${CRITICAL_ISSUES} critical issue(s) that MUST be resolved"
    echo "   before running SpecHLA successfully."
    echo ""
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "RECOMMENDATIONS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ $CRITICAL_ISSUES -eq 0 ] && [ $HLA_CONTIGS -gt 5 ] && [ "$SPECHLA_SUITABLE" = "YES" ]; then
    echo -e "${CHECK} ${GREEN}Your BAM is SUITABLE for SpecHLA analysis${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Ensure reference genome is accessible"
    echo "  2. Run: bash step3_prepare_bam_input.sh"
    echo "  3. Configure reference in step4_run_spechla_bam.sh"
    echo "  4. Submit job: sbatch step4_run_spechla_bam.sh"
    
elif [ $HLA_CONTIGS -eq 0 ]; then
    echo -e "${CROSS} ${RED}NOT RECOMMENDED for SpecHLA${NC}"
    echo ""
    echo "PRIMARY ISSUE: No HLA ALT contigs in reference"
    echo ""
    echo "Options:"
    echo "  ${GREEN}BEST:${NC} Realign your data to HLA-aware reference"
    echo "    • Download bwakit: https://github.com/lh3/bwa/tree/master/bwakit"
    echo "    • Use reference: hs38DH.fa"
    echo "    • Realign with BWA-MEM"
    echo ""
    echo "  ${YELLOW}ALTERNATIVE:${NC} Use different HLA typing tools"
    echo "    • OptiType (best for exome/RNA-seq)"
    echo "    • ArcasHLA (good for RNA-seq)"
    echo "    • Run both for validation"
    echo ""
    echo "  ${RED}NOT RECOMMENDED:${NC} Continue with SpecHLA"
    echo "    • Results will have limited accuracy"
    echo "    • Only do this if you understand the limitations"
    
elif [ "$SPECHLA_SUITABLE" = "NO" ]; then
    echo -e "${CROSS} ${RED}NOT SUITABLE for SpecHLA (Exome/Targeted data)${NC}"
    echo ""
    echo "Recommendations:"
    echo "  ${GREEN}USE INSTEAD:${NC}"
    echo "    • OptiType (designed for exome data)"
    echo "    • ArcasHLA (if RNA-seq)"
    echo ""
    echo "  ${INFO} These tools are specifically optimized for your data type"
    echo "     and will give better results than SpecHLA"
    
elif [ "$SPECHLA_SUITABLE" = "LIMITED" ]; then
    echo -e "${WARN} ${YELLOW}LIMITED suitability for SpecHLA (RNA-seq/Low-cov WGS)${NC}"
    echo ""
    echo "Recommendations:"
    echo "  ${GREEN}PREFERRED:${NC}"
    echo "    • OptiType (RNA mode for RNA-seq)"
    echo "    • ArcasHLA (designed for RNA-seq)"
    echo ""
    echo "  ${YELLOW}OPTIONAL:${NC}"
    echo "    • Run SpecHLA as validation"
    echo "    • Use majority voting with multiple tools"
    echo ""
    echo "  ${INFO} For RNA-seq data, OptiType and ArcasHLA typically"
    echo "     outperform SpecHLA due to expression variability"
else
    echo -e "${WARN} ${YELLOW}PROCEED WITH CAUTION${NC}"
    echo ""
    echo "Issues found:"
    [ $CRITICAL_ISSUES -gt 0 ] && echo "  • ${CRITICAL_ISSUES} critical issue(s)"
    [ $WARNINGS -gt 0 ] && echo "  • ${WARNINGS} warning(s)"
    echo ""
    echo "Review issues above and resolve before proceeding"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Configuration to use in pipeline:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "CHR6_NAME=\"${CHR6_NAME}\""
echo "DATA_TYPE=\"${DATA_TYPE}\""
echo "HLA_READS=${HLA_READS}"
echo "TOTAL_READS=${TOTAL_READS}"
echo ""

# Export results to file
REPORT_FILE="bam_validation_report.txt"
{
    echo "SpecHLA BAM Compatibility Report"
    echo "================================="
    echo "Date: $(date)"
    echo "BAM File: ${BAM_FILE}"
    echo ""
    echo "Results:"
    echo "  Score: ${SCORE}/${MAX_SCORE}"
    echo "  Critical Issues: ${CRITICAL_ISSUES}"
    echo "  Warnings: ${WARNINGS}"
    echo "  Data Type: ${DATA_TYPE}"
    echo "  SpecHLA Suitable: ${SPECHLA_SUITABLE}"
    echo "  Chr6 Name: ${CHR6_NAME}"
    echo "  HLA Contigs: ${HLA_CONTIGS}"
    echo "  HLA Reads: ${HLA_READS}"
    echo "  Total Reads: ${TOTAL_READS}"
} > ${REPORT_FILE}

echo "Report saved to: ${REPORT_FILE}"
echo ""
