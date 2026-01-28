# Script Consistency Review - Summary of Changes

## Overview
Reviewed and corrected all BAM preparation and pipeline scripts to ensure consistency across the entire workflow.

## Major Consistency Issues Fixed

### 1. Configuration File Naming
**Before**: Mixed usage of `config_batch_bam.sh` vs `config_bam_batch.sh`
**After**: Standardized to `config_bam_batch.sh` everywhere

### 2. Script Directory Structure
**Before**: Inconsistent script location references
**After**: All scripts now in `${SCRIPTS_DIR}` = `${WORK_DIR}/scripts/`
- All scripts source: `${WORK_DIR}/scripts/config_bam_batch.sh`

### 3. Samplesheet Naming and Usage
**Before**: Confusion about which samplesheet to use
**After**: Clear distinction:
- `${DETAILED_SHEET}` = `samplesheet.csv` (with statistics)
- `${SAMPLE_SHEET}` = `samplesheet_pipeline.csv` (simplified, used by pipeline)

### 4. Character Encoding
**Before**: Mixed unicode characters (✓, ✓, âœ—, etc.)
**After**: Consistent symbols:
- ✓ for success
- ✗ for error
- ⚠ for warning

### 5. Tool Configuration
**Before**: Inconsistent tool lists (optitype vs arcashla)
**After**: DNA-optimized tools: `arcashla,xhla,hlahd`

## File-by-File Changes

### config_bam_batch.sh
**Changes**:
- Added `SCRIPTS_DIR` variable for script location
- Split samplesheet into `SAMPLE_SHEET` (simplified) and `DETAILED_SHEET` (with stats)
- Fixed tool list for DNA compatibility
- Added comprehensive comments
- Standardized export statements
- Added Crypt4GH settings for encrypted files

**New Variables**:
```bash
export SCRIPTS_DIR="${BASE_DIR}/scripts"
export SAMPLE_SHEET="${INPUT_DIR}/samplesheet_pipeline.csv"
export DETAILED_SHEET="${INPUT_DIR}/samplesheet.csv"
export SECRET_KEY="${HLA_REFERENCES}/crypt4gh_keys/ozcanumu.sec"
```

### step1_transfer_bam.sh
**Changes**:
- Updated config sourcing path: `${WORK_DIR}/scripts/config_bam_batch.sh`
- Fixed all checkmark symbols
- Corrected "NEXT STEP" to use `${SCRIPTS_DIR}`
- Improved error messages with consistent formatting
- Added proper passphrase handling

**Key Fixes**:
```bash
source ${WORK_DIR}/scripts/config_bam_batch.sh  # Consistent path
log_msg "  bash ${SCRIPTS_DIR}/step2_submit_array_job.sh"  # Correct next step
```

### step2_prepare_bam_array.sh
**Changes**:
- Updated config sourcing to `${WORK_DIR}/scripts/config_bam_batch.sh`
- Fixed all unicode symbols (✓, ✗, ⚠)
- Added `${LOGS_DIR}` directory creation
- Consistent echo formatting
- Improved error messages

**Key Improvements**:
```bash
mkdir -p "${LOGS_DIR}"  # Create logs directory
echo "✓ BAM file found"  # Consistent symbols
```

### step2b_aggregate_results.sh
**Changes**:
- Updated config sourcing path
- Creates both detailed and simplified samplesheets
- Uses `${DETAILED_SHEET}` for full data
- Uses `${SAMPLE_SHEET}` for simplified pipeline input
- Updated "NEXT STEP" reference to use `${SCRIPTS_DIR}`
- Fixed all symbols

**Critical Changes**:
```bash
# Create detailed samplesheet
echo "sample,bam,bai,total_reads,mapped_reads,hla_reads,hla_status" > "${DETAILED_SHEET}"

# Create simplified samplesheet for pipeline
echo "sample,bam,bai" > "${SAMPLE_SHEET}"
tail -n +2 "${DETAILED_SHEET}" | cut -d',' -f1,2,3 >> "${SAMPLE_SHEET}"
```

### step2_submit_array_job.sh
**Changes**:
- Updated config sourcing path
- Uses `${SCRIPTS_DIR}` for all script references
- Uses `${LOGS_DIR}` in monitoring commands
- Fixed all symbols
- Proper temporary file handling
- Better error messages

**Key Changes**:
```bash
ARRAY_SCRIPT="${SCRIPTS_DIR}/step2_prepare_bam_array.sh"
AGG_SCRIPT="${SCRIPTS_DIR}/step2b_aggregate_results.sh"
echo "  ls -lh ${LOGS_DIR}/prepare_bam_array_${ARRAY_JOB_ID}_*.log"
```

### step3_run_pipeline_bam.sh
**Changes**:
- Updated config sourcing path
- Explicitly uses `${SAMPLE_SHEET}` (simplified samplesheet)
- Added "Please run step2_submit_array_job.sh first" error message
- Uses `${SCRIPTS_DIR}` in all script references
- Fixed all symbols
- Improved output messages
- Added `${DATA_TYPE}` to configuration display

**Critical Changes**:
```bash
source ${WORK_DIR}/scripts/config_bam_batch.sh
if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet not found: ${SAMPLE_SHEET}"
    echo "Please run step2_submit_array_job.sh first"
    exit 1
fi
echo "  sbatch ${SCRIPTS_DIR}/step4_check_results.sh"
```

### step4_check_results.sh
**Changes**:
- Updated config sourcing path
- Uses `${SAMPLE_SHEET}` for sample counting
- Added `${DATA_TYPE}` to summary report
- Fixed all symbols
- Consistent formatting
- Better error messages

**Key Changes**:
```bash
source ${WORK_DIR}/scripts/config_bam_batch.sh
Data Type: ${DATA_TYPE}
Sequence type: ${SEQ_TYPE}
```

### check_bam.sh
**Changes**:
- Updated config sourcing path
- Uses `${SCRIPTS_DIR}` in recommendations
- Fixed all symbols
- Added helpful recommendation at end

**Key Changes**:
```bash
source ${WORK_DIR}/scripts/config_bam_batch.sh
if [ ${without_index} -gt 0 ]; then
    echo "To create missing indices, run:"
    echo "  bash ${SCRIPTS_DIR}/step2_submit_array_job.sh"
fi
```

### README_ARRAY_JOB.md
**Changes**:
- Complete rewrite for consistency
- Added configuration variables section
- Clarified samplesheet usage
- Updated all paths to use proper variables
- Added troubleshooting for samplesheet issues
- Included consistency features section
- Updated workflow examples

## Testing Checklist

Before running in production, verify:

- [ ] Config file is at `${WORK_DIR}/scripts/config_bam_batch.sh`
- [ ] All scripts are in `${WORK_DIR}/scripts/` directory
- [ ] BAM files exist in `${RAW_BAM_DIR}`
- [ ] File list exists: `${BAM_FILE_LIST}`
- [ ] Logs directory will be created: `${LOGS_DIR}`
- [ ] Input directory structure is correct

## Variable Reference Table

| Variable | Value | Used By |
|----------|-------|---------|
| `WORK_DIR` | `/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis` | All scripts |
| `SCRIPTS_DIR` | `${WORK_DIR}/scripts` | All scripts |
| `RAW_BAM_DIR` | `${BASE_DIR}/raw_bam/batch1_VenEx_DNA` | Steps 1, 2 |
| `INPUT_DIR` | `${BASE_DIR}/pipeline_input_bam/${BATCH_NAME}` | Steps 2, 3 |
| `RESULTS_DIR` | `${BASE_DIR}/results_bam/${BATCH_NAME}` | Steps 3, 4 |
| `LOGS_DIR` | `${BASE_DIR}/logs` | All scripts |
| `SAMPLE_SHEET` | `${INPUT_DIR}/samplesheet_pipeline.csv` | Steps 2b, 3, 4 |
| `DETAILED_SHEET` | `${INPUT_DIR}/samplesheet.csv` | Step 2b |
| `BAM_FILE_LIST` | `${BASE_DIR}/batch1_bam_file_list.txt` | Step 1 |

## Pipeline Flow with Correct Files

```
Step 1: Transfer
  └─> BAM files → ${RAW_BAM_DIR}/
  
Step 2: Preparation (Array)
  ├─> Array tasks process BAMs in parallel
  └─> Aggregation creates:
      ├─> ${DETAILED_SHEET}  (samplesheet.csv with stats)
      └─> ${SAMPLE_SHEET}   (samplesheet_pipeline.csv simplified)
  
Step 3: HLA Typing
  ├─> Reads: ${SAMPLE_SHEET} (samplesheet_pipeline.csv)
  └─> Outputs: ${RESULTS_DIR}/
  
Step 4: Analysis
  └─> Summarizes results from ${RESULTS_DIR}/
```

## Benefits of Consistency

1. **No Path Confusion**: All scripts use same variable names
2. **Clear Data Flow**: Know which samplesheet is used where
3. **Easy Debugging**: Consistent logging and symbols
4. **Maintainable**: Changes to config propagate everywhere
5. **Professional**: Clean, readable code
6. **Reliable**: Reduced chance of file not found errors

## Recommended Deployment

```bash
# 1. Create scripts directory
mkdir -p /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# 2. Copy all scripts to scripts directory
cp config_bam_batch.sh step*.sh check_bam.sh \
   /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/

# 3. Make scripts executable
chmod +x /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/*.sh

# 4. Verify configuration
source /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh
echo "Scripts dir: ${SCRIPTS_DIR}"
echo "RAW BAM dir: ${RAW_BAM_DIR}"
echo "Sample sheet: ${SAMPLE_SHEET}"

# 5. Test with check_bam.sh
bash /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/check_bam.sh
```

## Version Control Recommendations

When committing to Git:
1. Tag this as a "consistency update" or "v2.0"
2. Update version number in all scripts
3. Document breaking changes if any
4. Provide migration guide for users of old scripts

## Migration from Old Scripts

If you have results from old scripts:

```bash
# Old samplesheet might be at
OLD_SHEET="${INPUT_DIR}/samplesheet.csv"

# Check format
head ${OLD_SHEET}

# If it's the simple 3-column format, just rename it
mv ${OLD_SHEET} ${INPUT_DIR}/samplesheet_pipeline.csv

# If it has statistics, extract the simple version
head -1 ${OLD_SHEET} | cut -d',' -f1,2,3 > ${INPUT_DIR}/samplesheet_pipeline.csv
tail -n +2 ${OLD_SHEET} | cut -d',' -f1,2,3 >> ${INPUT_DIR}/samplesheet_pipeline.csv
```

## Conclusion

All scripts now:
- Use consistent configuration sourcing
- Reference correct directories
- Use appropriate samplesheets
- Display uniform symbols
- Follow same naming conventions
- Support proper workflow execution

The workflow is now production-ready and maintainable.
