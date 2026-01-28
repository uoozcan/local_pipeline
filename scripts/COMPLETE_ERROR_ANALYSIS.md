# Complete Error Analysis: Batch BAM HLA Typing Pipeline

## Error Progression Summary

You've encountered 3 sequential errors while running the HLA typing pipeline. Here's the complete journey:

---

## ❌ Error #1: Pipeline Not Found

**Job ID**: 31171872  
**Error**: `ERROR: Pipeline not found at .../hla_typing_pipeline/main.nf`

### Problem
The `PIPELINE_DIR` configuration pointed to a non-existent location.

### Solution
✅ **FIXED** - Pipeline location corrected

### Scripts Provided
- `auto_fix_pipeline.sh` - Automatically find and fix path
- `diagnose_pipeline_path.sh` - Search for pipeline
- `fix_pipeline_path.sh` - Update config with correct path

### Result
Pipeline now found at: `/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline/`

---

## ❌ Error #2: Input Type Parameter

**Job ID**: 31172415  
**Error**: `ERROR ~ No files match pattern *.{bam,cram} at path: .../samplesheet_pipeline.csv/`

### Problem
Wrong `--input_type` parameter. Pipeline treated CSV file as a directory.

**Wrong**:
```bash
--input samplesheet.csv
--input_type bam  ← Expects directory
```

**Correct**:
```bash
--input samplesheet.csv
--input_type csv  ← Reads CSV file
```

### Solution
✅ **FIXED** - Changed input_type from 'bam' to 'csv'

### Scripts Provided
- `quick_fix_input_type.sh` - One-click fix
- `step3_run_pipeline_bam_FIXED.sh` - Pre-corrected script
- `diagnose_samplesheet.sh` - Validate samplesheet
- `check_pipeline_input.sh` - Check pipeline requirements

### Result
Pipeline now correctly reads the CSV samplesheet with 24 samples.

---

## ❌ Error #3: Pipeline Code Error (CURRENT)

**Job ID**: 31172534  
**Error**: `ERROR ~ No such property: ch_input for class: nextflow.script.WorkflowBinding`

### Problem
The pipeline's internal Nextflow code has an error. A channel variable called `ch_input` is being used but was never defined.

**Root Causes**:
1. Pipeline code is incomplete
2. Wrong pipeline version/branch
3. Missing workflow files
4. Pipeline bug

### Solution
🔧 **NEEDS FIX** - Switch to stable pipeline version

### Scripts Provided
- `fix_pipeline_version.sh` - Interactive version switcher
- `diagnose_pipeline_code.sh` - Check pipeline structure
- `TROUBLESHOOTING_PIPELINE_CODE.md` - Complete guide

### Recommended Action
```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash fix_pipeline_version.sh
# Choose Option 1: Switch to latest stable release
sbatch step3_run_pipeline_bam.sh
```

---

## Progress Tracking

| Error | Type | Status | Time to Fix |
|-------|------|--------|-------------|
| #1 - Pipeline Location | Configuration | ✅ Fixed | ~1 min |
| #2 - Input Type | Configuration | ✅ Fixed | ~30 sec |
| #3 - Pipeline Code | Code/Version | 🔧 In Progress | ~1 min |

---

## Complete Fix Workflow

### Step-by-Step Solution

```bash
# Navigate to scripts directory
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Error #1: Already fixed (pipeline location)
# Error #2: Already fixed (input type)

# Error #3: Fix pipeline version
bash fix_pipeline_version.sh
# Select: Option 1 (Switch to latest stable release)

# Verify fix worked
cd ../pipeline
git describe --tags --always

# Resubmit job
cd ../scripts
sbatch step3_run_pipeline_bam.sh

# Monitor
squeue -u $USER
tail -f logs/hla_batch_bam_*.log
```

---

## What Each Error Taught Us

### Error #1 Lesson
✓ Always verify file paths exist before running
✓ Use absolute paths in configuration
✓ Create diagnostic scripts to locate resources

### Error #2 Lesson
✓ Understand pipeline input parameters
✓ `--input_type` depends on input format (directory vs file)
✓ CSV samplesheets need `input_type='csv'`

### Error #3 Lesson
✓ Pipeline code can have bugs or be incomplete
✓ Use stable release tags, not development branches
✓ Version control is crucial for reproducible workflows

---

## Configuration Status

### Current Working Configuration

```bash
# Pipeline
PIPELINE_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline"

# Input
--input samplesheet_pipeline.csv
--input_type csv

# Data
24 samples from batch1_VenEx_DNA_BAM

# Tools
optitype, arcashla, spechla

# Resources
Max CPUs: 40
Max Memory: 180GB
Max Time: 24h
```

### Pending Fix
- Pipeline version: Need to switch to stable release

---

## All Diagnostic Scripts Created

### Quick Fix Scripts
1. `quick_fix_input_type.sh` - Fix input type parameter
2. `auto_fix_pipeline.sh` - Find and fix pipeline path
3. `quick_fix_pipeline_code.sh` - Version info and commands

### Interactive Fix Scripts
1. `fix_pipeline_path.sh` - Interactive path fixer
2. `fix_step3_input_type.sh` - Interactive input type fixer
3. `fix_pipeline_version.sh` - Interactive version manager

### Diagnostic Scripts
1. `diagnose_pipeline_path.sh` - Find pipeline location
2. `diagnose_samplesheet.sh` - Validate samplesheet and BAMs
3. `diagnose_pipeline_code.sh` - Check pipeline structure
4. `check_pipeline_input.sh` - Check pipeline parameters
5. `check_bam.sh` - Quick BAM file check

### Fixed Pipeline Scripts
1. `step3_run_pipeline_bam_FIXED.sh` - Pre-corrected step3
2. `step3_run_pipeline_bam_v2.sh` - Improved with diagnostics

### Documentation
1. `QUICK_FIX_PIPELINE_LOCATION.md` - Error #1 guide
2. `TROUBLESHOOTING_PIPELINE_LOCATION.md` - Error #1 detailed
3. `QUICK_FIX_INPUT_TYPE.md` - Error #2 guide
4. `TROUBLESHOOTING_INPUT_TYPE.md` - Error #2 detailed
5. `QUICK_FIX_PIPELINE_CODE.md` - Error #3 guide
6. `TROUBLESHOOTING_PIPELINE_CODE.md` - Error #3 detailed
7. `ERROR2_COMPLETE_ANALYSIS.md` - Error #2 analysis
8. `CONSISTENCY_REVIEW_SUMMARY.md` - Script consistency fixes
9. `README_ARRAY_JOB.md` - Array job documentation

---

## Pipeline Architecture

```
hla_rnaseq_analysis/
├── scripts/
│   ├── config_bam_batch.sh          ✅ Fixed paths
│   ├── step1_transfer_bam.sh        ✅ Working
│   ├── step2_submit_array_job.sh    ✅ Working
│   ├── step3_run_pipeline_bam.sh    ✅ Fixed input_type
│   ├── step4_check_results.sh       ⏳ Ready
│   └── [all diagnostic scripts]     ✅ Created
├── pipeline/                         🔧 Needs version fix
│   ├── main.nf
│   ├── nextflow.config
│   ├── workflows/
│   └── modules/
├── pipeline_input_bam/
│   └── batch1_VenEx_DNA_BAM/
│       └── samplesheet_pipeline.csv  ✅ Valid (24 samples)
├── raw_bam/                          ✅ BAM files present
└── work_batch_bam/                   ⏳ Work directory
```

---

## Expected Next Steps

### 1. Fix Pipeline Version (~1 minute)
```bash
bash fix_pipeline_version.sh
# Select stable release
```

### 2. Run Pipeline (3-8 hours)
```bash
sbatch step3_run_pipeline_bam.sh
```

### 3. Monitor Progress
```bash
squeue -u $USER
tail -f logs/hla_batch_bam_*.log
```

### 4. Check Results
```bash
sbatch step4_check_results.sh
```

### 5. Download Results
```bash
scp -r ozcanumu@puhti.csc.fi:/.../results_bam/batch1_VenEx_DNA_BAM/ .
```

---

## Success Indicators

### When Everything Works:

```
✓ Pipeline found
✓ Sample sheet found: 24 samples
Input type: csv

N E X T F L O W   ~  version 25.10.0
Launching pipeline...

executor >  slurm (72)
[12/345abc] process > PIPELINE:SAMPLESHEET_CHECK    [100%] 1 of 1 ✔
[ab/cd1234] process > PIPELINE:OPTITYPE:RUN        [ 25%] 6 of 24
[cd/ef5678] process > PIPELINE:ARCASHLA:EXTRACT    [ 12%] 3 of 24
[ef/gh9012] process > PIPELINE:SPECHLA:TYPE        [  8%] 2 of 24
```

---

## Lessons Learned

### Best Practices for Future
1. ✓ Always use stable release tags
2. ✓ Verify all paths before running
3. ✓ Test with small samples first (2-3 samples)
4. ✓ Keep diagnostic scripts handy
5. ✓ Document pipeline version used
6. ✓ Use version control for configurations

### Common Pitfalls to Avoid
1. ✗ Using development branches in production
2. ✗ Relative paths in configuration
3. ✗ Wrong input_type parameter
4. ✗ Skipping validation steps
5. ✗ No version documentation

---

## Support Resources

### When to Use Each Script

| Problem | Script to Run |
|---------|---------------|
| Pipeline not found | `auto_fix_pipeline.sh` |
| Input type error | `quick_fix_input_type.sh` |
| Pipeline code error | `fix_pipeline_version.sh` |
| BAM files missing | `diagnose_samplesheet.sh` |
| Unknown parameter | `check_pipeline_input.sh` |
| Version issues | `fix_pipeline_version.sh` |
| Quick BAM check | `check_bam.sh` |

### Quick Reference Commands

```bash
# Fix everything quickly
cd scripts
bash auto_fix_pipeline.sh          # If pipeline not found
bash quick_fix_input_type.sh       # If input type wrong
bash fix_pipeline_version.sh       # If code error

# Diagnose issues
bash diagnose_pipeline_path.sh     # Find pipeline
bash diagnose_samplesheet.sh       # Check samplesheet
bash diagnose_pipeline_code.sh     # Check pipeline code

# Run pipeline
sbatch step3_run_pipeline_bam.sh   # Main execution
sbatch step4_check_results.sh      # After completion

# Monitor
squeue -u $USER                     # Job status
tail -f logs/hla_batch_bam_*.log   # Real-time log
```

---

## Summary

**Current Status**: 2 of 3 errors fixed, 1 remaining

**Next Action**: Fix pipeline version (1 minute)

**Expected Outcome**: Successful HLA typing of 24 BAM samples

**Total Time**: ~4-10 hours (mostly computation)

**Scripts Available**: 22 diagnostic/fix scripts + documentation

---

## Quick Decision Tree

```
Start
  │
  ├─► Pipeline not found?
  │   └─► bash auto_fix_pipeline.sh
  │
  ├─► Input type error?
  │   └─► bash quick_fix_input_type.sh
  │
  ├─► Pipeline code error (ch_input)?
  │   └─► bash fix_pipeline_version.sh
  │       └─► Select stable release
  │
  └─► Success!
      └─► Monitor and wait for results
```

---

**You are here**: Error #3 - Pipeline code error  
**Next step**: Run `bash fix_pipeline_version.sh`  
**ETA to fix**: ~1 minute  
**ETA to completion**: ~4-10 hours after fix
