# Pipeline Code Error: ch_input Not Found

## Error Message
```
ERROR ~ No such property: ch_input for class: nextflow.script.WorkflowBinding
```

## What This Means

This is a **pipeline code error**, not a configuration error. The pipeline's Nextflow code is trying to access a channel variable named `ch_input` that doesn't exist in the workflow.

### Progress So Far
✅ Pipeline location fixed  
✅ Input type parameter fixed (`csv`)  
✅ Samplesheet format correct  
❌ Pipeline code has internal error  

## Root Causes

### 1. Incomplete Pipeline Code (Most Likely)
The pipeline workflows may be missing or incomplete.

### 2. Wrong Pipeline Version/Branch
You may be on a development branch with incomplete features.

### 3. Missing Workflow Files
The main workflow file that defines `ch_input` is missing.

### 4. Pipeline Bug
There's a bug in the pipeline code itself.

## Quick Diagnostic

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash diagnose_pipeline_code.sh
```

This will check:
- Pipeline structure
- Missing files
- Git version/branch
- Channel definitions

## Solutions (Try in Order)

### Solution 1: Switch to Stable Release (RECOMMENDED)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash fix_pipeline_version.sh
# Choose option 1: Switch to latest stable release
```

This will:
- Check available release versions
- Switch to the latest stable tag
- Update the pipeline code

### Solution 2: Update Pipeline to Latest

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline
git pull origin main  # or master
sbatch ../scripts/step3_run_pipeline_bam.sh
```

### Solution 3: Clone Fresh Pipeline

If the pipeline is corrupted:

```bash
# Backup current pipeline
mv /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline \
   /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline.backup

# Clone fresh copy (replace URL with your pipeline repo)
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
git clone <YOUR_PIPELINE_REPO_URL> pipeline

# Try again
sbatch scripts/step3_run_pipeline_bam.sh
```

### Solution 4: Check Pipeline Documentation

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline
cat README.md
cat docs/usage.md
```

Check for:
- Required dependencies
- Special setup steps
- Known issues

## Manual Investigation

### Check Pipeline Structure

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline

# Essential files
ls -lh main.nf nextflow.config

# Workflow files (where ch_input should be defined)
ls -lh workflows/

# Module files
ls -lh modules/local/
```

### Search for ch_input Definition

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline

# Where is ch_input used?
grep -r "ch_input" *.nf workflows/*.nf

# Where should it be defined?
grep -r "ch_input =" *.nf workflows/*.nf
```

### Check Git Status

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline

# Current version
git describe --tags --always

# Current branch
git branch --show-current

# Available versions
git tag -l | sort -V | tail -10

# Recent changes
git log --oneline -10
```

## Understanding the Error

### What is `ch_input`?
In Nextflow, `ch_input` is typically a **channel** that contains input data. Channels are like pipes that carry data between processes.

### Why It's Missing
```groovy
// Somewhere in the pipeline, there's code like:
workflow {
    // ... some code ...
    some_process(ch_input)  // ← Trying to use ch_input
}

// But ch_input was never created:
// ch_input = Channel.fromPath(...)  // ← This line is missing!
```

### What Should Happen
```groovy
workflow {
    // Create the channel
    ch_input = SAMPLESHEET_CHECK(samplesheet)
    
    // Use the channel
    some_process(ch_input)  // ← Now it exists
}
```

## Common Pipeline Issues

### Issue 1: Incomplete Workflow File

**Check**: Does `workflows/hla_typing.nf` exist?
```bash
ls -lh pipeline/workflows/hla_typing.nf
```

**Fix**: Ensure all workflow files are present

### Issue 2: Wrong Pipeline Branch

**Check**: Are you on a development branch?
```bash
cd pipeline
git branch --show-current
```

**Fix**: Switch to main/master or a release tag

### Issue 3: Missing Dependencies

**Check**: Are required modules present?
```bash
ls pipeline/modules/local/
```

**Fix**: Ensure pipeline is complete

## Nextflow.log Details

For more information:
```bash
cat /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/work_batch_bam/.nextflow.log
```

Look for:
- Which file has the error
- Line number of the error
- Stack trace

## If Using Custom Pipeline

If this is a custom or modified pipeline:

### Check Your Code
```bash
cd pipeline

# Find where ch_input is used but not defined
grep -n "ch_input" main.nf workflows/*.nf | grep -v "="
```

### Add Channel Definition
You may need to add:
```groovy
// In main.nf or workflow file
ch_input = Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
```

## Alternative: Use Different HLA Pipeline

If this pipeline has persistent issues, consider:

1. **nf-core/hlatyping** - Well-maintained community pipeline
2. **Custom scripts** - Use individual HLA tools directly
3. **Different pipeline** - Ask for working pipeline from colleagues

## Prevention

### Before Running Pipeline:
1. Check Git status: `git describe --tags`
2. Verify structure: `ls main.nf workflows/ modules/`
3. Read documentation: `cat README.md`
4. Test with small sample: Use 2-3 samples first

### Use Stable Versions:
```bash
# Always use tagged releases
git tag -l
git checkout v1.0.0  # or latest stable
```

## Getting Help

### Gather This Information:

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline

# Pipeline version
git describe --tags --always

# Pipeline structure
ls -R | head -50

# Error context
cat ../work_batch_bam/.nextflow.log | grep -A 10 -B 10 "ch_input"

# Workflow files
ls -lh workflows/
cat workflows/*.nf | head -50
```

### Contact:
- Pipeline maintainer (check README for contact)
- Your bioinformatics team
- CSC support for Puhti-specific issues

## Expected Output After Fix

When working correctly:
```
N E X T F L O W   ~  version 25.10.0
Launching pipeline...

Input type          : csv
executor >  slurm (72)
[12/345abc] process > PIPELINE:SAMPLESHEET_CHECK    [100%] 1 of 1 ✔
[ab/cd1234] process > PIPELINE:OPTITYPE:INDEX      [  0%] 0 of 24
[cd/ef5678] process > PIPELINE:ARCASHLA:EXTRACT    [  0%] 0 of 24
```

## Summary

| Issue | Status | Solution |
|-------|--------|----------|
| Pipeline location | ✅ Fixed | - |
| Input type | ✅ Fixed | Changed to `csv` |
| Pipeline code | ❌ Error | Try stable release |

**Next Step**: Run `bash fix_pipeline_version.sh` and select option 1 (stable release)

## Quick Reference

```bash
# Diagnose
bash diagnose_pipeline_code.sh

# Fix version
bash fix_pipeline_version.sh

# Check structure
ls pipeline/main.nf pipeline/workflows/

# Check Git
cd pipeline && git status

# View detailed error
cat work_batch_bam/.nextflow.log
```

---

**TL;DR**: Pipeline code is incomplete/wrong version. Run `bash fix_pipeline_version.sh` and switch to latest stable release.
