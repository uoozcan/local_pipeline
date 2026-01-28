# QUICK FIX: Pipeline Code Error (ch_input)

## Your Error
```
ERROR ~ No such property: ch_input for class: nextflow.script.WorkflowBinding
```

## What It Means
The pipeline code itself has a bug or is incomplete. The variable `ch_input` doesn't exist in the workflow.

## Fastest Fix (1 minute)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash fix_pipeline_version.sh
```

Then choose **Option 1**: Switch to latest stable release

This will fix the pipeline code by switching to a working version.

---

## Alternative: Manual Fix

### Check Current Version
```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline
git describe --tags --always
git branch --show-current
```

### Switch to Stable Version
```bash
# See available versions
git tag -l | sort -V | tail -10

# Switch to latest stable (e.g., v1.0.0)
git checkout v1.0.0

# Or switch to main branch
git checkout main
git pull
```

### Try Pipeline Again
```bash
sbatch /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/step3_run_pipeline_bam.sh
```

---

## Diagnose First (Optional)

```bash
bash diagnose_pipeline_code.sh
```

This shows:
- Pipeline structure
- Missing files
- Git version info
- Where `ch_input` is referenced

---

## Your Progress

| Issue | Status |
|-------|--------|
| 1. Pipeline location | ✅ Fixed |
| 2. Input type parameter | ✅ Fixed |
| 3. Pipeline code | 🔧 Fixing now |

---

## Why This Happens

### Development vs Stable
You might be on a **development branch** with incomplete features.

### Solution
Switch to a **stable release tag** that has been tested.

---

## Common Scenarios

### Scenario 1: On Development Branch
```bash
# Current: dev or feature/something
git branch --show-current

# Fix: Switch to stable
git checkout v1.0.0
```

### Scenario 2: Outdated Code
```bash
# Fix: Update to latest
git checkout main
git pull
```

### Scenario 3: Corrupted Pipeline
```bash
# Fix: Clone fresh
mv pipeline pipeline.backup
git clone <REPO_URL> pipeline
```

---

## After Fixing

### Verify Version
```bash
cd pipeline
git describe --tags --always
# Should show a version tag like: v1.0.0
```

### Resubmit Job
```bash
sbatch scripts/step3_run_pipeline_bam.sh
```

### Monitor
```bash
tail -f logs/hla_batch_bam_*.log
```

### Success Looks Like
```
N E X T F L O W   ~  version 25.10.0
executor >  slurm (72)
[xx/xxxxxx] process > PIPELINE:SAMPLESHEET_CHECK  [100%] 1 of 1 ✔
[xx/xxxxxx] process > PIPELINE:OPTITYPE:RUN       [  5%] 1 of 24
```

---

## If Still Failing

### Get More Info
```bash
# Check pipeline structure
ls pipeline/main.nf pipeline/workflows/

# Check detailed error
cat work_batch_bam/.nextflow.log | tail -50

# Run full diagnostic
bash diagnose_pipeline_code.sh
```

### Get Help
Provide this info:
```bash
cd pipeline
git describe --tags --always
git remote get-url origin
ls -lh workflows/
```

---

## Quick Commands

```bash
# Fix version
bash fix_pipeline_version.sh

# Diagnose
bash diagnose_pipeline_code.sh

# Check structure
ls pipeline/main.nf pipeline/workflows/*.nf

# View error
cat work_batch_bam/.nextflow.log | grep -A 10 ch_input

# Resubmit
sbatch scripts/step3_run_pipeline_bam.sh
```

---

## Summary

**Problem**: Pipeline code bug  
**Cause**: Wrong version or incomplete code  
**Solution**: Switch to stable release  
**Time**: ~1 minute  

---

**Bottom Line**: Run `bash fix_pipeline_version.sh` and choose option 1
