# Quick Fix: Pipeline Not Found Error

## The Error
```
ERROR: Pipeline not found at /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_typing_pipeline/main.nf
```

## Fastest Fix (One Command)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash auto_fix_pipeline.sh
```

This script will:
1. ✓ Find your pipeline automatically
2. ✓ Update the configuration
3. ✓ Create a backup
4. ✓ Tell you what to do next

## Alternative: Manual Fix (3 Commands)

```bash
# 1. Find the pipeline
find /scratch/project_2008084 -name "main.nf" -type f 2>/dev/null

# 2. Update config (replace PATH with actual path)
bash fix_pipeline_path.sh /PATH/TO/PIPELINE

# 3. Run the pipeline
sbatch step3_run_pipeline_bam.sh
```

## Check Current Status

```bash
# See what's configured
source /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh
echo $PIPELINE_DIR

# Check if it exists
ls -la $PIPELINE_DIR/main.nf
```

## Common Pipeline Locations

Try these paths (in order of likelihood):
```bash
/scratch/project_2008084/hla_typing_pipeline/
/scratch/project_2008084/ozcanumu/hla_typing_pipeline/
/scratch/project_2008084/pipelines/hla_typing_pipeline/
~/hla_typing_pipeline/
```

## After Fixing

```bash
# Verify the fix
source /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh
ls ${PIPELINE_DIR}/main.nf

# Resubmit the job
sbatch /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/step3_run_pipeline_bam.sh

# Monitor
squeue -u $USER
```

## All Available Tools

| Script | Purpose | When to Use |
|--------|---------|-------------|
| `auto_fix_pipeline.sh` | Automatically find and fix | First choice - easiest |
| `diagnose_pipeline_path.sh` | Find pipeline location | When you want to see options |
| `fix_pipeline_path.sh` | Update config with path | When you know the path |
| `step3_run_pipeline_bam_v2.sh` | Run pipeline (improved) | After fixing the path |

## Quick Diagnostic

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Option 1: Auto-fix everything
bash auto_fix_pipeline.sh

# Option 2: Diagnose first, then fix
bash diagnose_pipeline_path.sh
# (follow the instructions it gives you)

# Option 3: Manual search
find /scratch/project_2008084 -name "main.nf" -type f 2>/dev/null
```

## Expected Output After Fix

```
✓ SUCCESS! Pipeline Configured
Pipeline location: /scratch/project_2008084/hla_typing_pipeline

You can now run the pipeline:
  sbatch /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/step3_run_pipeline_bam.sh
```

## Still Not Working?

1. Check permissions:
   ```bash
   ls -la /scratch/project_2008084/hla_typing_pipeline/
   ```

2. Check if pipeline is complete:
   ```bash
   ls /scratch/project_2008084/hla_typing_pipeline/main.nf
   ls /scratch/project_2008084/hla_typing_pipeline/nextflow.config
   ```

3. See full troubleshooting guide:
   ```bash
   cat TROUBLESHOOTING_PIPELINE_LOCATION.md
   ```

## Contact CSC Support

If all else fails, provide this info:
```bash
# System info
whoami
hostname

# Pipeline search
find /scratch/project_2008084 -name "main.nf" -type f 2>/dev/null

# Current config
cat /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh | grep PIPELINE
```

---

**TL;DR**: Run `bash auto_fix_pipeline.sh` and follow instructions.
