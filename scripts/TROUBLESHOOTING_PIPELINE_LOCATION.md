# Pipeline Location Issue - Troubleshooting Guide

## Problem
The pipeline execution failed with:
```
ERROR: Pipeline not found at /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/hla_typing_pipeline/main.nf
```

## Root Cause
The `PIPELINE_DIR` variable in `config_bam_batch.sh` points to a location where the HLA typing pipeline doesn't exist.

## Quick Fix (3 Steps)

### Step 1: Run Diagnostic
```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
bash diagnose_pipeline_path.sh
```

This will:
- Search for the pipeline in common locations
- Show you where the pipeline actually is
- Give you the exact command to fix it

### Step 2: Fix the Configuration
Once you know the correct path (from Step 1), run:
```bash
bash fix_pipeline_path.sh /path/to/actual/pipeline
```

For example, if the pipeline is at `/scratch/project_2008084/hla_typing_pipeline`:
```bash
bash fix_pipeline_path.sh /scratch/project_2008084/hla_typing_pipeline
```

### Step 3: Rerun the Pipeline
```bash
sbatch step3_run_pipeline_bam.sh
```

## Manual Fix (If Scripts Don't Work)

### 1. Find the Pipeline
```bash
# Search for main.nf in your project space
find /scratch/project_2008084 -name "main.nf" -type f 2>/dev/null

# Common locations to check:
ls -la /scratch/project_2008084/ozcanumu/hla_typing_pipeline/
ls -la /scratch/project_2008084/hla_typing_pipeline/
ls -la ~/hla_typing_pipeline/
```

### 2. Edit Config File
```bash
# Backup the config
cp /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh \
   /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh.backup

# Edit the config
nano /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh

# Find this line:
export PIPELINE_DIR="${BASE_DIR}/hla_typing_pipeline"

# Change it to the correct path, for example:
export PIPELINE_DIR="/scratch/project_2008084/hla_typing_pipeline"
```

### 3. Verify the Fix
```bash
# Source the config and check
source /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh
echo "Pipeline directory: ${PIPELINE_DIR}"
ls -la ${PIPELINE_DIR}/main.nf
```

## If Pipeline Doesn't Exist Anywhere

### Option A: Clone from GitHub
If you have the pipeline in a Git repository:
```bash
cd /scratch/project_2008084/ozcanumu
git clone <your-pipeline-repo-url> hla_typing_pipeline

# Then update config to point to it:
export PIPELINE_DIR="/scratch/project_2008084/ozcanumu/hla_typing_pipeline"
```

### Option B: Copy from Another Location
If the pipeline is on your local machine or another server:
```bash
# From your local machine:
scp -r /path/to/local/hla_typing_pipeline \
    ozcanumu@puhti.csc.fi:/scratch/project_2008084/ozcanumu/

# Then update the config
```

### Option C: Create Symbolic Link
If the pipeline is in a shared location:
```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
ln -s /path/to/actual/pipeline hla_typing_pipeline

# Config will now work as originally written
```

## Expected Pipeline Structure

The pipeline directory should contain:
```
hla_typing_pipeline/
├── main.nf                  # Main Nextflow workflow
├── nextflow.config          # Pipeline configuration
├── modules/                 # Individual tool modules
│   ├── arcashla/
│   ├── xhla/
│   └── hlahd/
├── workflows/               # Workflow definitions
└── bin/                     # Helper scripts
```

Verify with:
```bash
ls -la ${PIPELINE_DIR}
```

## Verification Checklist

After fixing, verify everything is in order:

```bash
# 1. Config file is correct
source /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh
echo "Pipeline: ${PIPELINE_DIR}"

# 2. Pipeline files exist
ls ${PIPELINE_DIR}/main.nf
ls ${PIPELINE_DIR}/nextflow.config

# 3. Samplesheet exists
ls ${SAMPLE_SHEET}

# 4. Can load modules
module load biokit
module load nextflow

# 5. Test nextflow
cd ${PIPELINE_DIR}
nextflow -version
```

## Common Pipeline Locations on Puhti

Based on typical CSC Puhti setups:

1. **Project workspace** (most common):
   ```
   /scratch/project_2008084/hla_typing_pipeline/
   /scratch/project_2008084/pipelines/hla_typing_pipeline/
   ```

2. **User workspace**:
   ```
   /scratch/project_2008084/ozcanumu/hla_typing_pipeline/
   /scratch/project_2008084/ozcanumu/pipelines/hla_typing_pipeline/
   ```

3. **Home directory** (less common, smaller quota):
   ```
   ~/hla_typing_pipeline/
   /users/ozcanumu/hla_typing_pipeline/
   ```

4. **Shared installations** (if available):
   ```
   /projappl/project_2008084/hla_typing_pipeline/
   ```

## Still Having Issues?

### Check Permissions
```bash
ls -la ${PIPELINE_DIR}/
# Ensure you have read access to all files
```

### Check Disk Space
```bash
# Check quota
csc-workspaces
```

### Check Module Availability
```bash
module spider nextflow
module spider singularity
```

### Get More Details
```bash
# Look at the full error log
cat /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/logs/hla_batch_bam_31171872.log

# Check if there are any permission issues
ls -la ${PIPELINE_DIR}/ 2>&1 | head -20
```

## After Fixing

Once you've located and configured the pipeline correctly:

1. **Test with dry run** (optional):
   ```bash
   cd ${PIPELINE_DIR}
   nextflow run main.nf --help
   ```

2. **Resubmit the job**:
   ```bash
   sbatch /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/step3_run_pipeline_bam.sh
   ```

3. **Monitor progress**:
   ```bash
   squeue -u $USER
   tail -f logs/hla_batch_bam_*.log
   ```

## Prevention for Future

To avoid this issue in the future:

1. **Document your pipeline location** in a README
2. **Use absolute paths** in configuration
3. **Create a setup script** that verifies all paths before running
4. **Use environment modules** if available for the pipeline

## Need Help?

If none of these solutions work, gather this information:

```bash
# System info
hostname
whoami
pwd

# Pipeline search results
find /scratch/project_2008084 -name "main.nf" -type f 2>/dev/null

# Config contents
cat /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts/config_bam_batch.sh | grep PIPELINE

# Disk usage
df -h /scratch/project_2008084
```

And reach out to CSC support or your bioinformatics team with this information.
