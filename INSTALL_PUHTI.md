# Installation Guide for CSC Puhti

This guide will help you set up the HLA Typing Pipeline on CSC Puhti.

## 📦 Prerequisites

You should already have:
- Access to CSC Puhti
- Project number (e.g., project_2008084)
- Singularity containers at `/scratch/project_2008084/hla_references/singularity_cache/containers/`

## 🚀 Installation Steps

### Step 1: Upload Pipeline to Puhti

From your local machine:

```bash
# Download the pipeline from GitHub (after release)
git clone https://github.com/yourusername/hla-typing-pipeline.git

# Or if you have the files locally, upload to Puhti
scp -r hla-typing-pipeline your-username@puhti.csc.fi:/scratch/project_2008084/
```

### Step 2: Set Up on Puhti

Login to Puhti and navigate to the pipeline:

```bash
ssh your-username@puhti.csc.fi
cd /scratch/project_2008084/hla-typing-pipeline
```

### Step 3: Load Required Modules

```bash
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0
```

Add these to your `.bashrc` for automatic loading:

```bash
echo "module load java/21" >> ~/.bashrc
echo "module load biopython-env/3.10.6" >> ~/.bashrc
echo "module load nextflow/25.10.0" >> ~/.bashrc
```

### Step 4: Verify Container Files

Check that container files exist:

```bash
ls -lh /scratch/project_2008084/hla_references/singularity_cache/containers/
```

You should see:
- `optitype.sif`
- `arcashla.sif`
- `spechla.sif`

### Step 5: Test the Pipeline

Run a quick test:

```bash
# Test help message
nextflow run main.nf --help
```

If this shows the help message, installation is successful!

## 🔧 Configuration

### Configure for Your Data

Edit `submit_slurm.sh`:

```bash
nano submit_slurm.sh
```

Update these critical parameters:

```bash
# Your sample directory
INPUT_DIR="/scratch/project_2008084/YOUR_SAMPLES_HERE"

# Input type (bam, cram, or fastq)
INPUT_TYPE="fastq"

# Tools to run
TOOLS="optitype,arcashla"

# Your project number
SLURM_ACCOUNT="project_2008084"
```

**For exome data**, also set:
```bash
SPECHLA_EXON_ONLY=1
TOOLS="spechla"  # or include with other tools
```

### Alternative: Use Config File

Copy and edit the example config:

```bash
cp params.yaml.example params.yaml
nano params.yaml
```

Update the parameters in `params.yaml`, then run:

```bash
nextflow run main.nf -params-file params.yaml -profile puhti,singularity
```

## 📊 Running the Pipeline

### Method 1: Interactive Session

For testing or small jobs:

```bash
# Start interactive session
sinteractive --account project_2008084 --time 02:00:00 --mem 32G --cores 8

# Load modules
module load java/21 biopython-env/3.10.6 nextflow/25.10.0

# Run pipeline
nextflow run main.nf \
    --input /scratch/project_2008084/samples/ \
    --input_type fastq \
    --tools optitype,arcashla \
    --outdir results/ \
    -profile puhti,singularity
```

### Method 2: SLURM Batch Job (Recommended)

For production runs:

```bash
# Create logs directory
mkdir -p logs

# Edit submit script
nano submit_slurm.sh

# Submit job
sbatch submit_slurm.sh

# Check job status
squeue -u $USER

# Monitor progress
tail -f logs/hla_pipeline_*.out
```

### Method 3: Direct Nextflow (For Experts)

```bash
nextflow run main.nf \
    --input /scratch/project_2008084/samples/ \
    --input_type bam \
    --tools optitype,arcashla,spechla \
    --spechla_exon_only 1 \
    --enable_majority_voting \
    --slurm_account project_2008084 \
    --outdir results/ \
    -profile puhti,singularity \
    -resume
```

## 🔍 Monitoring Jobs

### Check SLURM Queue

```bash
# View your jobs
squeue -u $USER

# Detailed job info
scontrol show job JOBID

# Cancel job
scancel JOBID
```

### Check Nextflow Progress

```bash
# View main log
tail -f logs/hla_pipeline_*.out

# View Nextflow log
tail -f .nextflow.log

# View specific sample logs
ls -lh results/*/logs/
cat results/sample_name/optitype/sample_name.optitype.log
```

### Monitor Resources

```bash
# Check disk usage
du -sh results/

# Check your quota
csc-workspaces
```

## 📁 Output Location

Results will be in:
```
results/
├── sample1/
│   ├── optitype/
│   ├── arcashla/
│   └── spechla/
├── sample2/
│   └── ...
├── majority_voting/
└── pipeline_info/
```

## 🆘 Troubleshooting

### Issue: Module Not Found

```bash
# Solution: Load modules
module load java/21 biopython-env/3.10.6 nextflow/25.10.0
```

### Issue: Container Not Found

```bash
# Check container path
ls /scratch/project_2008084/hla_references/singularity_cache/containers/

# Update nextflow.config if needed
nano nextflow.config
# Update: singularity_cache_dir = '/your/path/here'
```

### Issue: Permission Denied

```bash
# Make scripts executable
chmod +x submit_slurm.sh
chmod +x main.nf
```

### Issue: Disk Quota Exceeded

```bash
# Check usage
csc-workspaces

# Clean work directory
rm -rf work/

# Use scratch efficiently
# Request more quota from CSC if needed
```

### Issue: Job Fails Immediately

```bash
# Check SLURM logs
cat logs/hla_pipeline_*.err

# Check Nextflow log
cat .nextflow.log

# Verify input files exist
ls /scratch/project_2008084/samples/
```

## 📚 Common Workflows

### Workflow 1: RNA-seq Analysis

```bash
# Edit submit_slurm.sh
INPUT_DIR="/scratch/project_2008084/rnaseq_samples"
INPUT_TYPE="fastq"
TOOLS="optitype,arcashla"
OPTITYPE_SEQ_TYPE="rna"

# Submit
sbatch submit_slurm.sh
```

### Workflow 2: Exome Analysis

```bash
# Edit submit_slurm.sh
INPUT_DIR="/scratch/project_2008084/exome_bams"
INPUT_TYPE="bam"
TOOLS="spechla"
SPECHLA_EXON_ONLY=1

# Submit
sbatch submit_slurm.sh
```

### Workflow 3: Multi-tool with Consensus

```bash
# Edit submit_slurm.sh
INPUT_DIR="/scratch/project_2008084/wgs_samples"
INPUT_TYPE="fastq"
TOOLS="optitype,arcashla,spechla"
ENABLE_MAJORITY_VOTING="true"

# Submit
sbatch submit_slurm.sh
```

## 🔄 Updates

To update the pipeline:

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Pull latest changes (after GitHub release)
git pull origin main

# Or re-download
rm -rf hla-typing-pipeline
git clone https://github.com/yourusername/hla-typing-pipeline.git
```

## 💡 Best Practices

1. **Always test with 1-2 samples first**
2. **Use `-resume` to recover from failures**
3. **Monitor disk usage regularly**
4. **Keep work directory clean** (can get large)
5. **Use meaningful output directory names**
6. **Document your runs** (keep submission scripts)
7. **Back up important results**

## 📞 Getting Help

### CSC Support
- Email: servicedesk@csc.fi
- Docs: https://docs.csc.fi

### Pipeline Issues
- GitHub: https://github.com/yourusername/hla-typing-pipeline/issues
- Check logs first
- Include error messages
- Describe your data type

## ✅ Quick Reference

```bash
# Load modules
module load java/21 biopython-env/3.10.6 nextflow/25.10.0

# Submit job
sbatch submit_slurm.sh

# Check status
squeue -u $USER

# Monitor log
tail -f logs/hla_pipeline_*.out

# Cancel job
scancel JOBID

# Clean work directory
rm -rf work/

# Check disk usage
du -sh results/
```

## 🎉 You're All Set!

Your HLA typing pipeline is now installed and ready to use on CSC Puhti!

For more details, see:
- [README.md](README.md) - Complete documentation
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [DEPLOYMENT.md](DEPLOYMENT.md) - GitHub release guide
