# HLA Typing Pipeline - Download Guide

## 📦 What to Download

You have two options to get your pipeline files:

### Option 1: Download Complete Archive (Recommended)

**[hla-typing-pipeline-v2.0.0.tar.gz](computer:///mnt/user-data/outputs/hla-typing-pipeline-v2.0.0.tar.gz)** (28 KB)

This contains everything you need:
- Main pipeline script
- All modules
- Configuration files
- Complete documentation
- SLURM submission script

**To extract on your computer:**
```bash
tar -xzf hla-typing-pipeline-v2.0.0.tar.gz
cd hla-typing-pipeline
```

**To extract on CSC Puhti:**
```bash
# After uploading to Puhti
tar -xzf hla-typing-pipeline-v2.0.0.tar.gz
cd hla-typing-pipeline
```

### Option 2: Download Individual Files

If you prefer to download files separately:

#### Core Pipeline Files

1. **[main.nf](computer:///mnt/user-data/outputs/main.nf)** (11 KB)
   - Main pipeline script

2. **[nextflow.config](computer:///mnt/user-data/outputs/nextflow.config)** (8.3 KB)
   - Pipeline configuration with CSC Puhti settings

3. **[submit_slurm.sh](computer:///mnt/user-data/outputs/submit_slurm.sh)** (4.4 KB)
   - SLURM job submission script (executable)

4. **[modules.tar.gz](computer:///mnt/user-data/outputs/modules.tar.gz)** (8.4 KB)
   - All Nextflow modules (optitype, arcashla, spechla, etc.)

#### Documentation Files

5. **[README.md](computer:///mnt/user-data/outputs/README.md)** (11 KB)
   - Complete pipeline documentation

6. **[QUICKSTART.md](computer:///mnt/user-data/outputs/QUICKSTART.md)** (5.4 KB)
   - 5-minute quick start guide

7. **[INSTALL_PUHTI.md](computer:///mnt/user-data/outputs/INSTALL_PUHTI.md)** (7.2 KB)
   - CSC Puhti installation guide

8. **[RELEASE_SUMMARY.md](computer:///mnt/user-data/outputs/RELEASE_SUMMARY.md)** (8.1 KB)
   - Release overview and what's included

## 🚀 Quick Setup After Download

### If you downloaded the complete archive:

```bash
# 1. Extract the archive
tar -xzf hla-typing-pipeline-v2.0.0.tar.gz
cd hla-typing-pipeline

# 2. Make submit script executable
chmod +x submit_slurm.sh

# 3. Upload to CSC Puhti
scp -r hla-typing-pipeline ozcanumu@puhti.csc.fi:/scratch/project_2008084/
```

### If you downloaded individual files:

```bash
# 1. Create directory structure
mkdir -p hla-typing-pipeline
cd hla-typing-pipeline

# 2. Place main.nf, nextflow.config, submit_slurm.sh here

# 3. Extract modules
tar -xzf modules.tar.gz

# 4. Make submit script executable
chmod +x submit_slurm.sh

# 5. Upload to CSC Puhti
cd ..
scp -r hla-typing-pipeline ozcanumu@puhti.csc.fi:/scratch/project_2008084/
```

## 📝 Next Steps on CSC Puhti

Once uploaded to Puhti:

```bash
# 1. SSH to Puhti
ssh ozcanumu@puhti.csc.fi

# 2. Navigate to pipeline
cd /scratch/project_2008084/hla-typing-pipeline

# 3. Load required modules
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0

# 4. Edit submit_slurm.sh with your data paths
nano submit_slurm.sh

# 5. Create logs directory
mkdir -p logs

# 6. Submit test job (with 1-2 samples first!)
sbatch submit_slurm.sh
```

## 🔍 What You Need to Configure

Before running, edit `submit_slurm.sh` and update:

```bash
# Your input data location
INPUT_DIR="/scratch/project_2008084/YOUR_SAMPLES_HERE"

# Input type: bam, cram, or fastq
INPUT_TYPE="fastq"

# Tools to run
TOOLS="optitype,arcashla"

# For RNA-seq
OPTITYPE_SEQ_TYPE="rna"

# For exome data - CRITICAL!
SPECHLA_EXON_ONLY=1
```

## ⚠️ Important Notes

1. **Test First**: Always test with 1-2 samples before running on full dataset

2. **Exome Data**: If using exome data with SpecHLA, MUST set:
   ```bash
   SPECHLA_EXON_ONLY=1
   ```

3. **Container Files**: The pipeline is configured to use your pre-downloaded containers at:
   ```
   /scratch/project_2008084/hla_references/singularity_cache/containers/
   ```

4. **Resume Jobs**: Use `-resume` flag to continue failed runs

## 📚 Documentation to Read

**Start here:**
1. RELEASE_SUMMARY.md - Overview of what's included
2. QUICKSTART.md - 5-minute quick start
3. README.md - Complete documentation
4. INSTALL_PUHTI.md - Puhti-specific setup

## 🆘 Troubleshooting Downloads

If you can't download files:

1. Try clicking the links above
2. Check your browser's download settings
3. Try a different browser
4. Look for files in your Downloads folder

## ✅ Files You Should Have

After extraction, verify you have:

```
hla-typing-pipeline/
├── main.nf
├── nextflow.config
├── submit_slurm.sh
├── modules/
│   ├── optitype.nf
│   ├── arcashla.nf
│   ├── spechla.nf
│   ├── bam_to_fastq.nf
│   ├── aggregation.nf
│   └── majority_voting.nf
└── [documentation files]
```

## 🎯 Ready to Run?

Once you have:
- ✅ Downloaded and extracted files
- ✅ Uploaded to CSC Puhti
- ✅ Edited submit_slurm.sh with your paths
- ✅ Loaded required modules

Then:
```bash
sbatch submit_slurm.sh
```

Good luck! 🧬
