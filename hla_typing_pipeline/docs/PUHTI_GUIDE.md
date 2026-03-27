# HLA Typing Pipeline - CSC Puhti User Guide

This guide provides comprehensive instructions for running the HLA Typing Pipeline on CSC Puhti supercomputer.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Initial Setup](#initial-setup)
3. [Running the Pipeline](#running-the-pipeline)
4. [Batch Processing](#batch-processing)
5. [Output Files](#output-files)
6. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Access
- CSC user account with access to Puhti
- Project allocation (e.g., `project_2008084`)
- Sufficient billing units for computation

### Software Requirements (pre-installed on Puhti)
- Nextflow (load via module)
- Singularity (available by default)
- Samtools (for BAM processing)

### Data Requirements
- Input BAM files (aligned to hg38 or hg19) OR paired FASTQ files
- BAM index files (.bai) if using BAM input

---

## Initial Setup

### Step 1: Clone the Repository

```bash
# Connect to Puhti
ssh <username>@puhti.csc.fi

# Navigate to your project scratch space (IMPORTANT: use your actual project ID)
cd /scratch/<YOUR_PROJECT_ID>/<username>

# Clone the repository
git clone https://github.com/uoozcan/local_pipeline.git hla_analysis
cd hla_analysis
```

### Step 2: Configure Your Project

**Option A: Run from project directory (recommended)**
```bash
# The setup script auto-detects project ID from the path
cd /scratch/project_XXXXXXX/$USER/hla_analysis
./hla_typing_pipeline/scripts/puhti_setup.sh
```

**Option B: Set environment variable**
```bash
# Set your project ID
export CSC_PROJECT=project_XXXXXXX

# Run setup
./hla_typing_pipeline/scripts/puhti_setup.sh
```

### Step 3: Configure Container and Database Paths

```bash
# Copy the configuration template
cp hla_typing_pipeline/conf/user.config.template hla_typing_pipeline/conf/user.config

# Edit with your paths
nano hla_typing_pipeline/conf/user.config
```

**Required settings in user.config:**
```groovy
params {
    // Your CSC project ID
    project = 'project_XXXXXXX'

    // Path to Singularity containers
    container_dir = '/scratch/project_XXXXXXX/containers'

    // Path to HLA-HD database
    hlahd_db = '/scratch/project_XXXXXXX/databases/hlahd_db'
}
```

### Step 4: Run Setup Script

```bash
# Make setup script executable and run
chmod +x hla_typing_pipeline/scripts/puhti_setup.sh
./hla_typing_pipeline/scripts/puhti_setup.sh
```

This script will:
- Create necessary directory structure
- Set up Singularity cache
- Verify module availability
- Create a project-specific configuration

### Step 3: Prepare Input Data

#### For BAM Files:
```bash
# Create input directory
mkdir -p input_bam

# Copy or link your BAM files
cp /path/to/your/sample.bam input_bam/
cp /path/to/your/sample.bam.bai input_bam/

# Or create symbolic links
ln -s /path/to/your/sample.bam input_bam/
ln -s /path/to/your/sample.bam.bai input_bam/
```

#### For FASTQ Files:
```bash
# Create input directory
mkdir -p input_fastq

# Copy or link your FASTQ files
ln -s /path/to/your/sample_R1.fastq.gz input_fastq/
ln -s /path/to/your/sample_R2.fastq.gz input_fastq/
```

### Step 4: Create Samplesheet (for multiple samples)

#### BAM Samplesheet (`samples_bam.csv`):
```csv
sample_id,bam_path
Sample1,/scratch/project_id/username/input_bam/Sample1.bam
Sample2,/scratch/project_id/username/input_bam/Sample2.bam
Sample3,/scratch/project_id/username/input_bam/Sample3.bam
```

#### FASTQ Samplesheet (`samples_fastq.csv`):
```csv
sample_id,fastq_1,fastq_2
Sample1,/scratch/project_id/username/input_fastq/Sample1_R1.fastq.gz,/scratch/project_id/username/input_fastq/Sample1_R2.fastq.gz
Sample2,/scratch/project_id/username/input_fastq/Sample2_R1.fastq.gz,/scratch/project_id/username/input_fastq/Sample2_R2.fastq.gz
```

---

## Running the Pipeline

### Option 1: Interactive Run (for testing)

```bash
# Start an interactive session
sinteractive --account project_2008084 --time 04:00:00 --mem 32G --cores 8

# Load required modules
module load nextflow
module load singularity

# Run pipeline with single BAM
nextflow run hla_typing_pipeline/main.nf \
    --input_bam input_bam/sample.bam \
    --outdir results \
    --tools spechla,hlahd \
    -profile singularity
```

### Option 2: SLURM Batch Job (recommended)

#### Single Sample:
```bash
# Edit the submission script with your parameters
nano hla_typing_pipeline/scripts/submit_hla_single.sh

# Submit the job
sbatch hla_typing_pipeline/scripts/submit_hla_single.sh
```

#### Multiple Samples:
```bash
# Edit the batch submission script
nano hla_typing_pipeline/scripts/submit_hla_batch.sh

# Submit the job
sbatch hla_typing_pipeline/scripts/submit_hla_batch.sh
```

### Option 3: Quick Run Script

```bash
# For a single BAM file
./hla_typing_pipeline/scripts/run_hla_puhti.sh --bam input_bam/sample.bam

# For paired FASTQ files
./hla_typing_pipeline/scripts/run_hla_puhti.sh --fastq input_fastq/sample_R1.fastq.gz input_fastq/sample_R2.fastq.gz

# For multiple samples with samplesheet
./hla_typing_pipeline/scripts/run_hla_puhti.sh --samplesheet samples_bam.csv
```

---

## Pipeline Parameters

### Input Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_bam` | Single BAM file | `sample.bam` |
| `--input_fastq_1` | R1 FASTQ file | `sample_R1.fq.gz` |
| `--input_fastq_2` | R2 FASTQ file | `sample_R2.fq.gz` |
| `--input_samplesheet` | CSV with multiple samples | `samples.csv` |
| `--outdir` | Output directory | `./results` |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tools` | `spechla,hlahd` | HLA typing tools: spechla,hlahd,hlala,arcashla,optitype,xhla |
| `--reference` | `hg38` | Reference genome version |
| `--resolution` | `2-field` | Output resolution |
| `--seq_type` | `dna` | Sequence type (for OptiType) |

### Resource Parameters (Puhti-optimized)

| Parameter | Recommended | Description |
|-----------|-------------|-------------|
| `--max_cpus` | `40` | Max CPUs per process |
| `--max_memory` | `180.GB` | Max memory per process |

---

## Batch Processing

### Array Job for Many Samples

For processing many samples efficiently, use the array job script:

```bash
# Generate sample list
ls input_bam/*.bam | sed 's/.*\///' | sed 's/.bam$//' > sample_list.txt

# Submit array job
sbatch --array=1-$(wc -l < sample_list.txt) hla_typing_pipeline/scripts/submit_hla_array.sh
```

### Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# View job output
tail -f slurm-<jobid>.out

# Check specific sample log
cat results/<sample_id>/logs/nextflow.log
```

---

## Output Files

After successful completion, you'll find:

```
results/
├── <sample_id>/
│   ├── qc/
│   │   └── <sample_id>_qc_report.txt      # QC metrics
│   ├── fastqc/
│   │   └── <sample_id>_fastqc.html        # FastQC report
│   ├── spechla/
│   │   └── <sample_id>_spechla.txt        # SpecHLA results
│   ├── hlahd/
│   │   └── <sample_id>_hlahd.txt          # HLA-HD results
│   ├── visualizations/
│   │   ├── <sample_id>_confidence.png     # Confidence chart
│   │   ├── <sample_id>_coverage.png       # Read coverage
│   │   ├── <sample_id>_agreement.png      # Tool agreement
│   │   └── <sample_id>_report.html        # Interactive report
│   ├── <sample_id>_consensus.txt          # Final HLA calls
│   └── <sample_id>_comparison.txt         # Tool comparison
├── summary/
│   ├── hla_summary_report.html            # Multi-sample summary
│   └── hla_summary_statistics.tsv         # Statistics table
├── multiqc/
│   └── multiqc_report.html                # Aggregated QC
└── pipeline_info/
    ├── timeline.html                      # Execution timeline
    └── report.html                        # Pipeline report
```

---

## Resource Recommendations

### For WGS BAM Files (~30-50GB)

```bash
#SBATCH --time=12:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=40
#SBATCH --partition=small
```

### For RNA-seq BAM Files (~5-10GB)

```bash
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=20
#SBATCH --partition=small
```

### For FASTQ Files

```bash
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=40
#SBATCH --partition=small
```

---

## Troubleshooting

### GT-Only Calibration Triage

For calibration runs, treat `conf/wgs_samples_50.txt`, `conf/wes_samples_50.txt`, and
`conf/rna_samples_50.txt` as curated inputs:

- every listed sample must be present in the GT TSV generated by `python3 bin/calibrate_tool_weights.py download-gt`
- GT provenance: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/`
- every listed sample must also be present in the relevant source dataset metadata
- if the runtime GT filter drops a listed sample, fix the curated list in git rather than treating it as expected runtime behavior

Before any new Puhti submissions, sync the repo and check status:

```bash
git pull
bash scripts/run_wgs_calibration_puhti.sh --project project_2008084 --status
bash scripts/run_wes_calibration_puhti.sh --project project_2008084 --status
bash scripts/run_rna_calibration_puhti.sh --project project_2008084 --status
bash scripts/run_e2e_test_puhti.sh --type wgs --project project_2008084 --status
bash scripts/run_e2e_test_puhti.sh --type wes --project project_2008084 --status
bash scripts/run_e2e_test_puhti.sh --type rna --project project_2008084 --status
```

Confirm the GT file exists on scratch and compare it against the curated lists:

```bash
# Upstream GT directory:
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/
python3 bin/calibrate_tool_weights.py download-gt \
  --output /scratch/project_2008084/hla_tools/hla_typing_pipeline/conf/1kgp_hla_gt.tsv
```

If WGS is pending, pull the curated WGS list correction before rerunning extraction or typing.
For WES and RNA, do not resubmit active jobs; inspect existing `typing_*.out/.err`,
`collect_*.out/.err`, and Nextflow trace/report files first, then rerun only if the jobs are stalled or failed.
For a compact non-destructive summary across all three E2E runs, use:
`bash scripts/check_e2e_puhti_status.sh --project project_2008084`

### Common Issues

#### 1. "Module not found" Error
```bash
# Solution: Load modules explicitly
module load nextflow/23.10.0
module load singularity
```

#### 2. "Out of Memory" Error
```bash
# Solution: Increase memory allocation
nextflow run ... --max_memory 256.GB
# Or use bigmem partition for very large files
#SBATCH --partition=bigmem
```

#### 3. "Container not found" Error
```bash
# Solution: Pull container manually
singularity pull --dir $SINGULARITY_CACHEDIR docker://quay.io/biocontainers/spechla:1.0.7
```

#### 4. "BAM index not found" Error
```bash
# Solution: Create BAM index
module load samtools
samtools index sample.bam
```

#### 5. Pipeline Hangs or Slow
```bash
# Check for issues
nextflow log
# Resume from last checkpoint
nextflow run ... -resume
```

### Getting Help

- Check Nextflow logs: `cat .nextflow.log`
- View process logs: `cat work/<hash>/.command.log`
- CSC documentation: https://docs.csc.fi
- Pipeline issues: https://github.com/uoozcan/local_pipeline/issues

---

## Best Practices

1. **Always use scratch space** for computation (not projappl)
2. **Create BAM index** before running the pipeline
3. **Use `-resume`** to continue interrupted runs
4. **Monitor billing units** with `csc-projects`
5. **Clean up work directory** after successful completion:
   ```bash
   rm -rf work/
   ```

---

## Example Workflow

```bash
# 1. Setup (first time only)
cd /scratch/project_2008084/$USER
git clone https://github.com/uoozcan/local_pipeline.git hla_analysis
cd hla_analysis
./hla_typing_pipeline/scripts/puhti_setup.sh

# 2. Prepare input
mkdir -p input_bam
ln -s /path/to/samples/*.bam input_bam/
ln -s /path/to/samples/*.bam.bai input_bam/

# 3. Create samplesheet
./hla_typing_pipeline/scripts/generate_samplesheet.sh input_bam > samples.csv

# 4. Submit job
sbatch hla_typing_pipeline/scripts/submit_hla_batch.sh

# 5. Monitor
squeue -u $USER
tail -f slurm-*.out

# 6. Check results
ls results/*/
cat results/summary/hla_summary_statistics.tsv
```

---

## Version Information

- Pipeline version: 1.2.0
- Nextflow version: >=22.10.0
- Last updated: January 2025
