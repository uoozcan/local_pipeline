# Scenario 4: RNA-seq FASTQ Analysis (10 Paired Samples)

This tutorial demonstrates HLA typing from **paired-end RNA-seq FASTQ files** using the multi-tool pipeline.

## Scenario Overview

| Parameter | Value |
|-----------|-------|
| **Data Type** | RNA-seq |
| **Input Format** | Paired FASTQ (R1 + R2) |
| **Number of Samples** | 10 (20 files total) |
| **Recommended Tools** | arcasHLA, OptiType, SpecHLA |
| **Expected Runtime** | 3-5 hours |

## Advantages of FASTQ Input

| Advantage | Description |
|-----------|-------------|
| **No alignment needed** | Skip BAM creation step |
| **Smaller files** | FASTQ.gz smaller than BAM |
| **Direct processing** | Tools align to HLA reference directly |

> **Note:** HLA*LA requires BAM input and will be skipped for FASTQ samples.

---

## Step 1: Prepare Your Data

### 1.1 Sample Information

Assume you have 10 RNA-seq samples with paired-end reads:

```
RNASEQ_001_R1.fastq.gz  RNASEQ_001_R2.fastq.gz  - Sample 1
RNASEQ_002_R1.fastq.gz  RNASEQ_002_R2.fastq.gz  - Sample 2
RNASEQ_003_R1.fastq.gz  RNASEQ_003_R2.fastq.gz  - Sample 3
RNASEQ_004_R1.fastq.gz  RNASEQ_004_R2.fastq.gz  - Sample 4
RNASEQ_005_R1.fastq.gz  RNASEQ_005_R2.fastq.gz  - Sample 5
RNASEQ_006_R1.fastq.gz  RNASEQ_006_R2.fastq.gz  - Sample 6
RNASEQ_007_R1.fastq.gz  RNASEQ_007_R2.fastq.gz  - Sample 7
RNASEQ_008_R1.fastq.gz  RNASEQ_008_R2.fastq.gz  - Sample 8
RNASEQ_009_R1.fastq.gz  RNASEQ_009_R2.fastq.gz  - Sample 9
RNASEQ_010_R1.fastq.gz  RNASEQ_010_R2.fastq.gz  - Sample 10
```

### 1.2 Verify FASTQ Format

Before upload, check your files:

```bash
# On local machine
# Check file format
zcat RNASEQ_001_R1.fastq.gz | head -4

# Expected output:
# @SEQ_ID
# GATCGATCGATCGATC...
# +
# FFFFFFFF:FFFFFFF...

# Check read count
zcat RNASEQ_001_R1.fastq.gz | wc -l | awk '{print $1/4 " reads"}'
```

### 1.3 Upload to Puhti

```bash
# Create directory on Puhti first
ssh username@puhti.csc.fi "mkdir -p /scratch/project_XXXXX/$USER/hla_analysis/input_fastq"

# Upload files
scp /path/to/fastq/*.fastq.gz \
    username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_fastq/
```

### 1.4 Verify Upload

```bash
ssh username@puhti.csc.fi
cd /scratch/project_XXXXX/$USER/hla_analysis

# Check files (should have 20 files: 10 R1 + 10 R2)
ls -lh input_fastq/
ls input_fastq/*.fastq.gz | wc -l  # Should be 20

# Verify pairs exist
for r1 in input_fastq/*_R1.fastq.gz; do
    r2="${r1/_R1/_R2}"
    if [ ! -f "$r2" ]; then
        echo "WARNING: Missing R2 for $r1"
    fi
done
```

---

## Step 2: Create Samplesheet

```bash
# Create FASTQ samplesheet
cat > samples_rnaseq_fastq.csv << 'EOF'
sample_id,fastq_1,fastq_2
RNASEQ_001,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_001_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_001_R2.fastq.gz
RNASEQ_002,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_002_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_002_R2.fastq.gz
RNASEQ_003,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_003_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_003_R2.fastq.gz
RNASEQ_004,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_004_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_004_R2.fastq.gz
RNASEQ_005,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_005_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_005_R2.fastq.gz
RNASEQ_006,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_006_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_006_R2.fastq.gz
RNASEQ_007,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_007_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_007_R2.fastq.gz
RNASEQ_008,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_008_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_008_R2.fastq.gz
RNASEQ_009,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_009_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_009_R2.fastq.gz
RNASEQ_010,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_010_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq/RNASEQ_010_R2.fastq.gz
EOF

# Or auto-generate
echo "sample_id,fastq_1,fastq_2" > samples_rnaseq_fastq.csv
for r1 in input_fastq/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${r1/_R1/_R2}"
    echo "${sample},$(realpath $r1),$(realpath $r2)" >> samples_rnaseq_fastq.csv
done
```

---

## Step 3: Create SLURM Script

```bash
cat > submit_rnaseq_fastq_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_rnaseq_fq
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=hla_rnaseq_fq_%j.out
#SBATCH --error=hla_rnaseq_fq_%j.err

# Load modules
module load nextflow
module load singularity

# Set cache
export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline with FASTQ input
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_rnaseq_fastq.csv \
    --outdir results_rnaseq_fastq \
    --tools arcashla,optitype,spechla \
    --seq_type rna \
    --max_cpus 20 \
    --max_memory 64.GB \
    -profile singularity \
    -resume

echo "RNA-seq FASTQ HLA typing completed at $(date)"
EOF
```

### Key FASTQ Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--input_samplesheet` | CSV with `fastq_1,fastq_2` columns | Pipeline auto-detects FASTQ format |
| `--seq_type` | `rna` | Required for OptiType RNA mode |
| `--tools` | No `hlala` | HLA*LA requires BAM input |

---

## Step 4: Submit and Monitor

```bash
# Submit
sbatch submit_rnaseq_fastq_hla.sh

# Monitor
squeue -u $USER
tail -f hla_rnaseq_fq_*.out
```

---

## Step 5: Results Interpretation

### 5.1 View Results

```bash
# Quick summary
for f in results_rnaseq_fastq/RNASEQ_*/RNASEQ_*_consensus.txt; do
    echo "=== $(basename $f _consensus.txt) ==="
    grep -v "^#" "$f"
    echo ""
done
```

### 5.2 Compare with BAM-based Results

If you have both BAM and FASTQ from the same samples:

```bash
# Compare consensus calls
diff results_rnaseq/RNASEQ_001/RNASEQ_001_consensus.txt \
     results_rnaseq_fastq/RNASEQ_001/RNASEQ_001_consensus.txt

# Should be identical or very similar
```

---

## FASTQ-Specific Considerations

### Read Length Requirements

| Read Length | Suitability |
|-------------|-------------|
| ≥100 bp | Optimal |
| 75-99 bp | Good |
| 50-74 bp | Acceptable |
| <50 bp | May reduce accuracy |

### Paired-End vs Single-End

This pipeline requires **paired-end** reads. For single-end data:
- Consider aligning to reference first
- Use BAM workflow instead

### File Naming Convention

The pipeline expects paired files with matching prefixes:

```
✓ Sample_R1.fastq.gz + Sample_R2.fastq.gz
✓ Sample_1.fastq.gz + Sample_2.fastq.gz
✓ Sample.1.fq.gz + Sample.2.fq.gz
✗ Sample_forward.fastq.gz + Sample_reverse.fastq.gz (non-standard)
```

---

## Troubleshooting FASTQ Analysis

| Issue | Cause | Solution |
|-------|-------|----------|
| "R2 file not found" | Mismatched naming | Rename to match R1/R2 pattern |
| Low read count | Poor library | Check FastQC report |
| arcasHLA fails | Corrupted FASTQ | Verify with `gzip -t file.gz` |
| OptiType no results | Wrong seq_type | Use `--seq_type rna` |
