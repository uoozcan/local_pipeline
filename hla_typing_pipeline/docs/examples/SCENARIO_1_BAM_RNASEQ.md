# Scenario 1: RNA-seq BAM Analysis (10 Samples)

This tutorial demonstrates HLA typing from **RNA-seq BAM files** using the multi-tool pipeline.

## Scenario Overview

| Parameter | Value |
|-----------|-------|
| **Data Type** | RNA-seq |
| **Input Format** | BAM (aligned to hg38) |
| **Number of Samples** | 10 |
| **Recommended Tools** | arcasHLA, OptiType, SpecHLA |
| **Expected Runtime** | 4-6 hours |

## Why These Tools for RNA-seq?

| Tool | Why Suitable for RNA-seq |
|------|-------------------------|
| **arcasHLA** | Specifically designed for RNA-seq; uses Kallisto for fast pseudoalignment |
| **OptiType** | Excellent Class I accuracy; works well with expression data |
| **SpecHLA** | High resolution; handles variable expression levels |

> **Note:** HLA*LA and xHLA are optimized for DNA data and may produce suboptimal results with RNA-seq.

---

## Step 1: Prepare Your Data

### 1.1 Sample Information

Assume you have 10 RNA-seq samples from a tumor study:

```
Sample_001.bam  - Patient 1, Tumor tissue
Sample_002.bam  - Patient 2, Tumor tissue
Sample_003.bam  - Patient 3, Tumor tissue
Sample_004.bam  - Patient 4, Tumor tissue
Sample_005.bam  - Patient 5, Tumor tissue
Sample_006.bam  - Patient 6, Normal tissue
Sample_007.bam  - Patient 7, Normal tissue
Sample_008.bam  - Patient 8, Normal tissue
Sample_009.bam  - Patient 9, Normal tissue
Sample_010.bam  - Patient 10, Normal tissue
```

### 1.2 Upload to Puhti

```bash
# From your local machine
scp -r /path/to/rnaseq_bams/*.bam username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam/
scp -r /path/to/rnaseq_bams/*.bam.bai username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam/
```

### 1.3 Verify Upload (on Puhti)

```bash
ssh username@puhti.csc.fi
cd /scratch/project_XXXXX/$USER/hla_analysis

# Check files
ls -lh input_bam/
# Should show ~2-8 GB per BAM file for RNA-seq
```

---

## Step 2: Create Samplesheet

```bash
# Auto-generate samplesheet
cat > samples_rnaseq.csv << 'EOF'
sample_id,bam_path
Sample_001,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_001.bam
Sample_002,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_002.bam
Sample_003,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_003.bam
Sample_004,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_004.bam
Sample_005,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_005.bam
Sample_006,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_006.bam
Sample_007,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_007.bam
Sample_008,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_008.bam
Sample_009,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_009.bam
Sample_010,/scratch/project_XXXXX/username/hla_analysis/input_bam/Sample_010.bam
EOF

# Or auto-generate from files
echo "sample_id,bam_path" > samples_rnaseq.csv
for bam in input_bam/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "${sample},$(realpath $bam)" >> samples_rnaseq.csv
done
```

---

## Step 3: Create SLURM Script

```bash
cat > submit_rnaseq_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_rnaseq
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=hla_rnaseq_%j.out
#SBATCH --error=hla_rnaseq_%j.err

# Load modules
module load nextflow
module load singularity

# Set cache directories
export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline with RNA-seq optimized settings
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_rnaseq.csv \
    --outdir results_rnaseq \
    --tools arcashla,optitype,spechla \
    --reference hg38 \
    --seq_type rna \
    --max_cpus 20 \
    --max_memory 64.GB \
    -profile singularity \
    -resume

echo "RNA-seq HLA typing completed at $(date)"
EOF
```

### Key RNA-seq Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--tools` | `arcashla,optitype,spechla` | Best tools for RNA-seq data |
| `--seq_type` | `rna` | Tells OptiType to use RNA mode |
| `--mem` | `64G` | RNA-seq BAMs are smaller, less memory needed |
| `--time` | `08:00:00` | RNA-seq processes faster |

---

## Step 4: Submit and Monitor

```bash
# Submit job
sbatch submit_rnaseq_hla.sh

# Monitor progress
squeue -u $USER
tail -f hla_rnaseq_*.out
```

---

## Step 5: Interpret Results

### 5.1 Check Consensus Results

```bash
# View all results
for f in results_rnaseq/Sample_*/Sample_*_consensus.txt; do
    echo "=== $(basename $(dirname $f)) ==="
    grep -v "^#" "$f" | head -10
done
```

### 5.2 Example Output

```
=== Sample_001 ===
Gene      Allele1      Allele2      Confidence
HLA-A     A*02:01      A*24:02      0.95
HLA-B     B*07:02      B*40:01      0.88
HLA-C     C*07:02      C*03:04      0.92
HLA-DRB1  DRB1*15:01   DRB1*04:01   0.85
HLA-DQB1  DQB1*06:02   DQB1*03:02   0.78
```

### 5.3 RNA-seq Specific Considerations

1. **Variable Expression**: HLA genes have different expression levels
   - Class I (A, B, C): Usually highly expressed
   - Class II (DR, DQ, DP): Variable, tissue-dependent

2. **Low Confidence Warnings**: If you see low confidence for Class II genes, this is normal for some tissues

3. **Check Read Counts**: Look at QC reports for HLA read counts
   ```bash
   cat results_rnaseq/Sample_001/qc/Sample_001_qc_report.txt
   ```

---

## Step 6: Download Results

```bash
# From your local machine
scp -r username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/results_rnaseq ./

# View HTML report in browser
open results_rnaseq/summary/hla_summary_report.html
```

---

## Expected Output Structure

```
results_rnaseq/
├── Sample_001/
│   ├── Sample_001_consensus.txt      # Final HLA calls
│   ├── Sample_001_comparison.txt     # Tool comparison
│   ├── arcashla/                     # arcasHLA results
│   ├── optitype/                     # OptiType results
│   ├── spechla/                      # SpecHLA results
│   └── visualizations/
│       └── Sample_001_report.html    # Interactive report
├── Sample_002/
│   └── ...
├── summary/
│   ├── hla_summary_report.html       # Multi-sample summary
│   └── hla_summary_statistics.tsv
└── multiqc/
    └── multiqc_report.html
```

---

## Troubleshooting RNA-seq Analysis

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| Low Class II confidence | Low expression in tissue | Normal for some tissues; use Class I results |
| arcasHLA fails | Few HLA reads | Check if tissue expresses HLA |
| OptiType empty results | RNA mode not set | Ensure `--seq_type rna` is used |
| Very few reads | Poor library prep | Check FastQC report |
