# Scenario 5: WES/WGS FASTQ Analysis (10 Paired Samples)

This tutorial demonstrates HLA typing from **paired-end DNA FASTQ files** (WES or WGS) using the multi-tool pipeline.

## Scenario Overview

| Parameter | Value |
|-----------|-------|
| **Data Type** | WES or WGS (DNA) |
| **Input Format** | Paired FASTQ (R1 + R2) |
| **Number of Samples** | 10 (20 files total) |
| **Recommended Tools** | SpecHLA, HLA-HD, arcasHLA, OptiType |
| **Expected Runtime** | 6-12 hours |

## When to Use FASTQ vs BAM

| Use FASTQ When | Use BAM When |
|----------------|--------------|
| Raw data from sequencer | Already have aligned data |
| Want to skip alignment | Need HLA*LA (requires BAM) |
| Storage is limited | Have existing BAM files |
| Starting fresh analysis | Re-analyzing previous alignments |

---

## Step 1: Prepare Your Data

### 1.1 Sample Information

Assume you have 10 WES samples from a clinical study:

```
Clinical_WES_P001_R1.fastq.gz  Clinical_WES_P001_R2.fastq.gz
Clinical_WES_P002_R1.fastq.gz  Clinical_WES_P002_R2.fastq.gz
Clinical_WES_P003_R1.fastq.gz  Clinical_WES_P003_R2.fastq.gz
Clinical_WES_P004_R1.fastq.gz  Clinical_WES_P004_R2.fastq.gz
Clinical_WES_P005_R1.fastq.gz  Clinical_WES_P005_R2.fastq.gz
Clinical_WES_P006_R1.fastq.gz  Clinical_WES_P006_R2.fastq.gz
Clinical_WES_P007_R1.fastq.gz  Clinical_WES_P007_R2.fastq.gz
Clinical_WES_P008_R1.fastq.gz  Clinical_WES_P008_R2.fastq.gz
Clinical_WES_P009_R1.fastq.gz  Clinical_WES_P009_R2.fastq.gz
Clinical_WES_P010_R1.fastq.gz  Clinical_WES_P010_R2.fastq.gz
```

### 1.2 Check FASTQ Quality

```bash
# On local machine before upload
# Check compression
file Clinical_WES_P001_R1.fastq.gz
# Should show: gzip compressed data

# Check read count
zcat Clinical_WES_P001_R1.fastq.gz | head -4

# Count reads
echo "$(zcat Clinical_WES_P001_R1.fastq.gz | wc -l) / 4" | bc
# WES: typically 50-100 million reads
# WGS: typically 300-500 million reads
```

### 1.3 Upload to Puhti

```bash
# Create directory
ssh username@puhti.csc.fi "mkdir -p /scratch/project_XXXXX/$USER/hla_analysis/input_fastq_dna"

# Upload (use rsync for large WGS files)
rsync -avP --progress /path/to/fastq/*.fastq.gz \
    username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_fastq_dna/
```

---

## Step 2: Create Samplesheet

```bash
ssh username@puhti.csc.fi
cd /scratch/project_XXXXX/$USER/hla_analysis

# Create samplesheet
cat > samples_dna_fastq.csv << 'EOF'
sample_id,fastq_1,fastq_2
Clinical_WES_P001,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P001_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P001_R2.fastq.gz
Clinical_WES_P002,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P002_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P002_R2.fastq.gz
Clinical_WES_P003,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P003_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P003_R2.fastq.gz
Clinical_WES_P004,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P004_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P004_R2.fastq.gz
Clinical_WES_P005,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P005_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P005_R2.fastq.gz
Clinical_WES_P006,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P006_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P006_R2.fastq.gz
Clinical_WES_P007,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P007_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P007_R2.fastq.gz
Clinical_WES_P008,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P008_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P008_R2.fastq.gz
Clinical_WES_P009,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P009_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P009_R2.fastq.gz
Clinical_WES_P010,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P010_R1.fastq.gz,/scratch/project_XXXXX/username/hla_analysis/input_fastq_dna/Clinical_WES_P010_R2.fastq.gz
EOF

# Or auto-generate
echo "sample_id,fastq_1,fastq_2" > samples_dna_fastq.csv
for r1 in input_fastq_dna/*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="${r1/_R1/_R2}"
    echo "${sample},$(realpath $r1),$(realpath $r2)" >> samples_dna_fastq.csv
done
```

---

## Step 3: Create SLURM Script

### For WES Data

```bash
cat > submit_wes_fastq_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_wes_fq
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --output=hla_wes_fq_%j.out
#SBATCH --error=hla_wes_fq_%j.err

module load nextflow
module load singularity

export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline with DNA FASTQ settings
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_dna_fastq.csv \
    --outdir results_wes_fastq \
    --tools spechla,hlahd,arcashla,optitype \
    --seq_type dna \
    --max_cpus 40 \
    --max_memory 128.GB \
    -profile singularity \
    -resume

echo "WES FASTQ HLA typing completed at $(date)"
EOF
```

### For WGS Data

```bash
cat > submit_wgs_fastq_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_wgs_fq
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --output=hla_wgs_fq_%j.out
#SBATCH --error=hla_wgs_fq_%j.err

module load nextflow
module load singularity

export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR
export TMPDIR=/scratch/project_XXXXX/$USER/tmp
mkdir -p $TMPDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline - NOTE: No HLA*LA for FASTQ input
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_dna_fastq.csv \
    --outdir results_wgs_fastq \
    --tools spechla,hlahd,arcashla,optitype,xhla \
    --seq_type dna \
    --max_cpus 40 \
    --max_memory 180.GB \
    -profile singularity \
    -resume

echo "WGS FASTQ HLA typing completed at $(date)"
EOF
```

---

## Step 4: Submit and Monitor

```bash
# Submit appropriate script
sbatch submit_wes_fastq_hla.sh  # For WES
# OR
sbatch submit_wgs_fastq_hla.sh  # For WGS

# Monitor
squeue -u $USER
tail -f hla_*_fq_*.out
```

---

## Step 5: Compare Tool Performance

### 5.1 Tool Agreement Analysis

```bash
# Check how tools agree
for f in results_wes_fastq/Clinical_*/Clinical_*_comparison.txt; do
    echo "=== $(basename $(dirname $f)) ==="
    cat "$f" | column -t
    echo ""
done | head -100
```

### 5.2 Confidence Score Summary

```bash
# Extract confidence scores
echo "Sample,Gene,Confidence" > confidence_summary.csv
for f in results_wes_fastq/Clinical_*/Clinical_*_consensus.txt; do
    sample=$(basename "$f" _consensus.txt)
    grep -v "^#" "$f" | while IFS=$'\t' read gene a1 a2 conf rest; do
        echo "$sample,$gene,$conf"
    done >> confidence_summary.csv
done

# Calculate average confidence per sample
awk -F',' 'NR>1 {sum[$1]+=$3; count[$1]++}
    END {for(s in sum) print s, sum[s]/count[s]}' confidence_summary.csv
```

---

## WES vs WGS FASTQ Comparison

| Aspect | WES FASTQ | WGS FASTQ |
|--------|-----------|-----------|
| **File size** | 5-15 GB per pair | 50-150 GB per pair |
| **Read count** | 50-100M | 300-500M |
| **HLA coverage** | Variable (depends on kit) | Uniform ~30x |
| **Processing time** | 6-12 hours | 18-36 hours |
| **Memory needed** | 128 GB | 180 GB |
| **Tools available** | All except HLA*LA | All except HLA*LA |

---

## Step 6: Clinical Report Generation

### For Transplant Matching

```bash
# Generate HLA typing report for transplant
cat > generate_transplant_report.sh << 'EOF'
#!/bin/bash
sample=$1
consensus="results_wes_fastq/${sample}/${sample}_consensus.txt"

echo "============================================"
echo "HLA TYPING REPORT - TRANSPLANT MATCHING"
echo "============================================"
echo "Patient ID: $sample"
echo "Report Date: $(date +%Y-%m-%d)"
echo ""
echo "HLA GENOTYPE:"
echo "---------------------------------------------"
grep -v "^#" "$consensus" | awk -F'\t' '{
    printf "  %-10s %s / %s\n", $1":", $2, $3
}'
echo ""
echo "============================================"
EOF

chmod +x generate_transplant_report.sh
./generate_transplant_report.sh Clinical_WES_P001
```

### Example Output

```
============================================
HLA TYPING REPORT - TRANSPLANT MATCHING
============================================
Patient ID: Clinical_WES_P001
Report Date: 2025-01-30

HLA GENOTYPE:
---------------------------------------------
  HLA-A:     A*02:01 / A*24:02
  HLA-B:     B*35:01 / B*44:02
  HLA-C:     C*04:01 / C*05:01
  HLA-DRB1:  DRB1*07:01 / DRB1*15:01
  HLA-DQB1:  DQB1*02:02 / DQB1*06:02
  HLA-DPB1:  DPB1*04:01 / DPB1*04:02

============================================
```

---

## Important Notes for FASTQ Input

### 1. HLA*LA Not Available

HLA*LA requires pre-aligned BAM files and cannot be used with FASTQ input:

```bash
# This will generate a warning:
--tools spechla,hlahd,hlala  # hlala will be skipped

# Use this instead:
--tools spechla,hlahd,arcashla,optitype
```

### 2. xHLA Performance

xHLA has limited FASTQ support:

```bash
# For BAM: Full support
# For FASTQ: May have reduced accuracy
# Consider using other tools as primary for FASTQ input
```

### 3. Quality Requirements

| Metric | Minimum | Recommended |
|--------|---------|-------------|
| Read length | 50 bp | ≥100 bp |
| Q30 bases | 75% | ≥85% |
| Read pairs | 10M | ≥50M |

---

## Troubleshooting FASTQ DNA Analysis

| Issue | Cause | Solution |
|-------|-------|----------|
| Very slow processing | Large WGS files | Increase time/memory |
| HLA-HD fails | Low HLA reads | Normal for some kits |
| SpecHLA timeout | Complex HLA region | Increase `--time` |
| Out of disk space | Large temp files | Set `TMPDIR` to scratch |
| Missing pairs | Upload incomplete | Re-upload R2 files |
