# Scenario 3: Whole Genome Sequencing (WGS) BAM Analysis (10 Samples)

This tutorial demonstrates HLA typing from **WGS BAM files** using the multi-tool pipeline.

## Scenario Overview

| Parameter | Value |
|-----------|-------|
| **Data Type** | Whole Genome Sequencing (WGS) |
| **Input Format** | BAM (aligned to hg38) |
| **Number of Samples** | 10 |
| **Recommended Tools** | SpecHLA, HLA-HD, HLA*LA, xHLA |
| **Expected Runtime** | 16-24 hours |

## Why These Tools for WGS?

| Tool | Why Suitable for WGS |
|------|---------------------|
| **SpecHLA** | High resolution; handles deep coverage well |
| **HLA-HD** | Exhaustive search; excellent for rare alleles |
| **HLA\*LA** | Graph-based; designed specifically for WGS |
| **xHLA** | Fast k-mer approach; good for large WGS files |

> **Note:** HLA*LA works **only** with BAM input and is specifically optimized for WGS data.

---

## Step 1: Prepare Your Data

### 1.1 Sample Information

Assume you have 10 WGS samples from a population genetics study:

```
WGS_Sample_AFR_01.bam  - African ancestry, 30x coverage
WGS_Sample_AFR_02.bam  - African ancestry, 30x coverage
WGS_Sample_EUR_01.bam  - European ancestry, 30x coverage
WGS_Sample_EUR_02.bam  - European ancestry, 30x coverage
WGS_Sample_EAS_01.bam  - East Asian ancestry, 30x coverage
WGS_Sample_EAS_02.bam  - East Asian ancestry, 30x coverage
WGS_Sample_SAS_01.bam  - South Asian ancestry, 30x coverage
WGS_Sample_SAS_02.bam  - South Asian ancestry, 30x coverage
WGS_Sample_AMR_01.bam  - American ancestry, 30x coverage
WGS_Sample_AMR_02.bam  - American ancestry, 30x coverage
```

### 1.2 File Size Considerations

WGS BAM files are large (30-100+ GB each). Plan for:
- **Storage**: ~500 GB - 1 TB for 10 samples
- **Transfer time**: Several hours for upload
- **Processing**: More memory and CPU needed

### 1.3 Upload Files

```bash
# For large files, use rsync with resume capability
rsync -avP --progress /path/to/wgs_bams/*.bam \
    username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam_wgs/

rsync -avP --progress /path/to/wgs_bams/*.bam.bai \
    username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam_wgs/
```

### 1.4 Verify HLA Coverage

```bash
ssh username@puhti.csc.fi
cd /scratch/project_XXXXX/$USER/hla_analysis

module load samtools

# Check HLA region coverage (should be high for 30x WGS)
for bam in input_bam_wgs/*.bam; do
    sample=$(basename "$bam" .bam)
    count=$(samtools view -c "$bam" chr6:28510120-33480577 2>/dev/null)
    echo "$sample: $count HLA reads"
done

# Expected: 50,000-200,000 reads per sample for 30x WGS
```

---

## Step 2: Create Samplesheet

```bash
cat > samples_wgs.csv << 'EOF'
sample_id,bam_path
WGS_Sample_AFR_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_AFR_01.bam
WGS_Sample_AFR_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_AFR_02.bam
WGS_Sample_EUR_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_EUR_01.bam
WGS_Sample_EUR_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_EUR_02.bam
WGS_Sample_EAS_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_EAS_01.bam
WGS_Sample_EAS_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_EAS_02.bam
WGS_Sample_SAS_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_SAS_01.bam
WGS_Sample_SAS_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_SAS_02.bam
WGS_Sample_AMR_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_AMR_01.bam
WGS_Sample_AMR_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wgs/WGS_Sample_AMR_02.bam
EOF
```

---

## Step 3: Create SLURM Script

```bash
cat > submit_wgs_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_wgs
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --output=hla_wgs_%j.out
#SBATCH --error=hla_wgs_%j.err

# Load modules
module load nextflow
module load singularity

# Set cache directories
export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

# Increase temp space for large files
export TMPDIR=/scratch/project_XXXXX/$USER/tmp
mkdir -p $TMPDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline with WGS optimized settings
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_wgs.csv \
    --outdir results_wgs \
    --tools spechla,hlahd,hlala,xhla \
    --reference hg38 \
    --hlala_graph PRG_MHC_GRCh38_withIMGT \
    --max_cpus 40 \
    --max_memory 180.GB \
    -profile singularity \
    -resume

echo "WGS HLA typing completed at $(date)"
EOF
```

### Key WGS Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--tools` | `spechla,hlahd,hlala,xhla` | Optimal for WGS data |
| `--hlala_graph` | `PRG_MHC_GRCh38_withIMGT` | HLA*LA reference graph |
| `--time` | `48:00:00` | WGS takes longer |
| `--mem` | `180G` | Large BAMs need more memory |

---

## Step 4: Submit and Monitor

```bash
# Submit
sbatch submit_wgs_hla.sh

# Monitor with estimated time
squeue -u $USER -o "%.10i %.9P %.20j %.8u %.2t %.10M %.6D %.10l"

# Check Nextflow progress
watch -n 60 'tail -20 .nextflow.log | grep -E "process|Submitted|Completed"'
```

---

## Step 5: Population Genetics Analysis

### 5.1 Allele Frequency by Population

```bash
# Create population frequency summary
echo "Population,Gene,Allele,Count" > population_hla_freq.csv

for pop in AFR EUR EAS SAS AMR; do
    for gene in A B C DRB1 DQB1 DPB1; do
        # Extract alleles for this population and gene
        grep "HLA-${gene}" results_wgs/WGS_Sample_${pop}_*/WGS_*_consensus.txt | \
        awk -F'\t' '{print $2; print $3}' | \
        sort | uniq -c | \
        while read count allele; do
            echo "${pop},HLA-${gene},${allele},${count}"
        done >> population_hla_freq.csv
    done
done
```

### 5.2 Compare Population Differences

```bash
# Find population-specific alleles
echo "=== HLA-B Allele Distribution ==="
for pop in AFR EUR EAS SAS AMR; do
    echo "--- $pop ---"
    grep "HLA-B" results_wgs/WGS_Sample_${pop}_*/WGS_*_consensus.txt | \
    awk -F'\t' '{print $2"\n"$3}' | sort | uniq -c | sort -rn | head -5
done
```

### 5.3 Example Population Results

```
=== HLA-B Allele Distribution ===
--- AFR ---
   2 B*58:01
   2 B*42:01
   1 B*15:10
   1 B*53:01

--- EUR ---
   2 B*07:02
   2 B*44:02
   1 B*08:01
   1 B*51:01

--- EAS ---
   2 B*46:01
   2 B*40:01
   1 B*51:01
   1 B*54:01
```

---

## Step 6: High-Resolution Typing

WGS enables 4-field resolution typing:

```bash
# Re-run with 4-field resolution
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_wgs.csv \
    --outdir results_wgs_4field \
    --tools spechla,hlahd,hlala \
    --reference hg38 \
    --resolution 4-field \
    --max_cpus 40 \
    --max_memory 180.GB \
    -profile singularity \
    -resume
```

### 4-Field vs 2-Field Resolution

| Resolution | Example | Use Case |
|------------|---------|----------|
| 2-field | A*02:01 | Standard clinical use |
| 4-field | A*02:01:01:01 | Research, rare variant detection |

---

## Step 7: Rare Allele Detection

WGS is excellent for detecting rare alleles:

```bash
# Check for rare/novel alleles (those not in top 100)
# These are indicated by HLA-HD's "new" flag or low frequency

grep -r "new\|novel\|rare" results_wgs/*/hlahd/ 2>/dev/null

# Check allele frequencies against reference databases
# Low-frequency alleles may need validation
```

---

## WGS Quality Metrics

### Expected Values for 30x WGS

| Metric | Expected Value |
|--------|---------------|
| HLA region reads | 50,000 - 200,000 |
| Average coverage | 25-35x |
| Confidence scores | >0.90 for most genes |
| Tool agreement | 4/4 expected |

### Check Quality

```bash
# Generate quality summary
for f in results_wgs/*/qc/*_qc_report.txt; do
    sample=$(basename "$f" _qc_report.txt)
    reads=$(grep "HLA region reads" "$f" | awk '{print $NF}')
    echo "$sample: $reads reads"
done | sort -t: -k2 -n
```

---

## WGS-Specific Considerations

### 1. HLA*LA Graph Selection

| Reference | Graph |
|-----------|-------|
| GRCh38/hg38 | `PRG_MHC_GRCh38_withIMGT` |
| GRCh37/hg19 | `PRG_MHC_GRCh37_withIMGT` |

### 2. Deep Coverage Benefits

- Better phasing of heterozygous positions
- Higher confidence in rare allele detection
- More accurate 4-field resolution

### 3. Computational Requirements

| Samples | CPUs | Memory | Time |
|---------|------|--------|------|
| 10 | 40 | 180 GB | 24-48 hours |
| 50 | 40 | 180 GB | 5-7 days |
| 100 | 40 | 180 GB | 10-14 days |

---

## Troubleshooting WGS Analysis

| Issue | Cause | Solution |
|-------|-------|----------|
| HLA*LA timeout | Large BAM file | Increase `--time` to 72:00:00 |
| Out of memory | High coverage | Use `--partition=hugemem` with 256GB |
| Slow processing | Too many samples | Split into batches of 5 |
| Temp space full | Large intermediates | Set `TMPDIR` to scratch |

---

## Research Output Template

```
=================================================================
HLA TYPING REPORT - Population Genetics Study
=================================================================
Sample ID: WGS_Sample_EUR_01
Population: European (EUR)
Coverage: 32x genome-wide
Analysis Date: 2025-01-30

HLA GENOTYPE (4-field resolution):
  HLA-A:    A*01:01:01:01 / A*02:01:01:01   (Confidence: 0.98)
  HLA-B:    B*08:01:01:01 / B*44:02:01:01   (Confidence: 0.96)
  HLA-C:    C*07:01:01:01 / C*05:01:01:01   (Confidence: 0.97)
  HLA-DRB1: DRB1*03:01:01:01 / DRB1*07:01:01:01 (Confidence: 0.94)
  HLA-DQB1: DQB1*02:01:01 / DQB1*03:03:02   (Confidence: 0.92)
  HLA-DPB1: DPB1*04:01:01:01 / DPB1*04:02:01:01 (Confidence: 0.89)

TOOL AGREEMENT:
  SpecHLA:  ✓ All alleles confirmed
  HLA-HD:   ✓ All alleles confirmed
  HLA*LA:   ✓ All alleles confirmed
  xHLA:     ✓ All alleles confirmed (2-field)

QUALITY METRICS:
  HLA Region Reads: 156,234
  Mean Coverage: 32x
  Tool Agreement: 4/4
  Phasing Quality: HIGH
=================================================================
```
