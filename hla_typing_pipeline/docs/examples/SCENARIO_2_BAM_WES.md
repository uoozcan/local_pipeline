# Scenario 2: Whole Exome Sequencing (WES) BAM Analysis (10 Samples)

This tutorial demonstrates HLA typing from **WES BAM files** using the multi-tool pipeline.

## Scenario Overview

| Parameter | Value |
|-----------|-------|
| **Data Type** | Whole Exome Sequencing (WES) |
| **Input Format** | BAM (aligned to hg38) |
| **Number of Samples** | 10 |
| **Recommended Tools** | SpecHLA, HLA-HD, arcasHLA, OptiType |
| **Expected Runtime** | 8-12 hours |

## Why These Tools for WES?

| Tool | Why Suitable for WES |
|------|---------------------|
| **SpecHLA** | Excellent accuracy with targeted sequencing; handles variable coverage |
| **HLA-HD** | Designed for exome data; exhaustive alignment approach |
| **arcasHLA** | Fast and accurate; works well with exome capture |
| **OptiType** | High Class I accuracy; good for clinical applications |

> **Note:** HLA*LA is designed for WGS and may not perform optimally with WES due to incomplete HLA region coverage.

---

## Step 1: Prepare Your Data

### 1.1 Sample Information

Assume you have 10 WES samples from a pharmacogenomics study:

```
PGx_Patient_01.bam  - Patient 1, Blood sample
PGx_Patient_02.bam  - Patient 2, Blood sample
PGx_Patient_03.bam  - Patient 3, Blood sample
PGx_Patient_04.bam  - Patient 4, Blood sample
PGx_Patient_05.bam  - Patient 5, Blood sample
PGx_Patient_06.bam  - Patient 6, Blood sample
PGx_Patient_07.bam  - Patient 7, Blood sample
PGx_Patient_08.bam  - Patient 8, Blood sample
PGx_Patient_09.bam  - Patient 9, Blood sample
PGx_Patient_10.bam  - Patient 10, Blood sample
```

### 1.2 Check HLA Region Coverage

Before analysis, verify your exome capture kit covers the HLA region:

```bash
# On Puhti, check if HLA region has reads
module load samtools

# For hg38
samtools view -c input_bam/PGx_Patient_01.bam chr6:28510120-33480577

# Expected: >10,000 reads for good HLA typing
# If <1,000 reads, your capture kit may not cover HLA well
```

### 1.3 Upload Files

```bash
# From local machine
scp -r /path/to/wes_bams/*.bam username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam_wes/
scp -r /path/to/wes_bams/*.bam.bai username@puhti.csc.fi:/scratch/project_XXXXX/$USER/hla_analysis/input_bam_wes/
```

---

## Step 2: Create Samplesheet

```bash
ssh username@puhti.csc.fi
cd /scratch/project_XXXXX/$USER/hla_analysis

# Create samplesheet
cat > samples_wes.csv << 'EOF'
sample_id,bam_path
PGx_Patient_01,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_01.bam
PGx_Patient_02,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_02.bam
PGx_Patient_03,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_03.bam
PGx_Patient_04,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_04.bam
PGx_Patient_05,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_05.bam
PGx_Patient_06,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_06.bam
PGx_Patient_07,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_07.bam
PGx_Patient_08,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_08.bam
PGx_Patient_09,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_09.bam
PGx_Patient_10,/scratch/project_XXXXX/username/hla_analysis/input_bam_wes/PGx_Patient_10.bam
EOF
```

---

## Step 3: Create SLURM Script

```bash
cat > submit_wes_hla.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hla_wes
#SBATCH --account=project_XXXXX
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --output=hla_wes_%j.out
#SBATCH --error=hla_wes_%j.err

# Load modules
module load nextflow
module load singularity

# Set cache directories
export SINGULARITY_CACHEDIR=/scratch/project_XXXXX/$USER/.singularity
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

cd /scratch/project_XXXXX/$USER/hla_analysis

# Run pipeline with WES optimized settings
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_wes.csv \
    --outdir results_wes \
    --tools spechla,hlahd,arcashla,optitype \
    --reference hg38 \
    --seq_type dna \
    --min_hla_reads 500 \
    --max_cpus 40 \
    --max_memory 128.GB \
    -profile singularity \
    -resume

echo "WES HLA typing completed at $(date)"
EOF
```

### Key WES Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--tools` | `spechla,hlahd,arcashla,optitype` | Optimal for exome data |
| `--seq_type` | `dna` | DNA sequencing mode |
| `--min_hla_reads` | `500` | Lower threshold due to targeted capture |
| `--mem` | `128G` | WES BAMs are medium-sized |

---

## Step 4: Submit and Monitor

```bash
# Submit
sbatch submit_wes_hla.sh

# Monitor
watch -n 30 'squeue -u $USER'

# Check progress
tail -f hla_wes_*.out
```

---

## Step 5: Interpret Results

### 5.1 Pharmacogenomics Application

For drug hypersensitivity prediction, focus on specific alleles:

```bash
# Check for HLA-B*57:01 (Abacavir hypersensitivity)
grep "B\*57:01" results_wes/*/PGx_*_consensus.txt

# Check for HLA-B*15:02 (Carbamazepine hypersensitivity)
grep "B\*15:02" results_wes/*/PGx_*_consensus.txt

# Check for HLA-A*31:01 (Carbamazepine hypersensitivity in Europeans)
grep "A\*31:01" results_wes/*/PGx_*_consensus.txt
```

### 5.2 Create Pharmacogenomics Summary

```bash
# Generate PGx report
echo "Sample,HLA-B*57:01,HLA-B*15:02,HLA-A*31:01" > pgx_summary.csv
for f in results_wes/PGx_*/PGx_*_consensus.txt; do
    sample=$(basename "$f" _consensus.txt)
    b5701=$(grep -c "B\*57:01" "$f" || echo "0")
    b1502=$(grep -c "B\*15:02" "$f" || echo "0")
    a3101=$(grep -c "A\*31:01" "$f" || echo "0")

    # Convert to Yes/No
    [ "$b5701" -gt 0 ] && b5701="YES" || b5701="NO"
    [ "$b1502" -gt 0 ] && b1502="YES" || b1502="NO"
    [ "$a3101" -gt 0 ] && a3101="YES" || a3101="NO"

    echo "$sample,$b5701,$b1502,$a3101" >> pgx_summary.csv
done

cat pgx_summary.csv
```

### 5.3 Example Output

```
Sample,HLA-B*57:01,HLA-B*15:02,HLA-A*31:01
PGx_Patient_01,NO,NO,NO
PGx_Patient_02,YES,NO,NO        <- Abacavir risk!
PGx_Patient_03,NO,NO,YES        <- Carbamazepine risk (Europeans)
PGx_Patient_04,NO,NO,NO
PGx_Patient_05,NO,YES,NO        <- Carbamazepine risk (Asians)
...
```

---

## Step 6: Quality Considerations for WES

### 6.1 Coverage Variability

WES has variable coverage across HLA region due to capture probe design:

```bash
# Check per-sample coverage
for f in results_wes/*/qc/*_qc_report.txt; do
    echo "=== $(basename $(dirname $(dirname $f))) ==="
    grep "HLA region reads" "$f"
done
```

### 6.2 Interpret Confidence Scores

| Confidence | Interpretation for WES |
|------------|----------------------|
| >0.90 | High confidence, reliable call |
| 0.70-0.90 | Good, typical for WES |
| 0.50-0.70 | Moderate, may have coverage gaps |
| <0.50 | Low, check if HLA region was captured |

### 6.3 Tool Agreement

For clinical applications, require at least 3/4 tools to agree:

```bash
# Check tool agreement
cat results_wes/PGx_Patient_01/PGx_Patient_01_comparison.txt
```

---

## WES-Specific Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Very few HLA reads | HLA not in capture kit | Check kit specifications |
| Uneven confidence | Variable probe coverage | Normal for WES; use consensus |
| Missing Class II | Poor DRB/DQB capture | Some kits focus on Class I |
| Tool disagreement | Low coverage regions | Trust tools with higher read counts |

---

## Clinical Reporting Template

```
=================================================================
HLA TYPING REPORT - Pharmacogenomics Screening
=================================================================
Sample ID: PGx_Patient_02
Analysis Date: 2025-01-30
Pipeline Version: 1.2.0
Tools Used: SpecHLA, HLA-HD, arcasHLA, OptiType

HLA GENOTYPE:
  HLA-A:  A*02:01 / A*24:02  (Confidence: 0.95)
  HLA-B:  B*57:01 / B*44:02  (Confidence: 0.92)  *** ALERT ***
  HLA-C:  C*06:02 / C*05:01  (Confidence: 0.88)

PHARMACOGENOMIC ALERTS:
  [!] HLA-B*57:01 DETECTED
      - Associated with Abacavir hypersensitivity
      - Recommendation: Avoid Abacavir prescription
      - Reference: CPIC Guidelines

QUALITY METRICS:
  HLA Region Reads: 15,234
  Tool Agreement: 4/4
  Overall Confidence: HIGH
=================================================================
```
