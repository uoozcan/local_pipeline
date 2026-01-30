# HLA Typing Pipeline: Scenario Comparison Summary

This document provides a quick reference for choosing the right analysis approach based on your data type.

---

## Quick Decision Guide

```
What is your input format?
│
├── BAM files
│   │
│   └── What type of sequencing?
│       │
│       ├── RNA-seq ──────► Scenario 1: Use arcashla,optitype,spechla
│       │                   Time: 4-6 hours | Memory: 64GB
│       │
│       ├── WES ──────────► Scenario 2: Use spechla,hlahd,arcashla,optitype
│       │                   Time: 8-12 hours | Memory: 128GB
│       │
│       └── WGS ──────────► Scenario 3: Use spechla,hlahd,hlala,xhla
│                           Time: 16-24 hours | Memory: 180GB
│
└── FASTQ files
    │
    └── What type of sequencing?
        │
        ├── RNA-seq ──────► Scenario 4: Use arcashla,optitype,spechla
        │                   Time: 3-5 hours | Memory: 64GB
        │
        └── WES/WGS ──────► Scenario 5: Use spechla,hlahd,arcashla,optitype
                            Time: 6-24 hours | Memory: 128-180GB
```

---

## Tool Availability by Input Type

| Tool | BAM RNA-seq | BAM WES | BAM WGS | FASTQ RNA-seq | FASTQ DNA |
|------|:-----------:|:-------:|:-------:|:-------------:|:---------:|
| **SpecHLA** | ✓ | ✓ | ✓ | ✓ | ✓ |
| **HLA-HD** | ✓ | ✓ | ✓ | ✓ | ✓ |
| **HLA\*LA** | △ | △ | ✓ | ✗ | ✗ |
| **arcasHLA** | ✓ | ✓ | ✓ | ✓ | ✓ |
| **OptiType** | ✓ | ✓ | ✓ | ✓ | ✓ |
| **xHLA** | △ | ✓ | ✓ | △ | △ |

✓ = Recommended | △ = Supported but not optimal | ✗ = Not supported

---

## Recommended Tool Combinations

### By Data Type

| Data Type | Recommended Tools | Notes |
|-----------|-------------------|-------|
| **RNA-seq** | `arcashla,optitype,spechla` | arcasHLA is purpose-built for RNA-seq |
| **WES** | `spechla,hlahd,arcashla,optitype` | 4 tools for high confidence |
| **WGS** | `spechla,hlahd,hlala,xhla` | HLA*LA excels with WGS |

### By Use Case

| Use Case | Tools | Why |
|----------|-------|-----|
| **Speed priority** | `arcashla,optitype` | Fastest combination |
| **Accuracy priority** | `spechla,hlahd,hlala,arcashla` | Maximum consensus |
| **Rare allele detection** | `spechla,hlahd` | Best for novel variants |
| **Class I only** | `optitype,arcashla` | Focused typing |
| **Transplant matching** | `spechla,hlahd,arcashla,optitype` | High confidence required |
| **Pharmacogenomics** | `spechla,hlahd,optitype` | Clinical-grade accuracy |

---

## Resource Requirements

### By Scenario

| Scenario | CPUs | Memory | Time (10 samples) | Storage |
|----------|------|--------|-------------------|---------|
| BAM RNA-seq | 20 | 64 GB | 4-6 hours | 50 GB |
| BAM WES | 40 | 128 GB | 8-12 hours | 200 GB |
| BAM WGS | 40 | 180 GB | 16-24 hours | 500 GB |
| FASTQ RNA-seq | 20 | 64 GB | 3-5 hours | 50 GB |
| FASTQ DNA | 40 | 128-180 GB | 6-24 hours | 200-500 GB |

### SLURM Partition Guide (CSC Puhti)

| Data Type | Partition | Time Limit |
|-----------|-----------|------------|
| RNA-seq (10 samples) | `small` | 8 hours |
| WES (10 samples) | `small` | 16 hours |
| WGS (10 samples) | `small` | 48 hours |
| WGS (>10 samples) | `small` | 72 hours |

---

## Sample Command Reference

### BAM Input

```bash
# RNA-seq
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples.csv \
    --tools arcashla,optitype,spechla \
    --seq_type rna \
    -profile singularity

# WES
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples.csv \
    --tools spechla,hlahd,arcashla,optitype \
    --seq_type dna \
    -profile singularity

# WGS
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples.csv \
    --tools spechla,hlahd,hlala,xhla \
    --reference hg38 \
    -profile singularity
```

### FASTQ Input

```bash
# RNA-seq FASTQ
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_fastq.csv \
    --tools arcashla,optitype,spechla \
    --seq_type rna \
    -profile singularity

# DNA FASTQ (WES/WGS)
nextflow run hla_typing_pipeline/main.nf \
    --input_samplesheet samples_fastq.csv \
    --tools spechla,hlahd,arcashla,optitype \
    --seq_type dna \
    -profile singularity
```

---

## Samplesheet Formats

### BAM Samplesheet

```csv
sample_id,bam_path
Sample_001,/path/to/Sample_001.bam
Sample_002,/path/to/Sample_002.bam
```

### FASTQ Samplesheet

```csv
sample_id,fastq_1,fastq_2
Sample_001,/path/to/Sample_001_R1.fastq.gz,/path/to/Sample_001_R2.fastq.gz
Sample_002,/path/to/Sample_002_R1.fastq.gz,/path/to/Sample_002_R2.fastq.gz
```

---

## Expected Output Quality

### Confidence Score Expectations

| Data Type | Expected Confidence |
|-----------|-------------------|
| WGS 30x | 0.90-1.00 |
| WES 100x | 0.80-0.95 |
| RNA-seq (high expr) | 0.85-1.00 |
| RNA-seq (low expr) | 0.50-0.80 |

### Tool Agreement Expectations

| Scenario | Expected Agreement |
|----------|-------------------|
| WGS, 4 tools | 4/4 for most genes |
| WES, 4 tools | 3-4/4 typical |
| RNA-seq, 3 tools | 2-3/3 typical |

---

## Scenario Documents

| Scenario | Document | Best For |
|----------|----------|----------|
| 1 | [SCENARIO_1_BAM_RNASEQ.md](SCENARIO_1_BAM_RNASEQ.md) | RNA-seq BAM files |
| 2 | [SCENARIO_2_BAM_WES.md](SCENARIO_2_BAM_WES.md) | Exome BAM files |
| 3 | [SCENARIO_3_BAM_WGS.md](SCENARIO_3_BAM_WGS.md) | Genome BAM files |
| 4 | [SCENARIO_4_FASTQ_RNASEQ.md](SCENARIO_4_FASTQ_RNASEQ.md) | RNA-seq FASTQ files |
| 5 | [SCENARIO_5_FASTQ_DNA.md](SCENARIO_5_FASTQ_DNA.md) | WES/WGS FASTQ files |

---

## Quick Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|--------------|----------|
| Low confidence Class II | RNA-seq tissue | Normal, use Class I |
| HLA*LA skipped | FASTQ input | Expected, use BAM for HLA*LA |
| xHLA fails | hg38 coordinates | Check BAM reference version |
| Very slow | Large WGS files | Increase time/memory |
| Out of memory | Many samples | Process in batches |
| Tool disagreement | Low coverage | Trust tool with most reads |
