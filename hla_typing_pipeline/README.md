# HLA Typing Pipeline

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/Singularity-%E2%89%A53.0-blue.svg)](https://sylabs.io/singularity/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive multi-tool HLA typing pipeline with **weighted consensus voting** for reliable allele determination from WGS, WES, and RNA-seq data.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Input Formats](#input-formats)
- [Output Files](#output-files)
- [HLA Tools](#hla-tools)
- [Parameters](#parameters)
- [Consensus Algorithm](#consensus-algorithm)
- [LOH Analysis](#loh-analysis)
- [CSC Puhti Guide](#csc-puhti-guide)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Overview

Human Leukocyte Antigen (HLA) typing is critical for:
- **Transplantation matching** - Donor-recipient compatibility
- **Pharmacogenomics** - Drug hypersensitivity prediction (e.g., HLA-B*57:01 and abacavir)
- **Disease association studies** - Autoimmune disease risk assessment
- **Cancer immunotherapy** - Neoantigen prediction and vaccine design

This pipeline addresses the challenge of **inconsistent results** across different HLA typing tools by integrating **6 state-of-the-art tools** and applying a **weighted consensus voting algorithm** to produce reliable, high-confidence HLA calls.

### Pipeline Workflow

```
                                    ┌─────────────┐
                                    │   SpecHLA   │──┐
                                    └─────────────┘  │
┌───────────┐     ┌──────────┐     ┌─────────────┐  │     ┌─────────────┐     ┌─────────────┐
│  BAM or   │────▶│    QC    │────▶│   HLA-HD    │──┼────▶│  Consensus  │────▶│   Results   │
│  FASTQ    │     │  Check   │     └─────────────┘  │     │   Voting    │     │  & Reports  │
└───────────┘     └──────────┘     ┌─────────────┐  │     └─────────────┘     └─────────────┘
                                    │   HLA*LA    │──┤
                                    └─────────────┘  │
                                    ┌─────────────┐  │
                                    │  arcasHLA   │──┤
                                    └─────────────┘  │
                                    ┌─────────────┐  │
                                    │  OptiType   │──┤
                                    └─────────────┘  │
                                    ┌─────────────┐  │
                                    │    xHLA     │──┘
                                    └─────────────┘
```

---

## Features

| Feature | Description |
|---------|-------------|
| **Multi-tool Integration** | 6 HLA typing tools for comprehensive coverage |
| **Weighted Consensus** | Confidence-based voting for reliable calls |
| **Dual Input Support** | BAM files or paired-end FASTQ files |
| **Quality Control** | Automatic QC with warnings for low-quality samples |
| **Visual Reports** | Interactive HTML reports with charts and heatmaps |
| **LOH Detection** | Loss of Heterozygosity analysis for tumor samples |
| **HPC Ready** | SLURM cluster and CSC Puhti support |
| **Containerized** | Singularity/Docker for reproducibility |

---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/uoozcan/local_pipeline.git
cd local_pipeline/hla_typing_pipeline
```

### 2. Run with a Single BAM File

```bash
nextflow run main.nf \
    --input_bam /path/to/sample.bam \
    --outdir results \
    --tools spechla,hlahd \
    -profile singularity
```

### 3. Run with Paired FASTQ Files

```bash
nextflow run main.nf \
    --input_fastq_1 sample_R1.fastq.gz \
    --input_fastq_2 sample_R2.fastq.gz \
    --outdir results \
    --tools spechla,hlahd,arcashla \
    -profile singularity
```

### 4. View Results

```bash
# Check consensus HLA types
cat results/SampleName/SampleName_consensus.txt

# Open interactive report
firefox results/SampleName/visualizations/SampleName_report.html
```

---

## Installation

### Prerequisites

| Software | Version | Purpose |
|----------|---------|---------|
| Nextflow | ≥22.10.0 | Workflow engine |
| Singularity | ≥3.0 | Container runtime |
| Samtools | ≥1.10 | BAM processing |

### Step 1: Install Nextflow

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow ~/bin/  # or add to PATH

# Verify installation
nextflow -version
```

### Step 2: Clone Repository

```bash
git clone https://github.com/uoozcan/local_pipeline.git
cd local_pipeline
```

### Step 3: Setup Containers and Databases

#### Option A: Automated Setup (Recommended for CSC Puhti)

```bash
# Run the comprehensive setup script
./hla_typing_pipeline/scripts/puhti_full_setup.sh

# Or for a minimal SpecHLA-only setup
./hla_typing_pipeline/scripts/puhti_quick_setup.sh
```

#### Option B: Manual Setup

```bash
# Create directories
mkdir -p hla_references/containers
mkdir -p hla_references/databases

# Pull containers
singularity pull hla_references/containers/hlahd.sif docker://quay.io/biocontainers/hlahd:1.7.0
singularity pull hla_references/containers/arcashla.sif docker://quay.io/biocontainers/arcas-hla:0.5.0
singularity pull hla_references/containers/optitype.sif docker://fred2/optitype:latest
singularity pull hla_references/containers/xhla.sif docker://humanlongevity/hla:latest

# Download HLA-HD database (requires registration)
# Visit: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
```

### Step 4: Configure User Settings

```bash
# Copy template
cp hla_typing_pipeline/conf/user.config.template hla_typing_pipeline/conf/user.config

# Edit with your paths
nano hla_typing_pipeline/conf/user.config
```

**user.config:**
```groovy
params {
    container_dir = '/path/to/hla_references/containers'
    hlahd_db = '/path/to/hla_references/databases/hlahd_db'
    spechla_path = '/path/to/spechla_local'
}
```

---

## Usage

### Single Sample Analysis

#### From BAM File
```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --outdir results \
    --tools spechla,hlahd,arcashla \
    -profile singularity
```

#### From FASTQ Files
```bash
nextflow run main.nf \
    --input_fastq_1 sample_R1.fastq.gz \
    --input_fastq_2 sample_R2.fastq.gz \
    --outdir results \
    --tools spechla,hlahd \
    -profile singularity
```

### Batch Processing (Multiple Samples)

#### Create a Samplesheet

**For BAM files (`samples.csv`):**
```csv
sample_id,bam_path
Patient001,/data/bams/Patient001.bam
Patient002,/data/bams/Patient002.bam
Patient003,/data/bams/Patient003.bam
```

**For FASTQ files (`samples.csv`):**
```csv
sample_id,fastq_1,fastq_2
Patient001,/data/fastq/Patient001_R1.fq.gz,/data/fastq/Patient001_R2.fq.gz
Patient002,/data/fastq/Patient002_R1.fq.gz,/data/fastq/Patient002_R2.fq.gz
```

#### Run Batch Analysis
```bash
nextflow run main.nf \
    --input_samplesheet samples.csv \
    --outdir results \
    --tools spechla,hlahd,arcashla,optitype \
    -profile singularity
```

### Running on SLURM Cluster

```bash
nextflow run main.nf \
    --input_samplesheet samples.csv \
    --outdir results \
    -profile singularity,slurm
```

---

## Input Formats

### BAM Files

- Aligned to **hg38** or **hg19** reference genome
- Index file (`.bai`) recommended (auto-generated if missing)
- Can be from WGS, WES, or RNA-seq

### FASTQ Files

- Paired-end reads (R1 and R2)
- Can be gzipped (`.fastq.gz`) or uncompressed
- Minimum recommended: 50bp read length

### Samplesheet Format

The pipeline auto-detects input type from the CSV header:

| Header | Input Type |
|--------|------------|
| `sample_id,bam_path` | BAM input |
| `sample_id,fastq_1,fastq_2` | FASTQ input |

---

## Output Files

```
results/
├── SampleName/
│   ├── qc/
│   │   └── SampleName_qc_report.txt        # Quality metrics & warnings
│   ├── fastqc/
│   │   └── SampleName_fastqc.html          # Read quality report
│   ├── spechla/
│   │   └── SampleName_spechla.txt          # SpecHLA results
│   ├── hlahd/
│   │   └── SampleName_hlahd.txt            # HLA-HD results
│   ├── arcashla/
│   │   └── SampleName_arcashla.txt         # arcasHLA results
│   ├── optitype/
│   │   └── SampleName_optitype.txt         # OptiType results
│   ├── xhla/
│   │   └── SampleName_xhla.txt             # xHLA results
│   ├── visualizations/
│   │   ├── SampleName_confidence.png       # Confidence chart
│   │   ├── SampleName_coverage.png         # Read coverage plot
│   │   ├── SampleName_agreement.png        # Tool agreement heatmap
│   │   └── SampleName_report.html          # Interactive HTML report
│   ├── loh/                                 # (if LOH analysis enabled)
│   │   ├── SampleName_hla_loh.txt          # LOH results
│   │   └── SampleName_loh_plot.png         # LOH visualization
│   ├── SampleName_consensus.txt            # Final HLA calls
│   └── SampleName_comparison.txt           # Tool comparison matrix
├── summary/
│   ├── hla_summary_report.html             # Multi-sample summary
│   ├── hla_summary_statistics.tsv          # Statistics table
│   └── allele_frequency.png                # Population frequencies
├── multiqc/
│   └── multiqc_report.html                 # Aggregated QC report
└── pipeline_info/
    ├── timeline.html                        # Execution timeline
    └── report.html                          # Pipeline report
```

### Consensus Output Format

```
# HLA Consensus Results for Patient001
# Tools used: spechla, hlahd, arcashla
# Resolution: 2-field
# Weighting: read_confidence

Gene    Allele1    Allele2    Confidence    Reads1    Reads2
HLA-A   A*02:01    A*03:01    0.95          1250      980
HLA-B   B*07:02    B*44:02    0.88          850       720
HLA-C   C*07:01    C*05:01    1.00          1100      1050
HLA-DRB1 DRB1*15:01 DRB1*04:01 0.92         890       760
HLA-DQB1 DQB1*06:02 DQB1*03:02 0.85         650       580
HLA-DPB1 DPB1*04:01 DPB1*02:01 0.78         420       390
```

---

## HLA Tools

### Tool Comparison

| Tool | Algorithm | Class I | Class II | Input | Best For |
|------|-----------|---------|----------|-------|----------|
| **SpecHLA** | Assembly + Alignment | ✓ | ✓ | BAM/FASTQ | High resolution, phasing |
| **HLA-HD** | Exhaustive alignment | ✓ | ✓ | BAM/FASTQ | Accuracy, rare alleles |
| **HLA*LA** | Graph-based | ✓ | ✓ | BAM only | WGS data |
| **arcasHLA** | Kallisto pseudoalign | ✓ | ✓ | BAM/FASTQ | RNA-seq, speed |
| **OptiType** | Integer programming | ✓ | ✗ | BAM/FASTQ | Class I accuracy |
| **xHLA** | K-mer based | ✓ | ✓ | BAM (rec.) | Speed, WGS |

### Tool Selection Recommendations

| Data Type | Recommended Tools |
|-----------|-------------------|
| **WGS** | `spechla,hlahd,hlala` |
| **WES** | `spechla,hlahd,arcashla` |
| **RNA-seq** | `arcashla,optitype,spechla` |
| **Quick analysis** | `arcashla,optitype` |
| **High accuracy** | `spechla,hlahd,hlala,arcashla` |

### BAM Input: Multi-Tool Processing

When running multiple tools from a single BAM file, each tool uses its own **specialized preprocessing**:

```bash
# Run multiple tools from BAM
nextflow run main.nf \
    --input_bam sample.bam \
    --tools optitype,arcashla,spechla \
    -profile singularity
```

**How each tool handles BAM input:**

| Tool | Preprocessing | Description |
|------|--------------|-------------|
| **SpecHLA** | HLA region extraction → FASTQ | Extracts chr6:28510120-33480577 (hg38) or chr6:28477797-33448354 (hg19), converts to paired FASTQ |
| **arcasHLA** | `arcasHLA extract` | Uses its own HLA reference index to extract relevant reads |
| **OptiType** | `samtools sort -n` + `samtools fastq` | Standard name-sorted BAM to paired FASTQ conversion |
| **HLA-HD** | Chromosome 6 extraction → FASTQ | Extracts MHC region reads |
| **HLA\*LA** | Direct BAM processing | Works directly with BAM, no conversion needed |
| **xHLA** | HLA region extraction | Extracts HLA region, processes internally |

**Why specialized preprocessing?**
- Each tool is optimized for specific read subsets
- Generic conversion would include unnecessary reads
- Tools like arcasHLA have custom reference indices
- HLA*LA uses graph-based alignment directly on BAM

All tools run **in parallel**, maximizing efficiency while ensuring each tool gets optimal input.

---

## Parameters

### Input/Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_bam` | - | Single BAM file path |
| `--input_fastq_1` | - | R1 FASTQ file path |
| `--input_fastq_2` | - | R2 FASTQ file path |
| `--input_samplesheet` | - | CSV samplesheet path |
| `--outdir` | `./results` | Output directory |

### Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tools` | `spechla,hlahd` | Comma-separated tool list |
| `--reference` | `hg38` | Reference genome (hg38/hg19) |
| `--resolution` | `2-field` | Output resolution (2-field/4-field) |
| `--min_tools` | `1` | Minimum tools for consensus |
| `--weighting` | `read_confidence` | Voting weight method |

### Quality Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_hla_reads` | `1000` | Warning threshold for HLA reads |
| `--min_read_length` | `50` | Minimum average read length |
| `--expected_reads` | `1000` | Expected reads for confidence calc |

### LOH Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_loh` | `false` | Enable LOH detection |
| `--tumor_purity` | - | Tumor purity (0-1) |
| `--tumor_ploidy` | - | Tumor ploidy estimate |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `8` | Maximum CPUs per process |
| `--max_memory` | `32.GB` | Maximum memory per process |

---

## Consensus Algorithm

The pipeline uses a **weighted voting algorithm** to determine consensus HLA calls:

### How It Works

1. **Collect Results**: Gather allele calls from each tool
2. **Normalize**: Convert to specified resolution (2-field or 4-field)
3. **Calculate Weights**: Based on read support confidence
4. **Vote**: Sum weighted votes for each allele
5. **Select**: Choose top 2 alleles per gene
6. **Score**: Calculate final confidence score

### Weighting Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| `equal` | All tools weight = 1.0 | Equal trust in all tools |
| `read_confidence` | Weight = min(1.0, reads/expected) | **Default** - considers read support |
| `tool_quality` | Predefined tool accuracy weights | Established benchmarks |

### Confidence Interpretation

| Score | Interpretation | Action |
|-------|----------------|--------|
| **0.90-1.00** | High confidence | Reliable call |
| **0.75-0.89** | Good confidence | Generally reliable |
| **0.50-0.74** | Moderate | Consider manual review |
| **<0.50** | Low confidence | Manual review recommended |

---

## LOH Analysis

Loss of Heterozygosity (LOH) analysis detects HLA allele loss in tumor samples.

### When to Use

- Tumor samples with known purity estimates
- Investigating immune evasion mechanisms
- Patient stratification for immunotherapy

### Running LOH Analysis

```bash
nextflow run main.nf \
    --input_bam tumor.bam \
    --tools spechla,hlahd \
    --run_loh true \
    --tumor_purity 0.75 \
    --tumor_ploidy 2.1 \
    -profile singularity
```

### LOH Output

```
Sample    Gene    Allele1    Allele2    CopyRatio  LOH
Patient1  HLA-A   A*02:01    A*03:01    2:0        Y
Patient1  HLA-B   B*07:02    B*44:02    1:1        N
Patient1  HLA-C   C*07:01    C*05:01    1:1        N
```

---

## CSC Puhti Guide

For detailed instructions on running on CSC Puhti supercomputer, see [docs/PUHTI_GUIDE.md](docs/PUHTI_GUIDE.md).

### Quick Setup on Puhti

```bash
# Login to Puhti
ssh username@puhti.csc.fi

# Navigate to scratch space
cd /scratch/project_XXXXXXX/$USER

# Clone and setup
git clone https://github.com/uoozcan/local_pipeline.git hla_analysis
cd hla_analysis
./hla_typing_pipeline/scripts/puhti_full_setup.sh

# Submit a job
sbatch hla_typing_pipeline/scripts/submit_hla_batch.sh
```

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| **Low HLA reads warning** | Increase sequencing depth or use HLA-enriched data |
| **Container not found** | Check `container_dir` path in config |
| **Out of memory** | Increase `--max_memory` parameter |
| **BAM index missing** | Run `samtools index sample.bam` |
| **Tool timeout** | Increase time limits in `nextflow.config` |

### Getting Help

```bash
# View pipeline help
nextflow run main.nf --help

# Check Nextflow logs
cat .nextflow.log

# View process logs
cat work/xx/xxxxxx/.command.log
```

### Resume Failed Runs

```bash
# Resume from last checkpoint
nextflow run main.nf [options] -resume
```

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **SpecHLA**: Liu et al. (2023) *Bioinformatics*
- **HLA-HD**: Kawaguchi et al. (2017) *Human Mutation*
- **HLA*LA**: Dilthey et al. (2019) *Bioinformatics*
- **arcasHLA**: Orenbuch et al. (2020) *Bioinformatics*
- **OptiType**: Szolek et al. (2014) *Bioinformatics*
- **xHLA**: Xie et al. (2017) *Genome Medicine*

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

- **Issues & Feature Requests**: [GitHub Issues](https://github.com/uoozcan/local_pipeline/issues)
- **Author**: Umut Ozcan

---

<p align="center">
  <i>Developed for reliable HLA typing in clinical and research applications</i>
</p>
