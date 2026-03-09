# Supplementary Documentation: HLA Typing Pipeline

## Technical Reference for High-Performance Computing Environments

**Pipeline Version:** 1.2.0
**Documentation Version:** 1.0
**Last Updated:** February 2026

---

## Table of Contents

1. [Pipeline Overview and Architecture](#1-pipeline-overview-and-architecture)
2. [System Requirements](#2-system-requirements)
3. [Installation Prerequisites](#3-installation-prerequisites)
4. [Step-by-Step Installation Instructions](#4-step-by-step-installation-instructions)
5. [Configuration Guidelines](#5-configuration-guidelines)
6. [Systematic Troubleshooting](#6-systematic-troubleshooting)
7. [Frequently Asked Questions](#7-frequently-asked-questions)
8. [Performance Optimization](#8-performance-optimization)
9. [Testing and Validation](#9-testing-and-validation)
10. [Support and Contact Information](#10-support-and-contact-information)

---

## 1. Pipeline Overview and Architecture

### 1.1 Purpose

The HLA Typing Pipeline is a Nextflow-based bioinformatics workflow designed for high-resolution Human Leukocyte Antigen (HLA) typing from next-generation sequencing data. The pipeline integrates multiple state-of-the-art HLA typing tools and employs a weighted consensus voting algorithm to produce robust and reliable HLA genotypes.

### 1.2 Key Features

- **Multi-tool integration**: Supports six HLA typing algorithms (SpecHLA, HLA-HD, HLA*LA, arcasHLA, OptiType, xHLA)
- **Flexible input**: Accepts both BAM and paired-end FASTQ files
- **Batch processing**: Process multiple samples via CSV samplesheets
- **Weighted consensus**: Combines results using read-confidence or tool-quality weighting
- **Comprehensive QC**: Automated quality control with configurable thresholds
- **Loss of Heterozygosity (LOH) analysis**: Optional tumor HLA-LOH detection
- **Visualization**: Per-sample and cohort-level visualizations
- **HPC-ready**: Native support for SLURM, PBS, and SGE schedulers

### 1.3 Workflow Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           INPUT DATA                                        │
│                    BAM files / Paired FASTQ files                           │
│                         (Single or Batch)                                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        QUALITY CONTROL                                      │
│   ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                 │
│   │   QC_BAM /   │───▶│   FastQC     │───▶│  QC Report   │                 │
│   │  QC_FASTQ    │    │   Analysis   │    │  Generation  │                 │
│   └──────────────┘    └──────────────┘    └──────────────┘                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       HLA TYPING TOOLS                                      │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐  │
│  │SpecHLA │ │ HLA-HD  │ │ HLA*LA  │ │arcasHLA│ │OptiType │ │  xHLA   │  │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘  │
│       │           │           │           │           │           │        │
│       └───────────┴───────────┴─────┬─────┴───────────┴───────────┘        │
└─────────────────────────────────────┼───────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    CONSENSUS & REPORTING                                    │
│   ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                 │
│   │  Weighted    │───▶│ Visualization│───▶│   MultiQC    │                 │
│   │  Consensus   │    │   Reports    │    │   Summary    │                 │
│   └──────────────┘    └──────────────┘    └──────────────┘                 │
│                              │                                              │
│                              ▼                                              │
│                    ┌──────────────┐                                         │
│                    │  LOH Analysis│ (Optional)                              │
│                    └──────────────┘                                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          OUTPUT                                             │
│   • Per-sample consensus HLA types (2-field or 4-field resolution)         │
│   • Tool comparison matrices                                                │
│   • QC reports and visualizations                                           │
│   • MultiQC aggregated report                                               │
│   • Pipeline execution reports (timeline, trace, DAG)                       │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 1.4 Supported HLA Typing Tools

| Tool | Version | Input Support | HLA Genes | Reference |
|------|---------|---------------|-----------|-----------|
| SpecHLA | 1.0.7 | BAM, FASTQ | Class I & II | Bai et al., 2024 |
| HLA-HD | 1.4.0+ | BAM, FASTQ | Class I & II | Kawaguchi et al., 2017 |
| HLA*LA | 1.0.3 | BAM only | Class I & II | Dilthey et al., 2019 |
| arcasHLA | 0.5.0 | BAM, FASTQ | Class I & II | Orenbuch et al., 2020 |
| OptiType | 1.3.5 | BAM, FASTQ | Class I only | Szolek et al., 2014 |
| xHLA | 1.0 | BAM, FASTQ | Class I & II | Xie et al., 2017 |

### 1.5 Default HLA Genes Typed

The pipeline types the following classical HLA genes by default:
- **Class I**: HLA-A, HLA-B, HLA-C
- **Class II**: HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1, HLA-DPB1

---

## 2. System Requirements

### 2.1 Minimum Hardware Specifications

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU Cores | 8 cores | 16+ cores |
| RAM | 32 GB | 64 GB |
| Storage | 100 GB free space | 500 GB+ SSD |
| Network | Standard | High-bandwidth for container pulls |

### 2.2 Per-Process Resource Requirements

| Process | CPUs | Memory | Walltime |
|---------|------|--------|----------|
| QC_BAM / QC_FASTQ | 2 | 4 GB | 1 hour |
| SpecHLA | 8 | 16 GB | 4 hours |
| HLA-HD | 8 | 16 GB | 4 hours |
| HLA*LA | 8 | 20 GB | 6 hours |
| arcasHLA | 8 | 8 GB | 2 hours |
| OptiType | 4 | 8 GB | 2 hours |
| xHLA | 4 | 8 GB | 2 hours |
| Consensus | 1 | 2 GB | 30 min |
| Visualization | 1 | 4 GB | 30 min |
| MultiQC | 1 | 4 GB | 30 min |

### 2.3 Supported Operating Systems

- **Linux**: Ubuntu 18.04+, CentOS 7+, RHEL 7+, Debian 10+
- **HPC Distributions**: Scientific Linux, Rocky Linux

> **Note**: macOS and Windows are not officially supported for production use. Development and testing on these platforms should use Docker containerization.

### 2.4 HPC Scheduler Compatibility

| Scheduler | Profile | Tested Versions |
|-----------|---------|-----------------|
| SLURM | `slurm` | 18.08+, 20.x, 21.x, 22.x |
| PBS/Torque | `pbs` | OpenPBS 19+, PBS Pro 2021+ |
| SGE/UGE | `sge` | SGE 8.1+, UGE 8.6+ |
| LSF | `lsf` | 10.1+ |
| Local | `local` | N/A (for development) |

---

## 3. Installation Prerequisites

### 3.1 Required Software

#### 3.1.1 Nextflow

**Required Version**: >= 22.10.0

```bash
# Installation via curl (recommended)
curl -s https://get.nextflow.io | bash

# Move to a directory in your PATH
mkdir -p $HOME/bin
mv nextflow $HOME/bin/

# Verify installation
nextflow -version
```

**Expected output**:
```
      N E X T F L O W
      version 23.10.0 build 5889
      created 15-10-2023 15:07 UTC
```

#### 3.1.2 Java Runtime Environment

**Required Version**: Java 11 or later (Java 17 recommended)

```bash
# Check Java version
java -version

# Installation on Ubuntu/Debian
sudo apt-get update
sudo apt-get install openjdk-17-jre

# Installation on CentOS/RHEL
sudo yum install java-17-openjdk

# Installation via module (HPC)
module load java/17
```

**Expected output**:
```
openjdk version "17.0.6" 2023-01-17
OpenJDK Runtime Environment (build 17.0.6+10-Ubuntu-0ubuntu122.04)
OpenJDK 64-Bit Server VM (build 17.0.6+10-Ubuntu-0ubuntu122.04, mixed mode, sharing)
```

#### 3.1.3 Container Technology

**Option A: Singularity/Apptainer (Recommended for HPC)**

**Required Version**: Singularity >= 3.7.0 or Apptainer >= 1.0.0

```bash
# Check Singularity version
singularity --version

# On HPC systems, typically loaded via module
module load singularity/3.8.0

# Or with Apptainer
module load apptainer/1.1.0
```

**Option B: Docker (Development/Cloud)**

**Required Version**: Docker >= 20.10

```bash
# Check Docker version
docker --version

# Verify Docker daemon is running
docker info
```

### 3.2 Container Images

The pipeline requires the following Singularity containers:

| Container | Description | Size (approx.) |
|-----------|-------------|----------------|
| `basetools.sif` | samtools, bc, fastqc, python | 500 MB |
| `hla_postprocess.sif` | matplotlib, pandas, seaborn | 800 MB |
| `spechla_with_spechap.sif` | SpecHLA with SpecHap | 2.5 GB |
| `hlahd.sif` | HLA-HD | 1.2 GB |
| `hlala.sif` | HLA*LA | 1.8 GB |
| `arcashla.sif` | arcasHLA | 1.0 GB |
| `optitype.sif` | OptiType | 1.5 GB |
| `xhla.sif` | xHLA | 800 MB |
| `flow_optitype.sif` | Flow-OptiType (optional) | 1.5 GB |
| `bamqc.sif` | BAMQC (optional) | 600 MB |

### 3.3 Reference Databases

The following databases are required:

| Database | Tool | Size | Location |
|----------|------|------|----------|
| HLA-HD database | HLA-HD | ~500 MB | `hla_references/databases/hlahd_db/` |
| HLA*LA graphs | HLA*LA | ~2 GB | `hla_references/databases/hlala_graphs/` |
| IPD-IMGT/HLA | Multiple | Bundled | Within containers |

### 3.4 Directory Structure

Recommended directory structure for deployment:

```
hla_typing_pipeline/
├── main.nf                     # Main workflow
├── nextflow.config             # Configuration
├── modules/                    # Process definitions
│   ├── spechla.nf
│   ├── hlahd.nf
│   ├── hlala.nf
│   ├── arcashla.nf
│   ├── optitype.nf
│   ├── xhla.nf
│   ├── qc.nf
│   ├── consensus.nf
│   ├── fastqc.nf
│   ├── multiqc.nf
│   ├── visualize.nf
│   └── loh.nf
├── assets/
│   ├── multiqc_config.yaml
│   ├── samplesheet_template.csv
│   └── samplesheet_fastq_template.csv
├── bin/                        # Custom scripts
│   └── consensus_voting.py
└── ../hla_references/          # External references
    ├── containers/
    │   ├── basetools.sif
    │   ├── hla_postprocess.sif
    │   ├── spechla_with_spechap.sif
    │   ├── hlahd.sif
    │   ├── hlala.sif
    │   ├── arcashla.sif
    │   ├── optitype.sif
    │   └── xhla.sif
    └── databases/
        ├── hlahd_db/
        └── hlala_graphs/
```

---

## 4. Step-by-Step Installation Instructions

### 4.1 Method A: Container-Based Deployment (Recommended)

This method uses pre-built Singularity containers and is recommended for HPC environments.

#### Step 1: Clone or Download the Pipeline

```bash
# Clone from repository (if available)
git clone https://github.com/your-org/hla-typing-pipeline.git
cd hla-typing-pipeline

# Or download and extract release
wget https://github.com/your-org/hla-typing-pipeline/releases/download/v1.2.0/hla-typing-pipeline-1.2.0.tar.gz
tar -xzf hla-typing-pipeline-1.2.0.tar.gz
cd hla-typing-pipeline-1.2.0
```

#### Step 2: Set Up Reference Directory

```bash
# Create reference directory structure
mkdir -p ../hla_references/containers
mkdir -p ../hla_references/databases/hlahd_db
mkdir -p ../hla_references/databases/hlala_graphs
```

#### Step 3: Download or Build Container Images

**Option A: Download Pre-built Containers**

```bash
# Download containers from repository
cd ../hla_references/containers

# Example download commands (adjust URLs as needed)
wget https://your-container-registry/basetools.sif
wget https://your-container-registry/hla_postprocess.sif
wget https://your-container-registry/spechla_with_spechap.sif
wget https://your-container-registry/hlahd.sif
wget https://your-container-registry/hlala.sif
wget https://your-container-registry/arcashla.sif
wget https://your-container-registry/optitype.sif
wget https://your-container-registry/xhla.sif

cd ../../hla_typing_pipeline
```

**Option B: Build Containers from Definition Files**

```bash
# Build each container (requires root or fakeroot)
cd container_definitions/

singularity build --fakeroot basetools.sif basetools.def
singularity build --fakeroot hla_postprocess.sif hla_postprocess.def
singularity build --fakeroot spechla_with_spechap.sif spechla.def
# ... repeat for other containers

mv *.sif ../hla_references/containers/
```

#### Step 4: Set Up HLA Databases

**HLA-HD Database**:
```bash
# Download HLA-HD database
cd ../hla_references/databases/hlahd_db

# Follow HLA-HD documentation to download/setup database
# Typically involves downloading from IPD-IMGT/HLA
wget ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hlahd/hlahd_database.tar.gz
tar -xzf hlahd_database.tar.gz
```

**HLA*LA Graph Database**:
```bash
# Download HLA*LA graph
cd ../hlala_graphs

# Download appropriate graph for your reference genome
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xzf PRG_MHC_GRCh38_withIMGT.tar.gz
```

#### Step 5: Configure the Pipeline

Edit `nextflow.config` to set correct paths:

```groovy
params {
    // Container paths
    container_dir = "/path/to/hla_references/containers"
    basetools_container = "${container_dir}/basetools.sif"
    postprocess_container = "${container_dir}/hla_postprocess.sif"

    // Database paths
    hlahd_db = "/path/to/hla_references/databases/hlahd_db"
    hlala_graphs = "/path/to/hla_references/databases/hlala_graphs"
}
```

#### Step 6: Verify Installation

```bash
# Run with test profile
nextflow run main.nf -profile test,singularity --help

# Check that containers can be pulled/executed
singularity exec ../hla_references/containers/basetools.sif samtools --version
```

### 4.2 Method B: Local Installation (Advanced)

For environments where containers are not available.

> **Warning**: Local installation requires manual management of all dependencies and is more prone to version conflicts.

#### Step 1: Install System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-pip \
    samtools \
    bwa

# CentOS/RHEL
sudo yum groupinstall -y "Development Tools"
sudo yum install -y \
    zlib-devel \
    bzip2-devel \
    xz-devel \
    ncurses-devel \
    libcurl-devel \
    python3 \
    python3-pip
```

#### Step 2: Install Python Dependencies

```bash
pip3 install --user \
    pandas>=1.3.0 \
    numpy>=1.21.0 \
    matplotlib>=3.4.0 \
    seaborn>=0.11.0 \
    biopython>=1.79
```

#### Step 3: Install HLA Typing Tools

Each tool must be installed according to its documentation:

- **SpecHLA**: https://github.com/deepomicslab/SpecHLA
- **HLA-HD**: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
- **HLA*LA**: https://github.com/DiltheyLab/HLA-LA
- **arcasHLA**: https://github.com/RabadanLab/arcasHLA
- **OptiType**: https://github.com/FRED-2/OptiType
- **xHLA**: https://github.com/humanlongevity/HLA

#### Step 4: Configure for Local Execution

```groovy
// In nextflow.config, add or modify:
params {
    use_local_spechla = true
    spechla_path = "/path/to/SpecHLA"
    local_lib_path = "/path/to/local/lib"  // For htslib compatibility
}

profiles {
    local {
        process.executor = 'local'
        process {
            withName: '.*' {
                container = null
            }
        }
    }
}
```

### 4.3 HPC-Specific Configuration

#### 4.3.1 SLURM Configuration

Create or modify a SLURM profile in `nextflow.config`:

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.clusterOptions = '--account=your_project_account'

        // Optional: specify partition based on resource requirements
        process {
            withLabel: 'process_high' {
                queue = 'bigmem'
                clusterOptions = '--account=your_project_account --partition=bigmem'
            }
        }
    }

    // Combined profile for Singularity on SLURM
    slurm_singularity {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.clusterOptions = '--account=your_project_account'

        singularity.enabled = true
        singularity.autoMounts = true
        singularity.runOptions = '--writable-tmpfs --cleanenv'
    }
}
```

#### 4.3.2 PBS/Torque Configuration

```groovy
profiles {
    pbs {
        process.executor = 'pbs'
        process.queue = 'workq'
        process.clusterOptions = '-A your_project_account'

        process {
            withLabel: 'process_high' {
                queue = 'largemem'
            }
        }
    }
}
```

#### 4.3.3 SGE Configuration

```groovy
profiles {
    sge {
        process.executor = 'sge'
        process.queue = 'all.q'
        process.penv = 'smp'
        process.clusterOptions = '-P your_project'

        process {
            withLabel: 'process_high' {
                queue = 'bigmem.q'
            }
        }
    }
}
```

### 4.4 Environment Module Setup

For HPC systems using environment modules:

```bash
# Create a module file: /path/to/modules/hla-pipeline/1.2.0

#%Module1.0
proc ModulesHelp { } {
    puts stderr "HLA Typing Pipeline v1.2.0"
}

module-whatis "HLA Typing Pipeline for NGS data"

# Prerequisites
prereq java/17
prereq singularity/3.8

# Set environment variables
setenv HLA_PIPELINE_HOME /path/to/hla_typing_pipeline
setenv NXF_SINGULARITY_CACHEDIR /path/to/singularity_cache

# Add to PATH
prepend-path PATH $env(HLA_PIPELINE_HOME)
prepend-path PATH $env(HLA_PIPELINE_HOME)/bin
```

Usage:
```bash
module load hla-pipeline/1.2.0
```

---

## 5. Configuration Guidelines

### 5.1 Nextflow Configuration Syntax

The pipeline uses Nextflow's configuration system. Configuration can be specified in multiple ways (in order of precedence):

1. Command-line parameters (`--param value`)
2. `-params-file params.yaml`
3. `-c custom.config`
4. `nextflow.config` in the pipeline directory

### 5.2 Essential Parameters

#### 5.2.1 Input Parameters

```groovy
params {
    // Single BAM input
    input_bam = '/path/to/sample.bam'

    // Single paired FASTQ input
    input_fastq_1 = '/path/to/sample_R1.fastq.gz'
    input_fastq_2 = '/path/to/sample_R2.fastq.gz'

    // Multiple samples via samplesheet
    input_samplesheet = '/path/to/samplesheet.csv'

    // Output directory
    outdir = './results'
}
```

#### 5.2.2 Tool Selection

```groovy
params {
    // Comma-separated list of tools to run
    // Available: spechla, hlahd, hlala, arcashla, optitype, xhla
    tools = 'spechla,hlahd'

    // Note: hlala only works with BAM input
    // Note: optitype only types Class I genes
}
```

#### 5.2.3 Reference Genome

```groovy
params {
    // Reference genome version: hg38 or hg19
    reference = 'hg38'

    // HLA genes to type
    hla_genes = 'HLA-A,HLA-B,HLA-C,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DPA1,HLA-DPB1'
}
```

#### 5.2.4 Quality Control Thresholds

```groovy
params {
    // Minimum HLA reads before warning
    min_hla_reads = 1000

    // Minimum average read length
    min_read_length = 50

    // Maximum acceptable error rate
    max_error_rate = 0.05
}
```

#### 5.2.5 Consensus Settings

```groovy
params {
    // Output resolution: 2-field or 4-field
    resolution = '2-field'

    // Minimum tools required for consensus call
    min_tools = 1

    // Expected reads per allele for confidence calculation
    expected_reads = 1000

    // Weighting method: equal, read_confidence, tool_quality
    weighting = 'read_confidence'
}
```

### 5.3 Resource Allocation

#### 5.3.1 Global Resource Limits

```groovy
params {
    max_cpus = 8
    max_memory = '32.GB'
    max_time = '24.h'
}
```

#### 5.3.2 Per-Process Resource Configuration

```groovy
process {
    // Default error handling
    errorStrategy = 'retry'
    maxRetries = 2

    // Resource labels
    withLabel: 'process_low' {
        cpus = 2
        memory = '4.GB'
        time = '2.h'
    }

    withLabel: 'process_medium' {
        cpus = 4
        memory = '8.GB'
        time = '4.h'
    }

    withLabel: 'process_high' {
        cpus = 8
        memory = '16.GB'
        time = '8.h'
    }

    // Tool-specific overrides
    withName: 'HLALA' {
        memory = '20.GB'  // HLA*LA needs more memory
        time = '6.h'
    }

    withName: 'SPECHLA' {
        cpus = { params.max_cpus }
        memory = '16.GB'
        time = '4.h'
    }
}
```

### 5.4 Scheduler-Specific Directives

#### 5.4.1 SLURM-Specific Options

```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    clusterOptions = '--account=project123'

    withName: 'HLALA' {
        clusterOptions = '--account=project123 --constraint=bigmem'
    }
}
```

#### 5.4.2 PBS-Specific Options

```groovy
process {
    executor = 'pbs'
    queue = 'workq'
    clusterOptions = '-A project123 -l select=1:ncpus=8:mem=32gb'
}
```

### 5.5 Container Options

#### 5.5.1 Singularity Configuration

```groovy
singularity {
    enabled = true
    autoMounts = true
    runOptions = '--writable-tmpfs --cleanenv'
    cacheDir = '/path/to/singularity_cache'

    // Bind additional paths if needed
    runOptions = '--writable-tmpfs --cleanenv --bind /scratch:/scratch'
}
```

#### 5.5.2 Docker Configuration

```groovy
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
    temp = 'auto'
}
```

### 5.6 Input Samplesheet Format

#### 5.6.1 BAM Samplesheet

Create a CSV file with the following format:

```csv
sample_id,bam_path
Sample1,/path/to/Sample1.bam
Sample2,/path/to/Sample2.bam
Sample3,/path/to/Sample3.bam
```

#### 5.6.2 FASTQ Samplesheet

```csv
sample_id,fastq_1,fastq_2
Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
Sample2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
Sample3,/path/to/Sample3_R1.fastq.gz,/path/to/Sample3_R2.fastq.gz
```

> **Note**: Column names for FASTQ files are flexible. The pipeline accepts: `fastq_1`/`fastq_2`, `fastq1`/`fastq2`, `fq1`/`fq2`, `read1`/`read2`, or `R1`/`R2`.

### 5.7 LOH Analysis Configuration

For tumor samples with Loss of Heterozygosity analysis:

```groovy
params {
    // Enable LOH analysis (requires SpecHLA)
    run_loh = true

    // Tumor purity (0-1), required for LOH
    tumor_purity = 0.7

    // Tumor ploidy, required for LOH
    tumor_ploidy = 2.0

    // Minimum heterozygous SNPs for LOH call
    loh_het_cutoff = 5
}
```

### 5.8 Custom Configuration File Example

Create a custom configuration file for your environment:

```groovy
// my_hpc.config

params {
    // Paths specific to your HPC
    container_dir = '/project/shared/containers/hla'
    hlahd_db = '/project/shared/databases/hlahd'
    hlala_graphs = '/project/shared/databases/hlala'

    // Resource limits for your cluster
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '48.h'
}

process {
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--account=my_allocation'

    // Use scratch for work directory
    scratch = '/scratch/$USER'
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/project/shared/singularity_cache'
    runOptions = '--writable-tmpfs --cleanenv --bind /scratch:/scratch'
}

// Execution reports
timeline {
    enabled = true
    file = "${params.outdir}/pipeline_info/timeline.html"
}
report {
    enabled = true
    file = "${params.outdir}/pipeline_info/report.html"
}
trace {
    enabled = true
    file = "${params.outdir}/pipeline_info/trace.txt"
}
```

Usage:
```bash
nextflow run main.nf -c my_hpc.config --input_samplesheet samples.csv --tools spechla,hlahd,arcashla
```

---

## 6. Systematic Troubleshooting

### 6.1 Installation and Dependency Errors

#### Issue 6.1.1: Nextflow Version Incompatibility

**Symptoms**:
```
ERROR ~ Unknown config attribute `nextflow.enable.dsl`
```
or
```
Script compilation error
- Invalid DSL syntax
```

**Cause**: The pipeline requires Nextflow >= 22.10.0 with DSL2 support. Older versions do not support the syntax used.

**Solution**:
```bash
# Check current version
nextflow -version

# Update Nextflow
nextflow self-update

# Or reinstall
curl -s https://get.nextflow.io | bash
```

**Prevention**: Add version check to your job submission scripts:
```bash
NXF_VER=$(nextflow -version 2>&1 | grep -oP 'version \K[0-9.]+')
if [[ $(echo "$NXF_VER < 22.10" | bc -l) -eq 1 ]]; then
    echo "ERROR: Nextflow >= 22.10.0 required"
    exit 1
fi
```

---

#### Issue 6.1.2: Java Version Error

**Symptoms**:
```
ERROR: Cannot find Java or it's a wrong version -- please make sure that Java 11 or later is installed
```
or
```
Error: A JNI error has occurred, please check your installation
```

**Cause**: Java is not installed, not in PATH, or version is too old.

**Solution**:
```bash
# Check Java version
java -version

# On HPC, load appropriate module
module avail java
module load java/17

# Set JAVA_HOME if needed
export JAVA_HOME=/path/to/java
export PATH=$JAVA_HOME/bin:$PATH
```

**Prevention**: Add to your `.bashrc` or job script:
```bash
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk
export PATH=$JAVA_HOME/bin:$PATH
```

---

#### Issue 6.1.3: Singularity Container Not Found

**Symptoms**:
```
FATAL:   Unable to handle docker://path/to/container.sif: failed to get checksum
```
or
```
ERROR ~ Error executing process > 'SPECHLA'
Caused by:
  Failed to pull singularity image
```

**Cause**: Container file does not exist at specified path, or path in configuration is incorrect.

**Solution**:
```bash
# Verify container exists
ls -la /path/to/containers/spechla_with_spechap.sif

# Check configuration path
grep container_dir nextflow.config

# Test container manually
singularity exec /path/to/containers/basetools.sif echo "Container works"
```

**Prevention**: Use absolute paths in configuration:
```groovy
params {
    container_dir = "/absolute/path/to/containers"  // Not relative paths
}
```

---

#### Issue 6.1.4: Missing Database Files

**Symptoms**:
```
Error: HLA-HD database not found at /path/to/hlahd_db
```
or process fails with missing reference errors.

**Cause**: Database directories are empty or paths are incorrectly configured.

**Solution**:
```bash
# Check database directories
ls -la /path/to/databases/hlahd_db/
ls -la /path/to/databases/hlala_graphs/

# Verify expected files exist
find /path/to/databases -name "*.fa" -o -name "*.idx"

# Re-download if necessary
# Follow database setup instructions in Section 4.1
```

**Prevention**: Validate database integrity after download:
```bash
# For HLA-HD database
if [ ! -f "$HLAHD_DB/dictionary/dictionary.txt" ]; then
    echo "ERROR: HLA-HD database incomplete"
    exit 1
fi
```

---

### 6.2 Job Submission and Scheduler Issues

#### Issue 6.2.1: SLURM Job Submission Failure

**Symptoms**:
```
ERROR ~ Error executing process > 'QC_BAM'
Caused by:
  Process `QC_BAM` terminated with an error exit status (1)
  sbatch: error: Batch job submission failed: Invalid account
```

**Cause**: SLURM account/allocation not specified or invalid.

**Solution**:
```bash
# Check available accounts
sacctmgr show associations user=$USER

# Set correct account in config
process {
    clusterOptions = '--account=valid_account_name'
}

# Or via command line
nextflow run main.nf --clusterOptions '--account=valid_account_name'
```

**Prevention**: Verify account before running:
```bash
sshare -U | grep $USER
```

---

#### Issue 6.2.2: Queue/Partition Not Found

**Symptoms**:
```
sbatch: error: Batch job submission failed: Invalid partition name specified
```

**Cause**: Specified queue/partition does not exist on the cluster.

**Solution**:
```bash
# List available partitions
sinfo -s

# Update configuration with valid partition
process {
    queue = 'normal'  # Use an existing partition name
}
```

**Prevention**: Query cluster configuration before setting up pipeline.

---

#### Issue 6.2.3: Walltime Exceeded

**Symptoms**:
```
slurmstepd: error: *** JOB 12345 ON node01 CANCELLED AT 2024-01-15T12:00:00 DUE TO TIME LIMIT ***
```

**Cause**: Process runtime exceeded the allocated walltime.

**Solution**:
```bash
# Increase time for specific process
process {
    withName: 'HLALA' {
        time = '12.h'  # Increase from default 6h
    }
}

# Or increase globally
params {
    max_time = '48.h'
}
```

**Prevention**: Monitor typical runtimes and set appropriate limits with buffer:
```groovy
process {
    withName: 'SPECHLA' {
        time = { 4.h * task.attempt }  // Increase on retry
    }
}
```

---

### 6.3 Memory and Resource Allocation Problems

#### Issue 6.3.1: Out of Memory (OOM) Error

**Symptoms**:
```
Command error:
  /bin/bash: line 1: 12345 Killed
```
or in SLURM logs:
```
slurmstepd: error: Detected 1 oom-kill event(s)
```

**Cause**: Process exceeded allocated memory.

**Solution**:
```groovy
// Increase memory for affected process
process {
    withName: 'HLALA' {
        memory = '32.GB'  // Increase from default 20GB
    }
}

// Enable memory scaling on retry
process {
    withName: 'HLALA' {
        memory = { 20.GB * task.attempt }
        maxRetries = 3
    }
}
```

**Prevention**: Profile memory usage on test data:
```bash
# Monitor peak memory during test run
/usr/bin/time -v nextflow run main.nf -profile test,singularity 2>&1 | grep "Maximum resident"
```

---

#### Issue 6.3.2: Insufficient CPU Cores

**Symptoms**:
Process runs very slowly, or scheduler warns about resource availability.

**Cause**: Requested CPUs exceed available resources or node limits.

**Solution**:
```groovy
// Adjust CPU allocation
params {
    max_cpus = 8  // Match cluster node configuration
}

process {
    withName: 'SPECHLA' {
        cpus = { Math.min(8, params.max_cpus) }
    }
}
```

**Prevention**: Query cluster node specifications:
```bash
sinfo -N -l | head -10
```

---

#### Issue 6.3.3: Disk Space Exhaustion

**Symptoms**:
```
ERROR: No space left on device
```
or
```
IOException: No space left on device
```

**Cause**: Work directory or output directory on filesystem with insufficient space.

**Solution**:
```bash
# Check disk usage
df -h /path/to/workdir
df -h /path/to/results

# Clean old work directories
nextflow clean -f -before 7d

# Use scratch space
nextflow run main.nf -w /scratch/$USER/nf_work --outdir /scratch/$USER/results
```

**Prevention**: Set up automatic cleanup in configuration:
```groovy
cleanup = true  // Enable automatic work directory cleanup

// Or use scratch with cleanup
process {
    scratch = '/scratch/$USER'
    afterScript = 'rm -rf $NXF_SCRATCH/*'
}
```

---

### 6.4 Container/Singularity Execution Failures

#### Issue 6.4.1: Permission Denied in Container

**Symptoms**:
```
FATAL:   container creation failed: mount hook function failure:
  mount /path source doesn't exist
```
or
```
ERROR: PermissionError: [Errno 13] Permission denied
```

**Cause**: Singularity cannot access paths outside the container due to binding restrictions.

**Solution**:
```groovy
// Add explicit binds in configuration
singularity {
    runOptions = '--writable-tmpfs --cleanenv --bind /data:/data --bind /scratch:/scratch'
}

// Or set environment variable
export SINGULARITY_BIND="/data,/scratch,/home"
```

**Prevention**: Use `autoMounts = true` and ensure data paths are under mounted directories.

---

#### Issue 6.4.2: Singularity Cache Issues

**Symptoms**:
```
FATAL:   Unable to pull docker://image: failed to get checksum for image
```
or intermittent container pull failures.

**Cause**: Singularity cache is corrupted or on slow/unreliable storage.

**Solution**:
```bash
# Clear Singularity cache
rm -rf ~/.singularity/cache/*

# Set cache to faster storage
export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

# In config:
singularity {
    cacheDir = '/scratch/$USER/singularity_cache'
}
```

**Prevention**: Pre-pull containers before running pipeline:
```bash
singularity pull /path/to/containers/image.sif docker://registry/image:tag
```

---

#### Issue 6.4.3: Incompatible Singularity Version

**Symptoms**:
```
FATAL:   could not open image /path/to/container.sif: unknown image format/type
```

**Cause**: Container built with newer Singularity version than available on system.

**Solution**:
```bash
# Check Singularity version
singularity --version

# Rebuild container with compatible version
singularity build --fakeroot new_container.sif container.def

# Or use Docker hub image directly
singularity {
    runOptions = '--writable-tmpfs'
}
process {
    withName: 'SPECHLA' {
        container = 'docker://organization/spechla:latest'
    }
}
```

**Prevention**: Document and test with specific Singularity versions. Use `apptainer` (successor to Singularity) if available.

---

### 6.5 Data Input/Output Errors

#### Issue 6.5.1: BAM File Index Missing

**Symptoms**:
```
[E::idx_find_and_load] Could not retrieve index file for 'sample.bam'
```

**Cause**: BAM index file (.bai) is missing or has incorrect name.

**Solution**:
```bash
# Create BAM index
samtools index sample.bam

# Ensure naming convention matches
# Either: sample.bam + sample.bam.bai
# Or: sample.bam + sample.bai
```

**Prevention**: The pipeline auto-generates missing indices, but pre-indexing improves performance:
```bash
for bam in *.bam; do
    if [ ! -f "${bam}.bai" ]; then
        samtools index "$bam"
    fi
done
```

---

#### Issue 6.5.2: FASTQ File Corruption

**Symptoms**:
```
ERROR: gzip: stdin: unexpected end of file
```
or
```
Error: invalid quality score character
```

**Cause**: FASTQ files are truncated, corrupted, or not properly gzipped.

**Solution**:
```bash
# Verify FASTQ integrity
gzip -t sample_R1.fastq.gz && echo "OK" || echo "CORRUPTED"

# Check FASTQ format
zcat sample_R1.fastq.gz | head -8

# Re-download or re-transfer if corrupted
```

**Prevention**: Validate files after transfer:
```bash
# Use checksums
md5sum -c checksums.md5

# Validate FASTQ format
seqkit stats sample_R1.fastq.gz sample_R2.fastq.gz
```

---

#### Issue 6.5.3: Sample ID Conflicts

**Symptoms**:
```
ERROR ~ Workflow execution completed unsuccessfully
Caused by:
  Duplicate sample ID detected: Sample1
```
or output files are overwritten.

**Cause**: Multiple samples have the same ID in the samplesheet.

**Solution**:
```bash
# Check for duplicates
cut -d',' -f1 samplesheet.csv | sort | uniq -d

# Fix samplesheet to have unique IDs
sed -i 's/Sample1/Sample1_batch2/' samplesheet.csv
```

**Prevention**: Validate samplesheet before running:
```bash
# Validation script
awk -F',' 'NR>1 {seen[$1]++} END {for(id in seen) if(seen[id]>1) print "Duplicate:", id}' samplesheet.csv
```

---

### 6.6 Pipeline Execution Interruptions

#### Issue 6.6.1: Resume Failed Pipeline

**Symptoms**:
Pipeline stopped mid-execution due to error, timeout, or system issue.

**Cause**: Various - job failure, cluster maintenance, manual interruption.

**Solution**:
```bash
# Resume from last successful checkpoint
nextflow run main.nf -resume

# With specific session ID
nextflow run main.nf -resume [session-id]

# Check session history
nextflow log
```

**Prevention**: Always use `-resume` flag for production runs:
```bash
nextflow run main.nf -resume --input_samplesheet samples.csv
```

---

#### Issue 6.6.2: Orphaned Processes After Failure

**Symptoms**:
After pipeline failure, background processes continue running and consuming resources.

**Cause**: Nextflow was terminated abruptly without cleaning up executor processes.

**Solution**:
```bash
# Find and kill orphaned processes (SLURM)
squeue -u $USER | grep nf- | awk '{print $1}' | xargs -r scancel

# Clean work directory
nextflow clean -f

# Remove lock files
rm -f .nextflow/cache/*/lock
```

**Prevention**: Use proper termination:
```bash
# Graceful shutdown
kill -TERM $(cat .nextflow.pid)
# Wait for cleanup
sleep 30
```

---

#### Issue 6.6.3: Cache Corruption

**Symptoms**:
```
ERROR ~ Failed to load cached task:
Caused by:
  java.io.EOFException
```

**Cause**: Cache files corrupted due to interrupted write operations.

**Solution**:
```bash
# Remove corrupted cache
rm -rf .nextflow/cache/*

# Start fresh without resume
nextflow run main.nf --input_samplesheet samples.csv
```

**Prevention**: Ensure work directory is on reliable filesystem with journaling.

---

## 7. Frequently Asked Questions

### 7.1 Pipeline Execution

**Q: How do I resume a failed pipeline run?**

A: Use the `-resume` flag:
```bash
nextflow run main.nf -resume --input_samplesheet samples.csv
```
Nextflow will skip completed processes and restart from the last checkpoint. The work directory must still exist with cached results.

---

**Q: How do I run the pipeline with different input data types?**

A: The pipeline supports three input methods:

1. **Single BAM file**:
   ```bash
   nextflow run main.nf --input_bam /path/to/sample.bam
   ```

2. **Single paired FASTQ files**:
   ```bash
   nextflow run main.nf --input_fastq_1 /path/to/R1.fq.gz --input_fastq_2 /path/to/R2.fq.gz
   ```

3. **Multiple samples via samplesheet**:
   ```bash
   nextflow run main.nf --input_samplesheet /path/to/samples.csv
   ```

---

**Q: How do I modify resource allocations for specific processes?**

A: Create a custom configuration file:
```groovy
// custom.config
process {
    withName: 'HLALA' {
        cpus = 12
        memory = '48.GB'
        time = '12.h'
    }
}
```
Then include it:
```bash
nextflow run main.nf -c custom.config --input_bam sample.bam
```

---

**Q: Can I run only specific HLA typing tools?**

A: Yes, use the `--tools` parameter:
```bash
# Run only SpecHLA and HLA-HD
nextflow run main.nf --tools spechla,hlahd --input_bam sample.bam

# Run all available tools
nextflow run main.nf --tools spechla,hlahd,hlala,arcashla,optitype,xhla --input_bam sample.bam
```

---

### 7.2 Output Interpretation

**Q: How do I interpret pipeline outputs and logs?**

A: The pipeline generates several output directories:

```
results/
├── <sample_id>/
│   ├── qc/
│   │   └── <sample>_qc_report.txt       # Quality control metrics
│   ├── spechla/
│   │   └── <sample>_spechla.txt         # SpecHLA results
│   ├── hlahd/
│   │   └── <sample>_hlahd.txt           # HLA-HD results
│   ├── <sample>_consensus.txt           # Final consensus HLA types
│   ├── <sample>_comparison.txt          # Tool-by-tool comparison
│   └── visualizations/
│       └── <sample>_hla_plot.png        # Visualization
├── summary/
│   ├── hla_summary.tsv                  # Multi-sample summary
│   └── allele_frequencies.tsv           # Population frequencies
├── multiqc/
│   └── multiqc_report.html              # Aggregated QC report
└── pipeline_info/
    ├── timeline_*.html                  # Execution timeline
    ├── report_*.html                    # Execution report
    └── trace_*.txt                      # Process trace
```

**Consensus file format**:
```
# HLA Consensus Results for Sample1
# Resolution: 2-field
# Tools used: spechla, hlahd
Gene    Allele1     Allele2     Confidence  Support
A       A*02:01     A*24:02     0.95        2/2
B       B*07:02     B*44:02     0.90        2/2
C       C*07:02     C*05:01     0.88        2/2
...
```

---

**Q: What do the QC warnings mean?**

A: The QC module generates the following warnings:

| Warning | Meaning | Recommendation |
|---------|---------|----------------|
| LOW_HLA_READS | Fewer than 1000 reads in HLA region | Consider deeper sequencing |
| SHORT_READS | Average read length < 50bp | May reduce typing accuracy |
| LOW_MAPQ | Average mapping quality < 20 | Check alignment quality |
| CRITICAL | < 100 HLA reads | Results likely unreliable |

---

### 7.3 Configuration and Customization

**Q: How do I change the output resolution from 2-field to 4-field?**

A: Use the `--resolution` parameter:
```bash
nextflow run main.nf --resolution '4-field' --input_bam sample.bam
```

---

**Q: How do I enable LOH analysis for tumor samples?**

A: LOH analysis requires SpecHLA and tumor purity/ploidy estimates:
```bash
nextflow run main.nf \
    --tools spechla,hlahd \
    --run_loh true \
    --tumor_purity 0.7 \
    --tumor_ploidy 2.0 \
    --input_bam tumor_sample.bam
```

---

**Q: Can I use RNA-seq data with this pipeline?**

A: Yes, the pipeline supports RNA-seq data. For OptiType, set the sequence type:
```bash
nextflow run main.nf \
    --seq_type rna \
    --tools optitype,arcashla \
    --input_bam rnaseq_sample.bam
```
Note: Some tools (especially SpecHLA and HLA-HD) work best with WGS/WES data.

---

### 7.4 Updates and Maintenance

**Q: How do I update the pipeline to a newer version?**

A:
1. **Backup current configuration**:
   ```bash
   cp nextflow.config nextflow.config.backup
   cp -r conf/ conf.backup/
   ```

2. **Update pipeline files**:
   ```bash
   git pull origin main
   # Or download new release
   ```

3. **Review changelog** for breaking changes

4. **Update containers if needed**:
   ```bash
   # Pull updated containers
   singularity pull new_spechla.sif docker://registry/spechla:v2.0
   ```

5. **Merge configuration changes**:
   ```bash
   diff nextflow.config.backup nextflow.config
   ```

---

**Q: How do I cite the pipeline and its dependencies?**

A: Please cite the following:

**Pipeline**:
> HLA Typing Pipeline v1.2.0. [URL]. Accessed [date].

**Individual tools** (cite those used):
- SpecHLA: Bai Y, et al. (2024) *Journal Name*
- HLA-HD: Kawaguchi S, et al. (2017) *Hum Mutat* 38:788-797
- HLA*LA: Dilthey AT, et al. (2019) *Bioinformatics* 35:4394-4396
- arcasHLA: Orenbuch R, et al. (2020) *Bioinformatics* 36:33-40
- OptiType: Szolek A, et al. (2014) *Bioinformatics* 30:3310-3316
- xHLA: Xie C, et al. (2017) *Genome Med* 9:53

**Nextflow**:
> Di Tommaso P, et al. (2017) *Nat Biotechnol* 35:316-319

---

## 8. Performance Optimization

### 8.1 Parallelization Strategies

#### 8.1.1 Sample-Level Parallelization

By default, the pipeline processes multiple samples in parallel. Control parallelism with:

```groovy
// Limit concurrent samples
process {
    maxForks = 10  // Maximum concurrent processes of each type
}

// Or via executor configuration
executor {
    queueSize = 50  // Maximum queued jobs
    submitRateLimit = '10/1min'  // Job submission rate
}
```

#### 8.1.2 Tool-Level Parallelization

HLA typing tools run in parallel for each sample:
```
Sample1 ──┬── SpecHLA ──┐
          ├── HLA-HD  ──┤
          ├── arcasHLA ─┼── Consensus ── Output
          └── OptiType ─┘

Sample2 ──┬── SpecHLA ──┐
          ├── HLA-HD  ──┤
          ...
```

### 8.2 Resource Tuning

#### 8.2.1 Memory Optimization

```groovy
process {
    // Use memory scaling on retry
    withName: 'HLALA' {
        memory = { 20.GB * task.attempt }
        maxRetries = 3
        errorStrategy = { task.exitStatus in [137, 140] ? 'retry' : 'finish' }
    }

    // Use local scratch to reduce memory pressure
    scratch = true
}
```

#### 8.2.2 CPU Optimization

```groovy
params {
    // Match CPUs to HPC node architecture
    max_cpus = 16  // Typical node has 16-32 cores
}

process {
    // CPU-intensive tools
    withName: 'SPECHLA|HLAHD|HLALA' {
        cpus = { Math.min(8, params.max_cpus) }
    }

    // Light processes
    withName: 'CONSENSUS|VISUALIZE' {
        cpus = 1
    }
}
```

### 8.3 I/O Optimization

#### 8.3.1 Using Local Scratch

```groovy
process {
    // Use node-local scratch for I/O-intensive operations
    scratch = '/local/scratch'

    // Alternatively, use SLURM's TMPDIR
    scratch = true  // Uses $TMPDIR
}
```

#### 8.3.2 Staging Strategies

```groovy
process {
    // Copy files to scratch (default)
    stageInMode = 'copy'
    stageOutMode = 'copy'

    // For large files, use symlinks if on shared filesystem
    stageInMode = 'symlink'
}
```

### 8.4 Cluster-Specific Optimizations

#### 8.4.1 SLURM: Using Array Jobs

For very large sample sets, consider preprocessing into chunks:
```bash
# Split samplesheet into chunks of 50
split -l 50 samples.csv samples_chunk_

# Submit each chunk
for chunk in samples_chunk_*; do
    nextflow run main.nf --input_samplesheet $chunk -profile slurm -resume
done
```

#### 8.4.2 Reducing Job Scheduler Overhead

```groovy
// Group small tasks to reduce scheduler calls
process {
    withLabel: 'process_low' {
        executor = 'local'  // Run small tasks on head node
    }
}

// Or use job grouping
executor {
    queueSize = 100
    pollInterval = '30 sec'
    submitRateLimit = '20/1min'
}
```

### 8.5 Benchmark Reference

Expected runtimes for a typical WGS sample (~30x coverage):

| Tool | Walltime | Peak Memory | CPUs |
|------|----------|-------------|------|
| QC_BAM | 5-10 min | 2 GB | 2 |
| SpecHLA | 1-2 hours | 12 GB | 8 |
| HLA-HD | 30-60 min | 10 GB | 8 |
| HLA*LA | 2-4 hours | 18 GB | 8 |
| arcasHLA | 20-40 min | 6 GB | 8 |
| OptiType | 10-30 min | 6 GB | 4 |
| xHLA | 15-30 min | 6 GB | 4 |
| Consensus | 2-5 min | 1 GB | 1 |

---

## 9. Testing and Validation

### 9.1 Test Dataset

The pipeline includes a test profile for validation:

```bash
# Run test with minimal resources
nextflow run main.nf -profile test,singularity

# Test profile settings:
# - max_cpus = 2
# - max_memory = 6.GB
# - max_time = 1.h
```

### 9.2 Installation Verification

#### 9.2.1 Verify Nextflow Installation

```bash
# Check version
nextflow -version

# Run hello world
nextflow run hello
```

#### 9.2.2 Verify Container Access

```bash
# Test each container
for sif in /path/to/containers/*.sif; do
    echo "Testing: $sif"
    singularity exec $sif echo "OK" || echo "FAILED: $sif"
done
```

#### 9.2.3 Verify Database Setup

```bash
# Check HLA-HD database
ls -la /path/to/hlahd_db/dictionary/

# Check HLA*LA graphs
ls -la /path/to/hlala_graphs/PRG_MHC_GRCh38_withIMGT/
```

### 9.3 Validation with Known Samples

For validation, use samples with known HLA types (e.g., 1000 Genomes samples or cell lines):

| Sample | Known HLA Types |
|--------|-----------------|
| NA12878 | A*01:01, A*11:01, B*08:01, B*56:01, C*01:02, C*07:01 |
| HG00096 | A*01:01, A*29:02, B*08:01, B*44:03, C*07:01, C*16:01 |

```bash
# Run validation
nextflow run main.nf \
    --input_bam NA12878.bam \
    --tools spechla,hlahd,arcashla \
    -profile singularity

# Compare results
diff results/NA12878/NA12878_consensus.txt expected_hla.txt
```

### 9.4 Expected Test Outputs

After successful test run, verify:

1. **Output directories exist**:
   ```bash
   ls results/*/
   ```

2. **Consensus files generated**:
   ```bash
   cat results/*_consensus.txt
   ```

3. **QC reports present**:
   ```bash
   cat results/*/qc/*_qc_report.txt
   ```

4. **MultiQC report generated**:
   ```bash
   ls results/multiqc/multiqc_report.html
   ```

5. **Pipeline completed successfully**:
   ```bash
   # Check exit status
   grep "success" results/pipeline_info/report_*.html
   ```

### 9.5 Troubleshooting Test Failures

If tests fail:

1. **Check Nextflow log**:
   ```bash
   cat .nextflow.log
   ```

2. **Check work directory for failed process**:
   ```bash
   cat work/xx/yyyyyy/.command.log
   cat work/xx/yyyyyy/.command.err
   ```

3. **Run with verbose output**:
   ```bash
   nextflow run main.nf -profile test,singularity -with-trace -with-report
   ```

---

## 10. Support and Contact Information

### 10.1 Documentation Resources

- **Pipeline Documentation**: This document
- **Nextflow Documentation**: https://www.nextflow.io/docs/latest/
- **Singularity Documentation**: https://sylabs.io/guides/latest/user-guide/

### 10.2 Issue Reporting

When reporting issues, please include:

1. **System information**:
   ```bash
   nextflow info
   uname -a
   singularity --version
   ```

2. **Configuration used**:
   ```bash
   cat nextflow.config
   ```

3. **Error messages**:
   ```bash
   cat .nextflow.log | tail -100
   ```

4. **Work directory logs** (for failed processes):
   ```bash
   cat work/xx/yyyyyy/.command.log
   cat work/xx/yyyyyy/.command.err
   ```

### 10.3 Community Support

- **GitHub Issues**: https://github.com/your-org/hla-typing-pipeline/issues
- **Discussions**: https://github.com/your-org/hla-typing-pipeline/discussions

### 10.4 Tool-Specific Support

For issues with individual HLA typing tools:

| Tool | Support Channel |
|------|-----------------|
| SpecHLA | https://github.com/deepomicslab/SpecHLA/issues |
| HLA-HD | https://www.genome.med.kyoto-u.ac.jp/HLA-HD/ |
| HLA*LA | https://github.com/DiltheyLab/HLA-LA/issues |
| arcasHLA | https://github.com/RabadanLab/arcasHLA/issues |
| OptiType | https://github.com/FRED-2/OptiType/issues |
| xHLA | https://github.com/humanlongevity/HLA/issues |

### 10.5 Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.2.0 | 2024-XX | Added LOH analysis, improved consensus voting |
| 1.1.0 | 2024-XX | Added FASTQ input support, MultiQC integration |
| 1.0.0 | 2024-XX | Initial release |

---

## Appendix A: Quick Reference Commands

### Basic Usage

```bash
# Single BAM, default tools (SpecHLA + HLA-HD)
nextflow run main.nf --input_bam sample.bam -profile singularity

# Single FASTQ pair
nextflow run main.nf --input_fastq_1 R1.fq.gz --input_fastq_2 R2.fq.gz -profile singularity

# Multiple samples
nextflow run main.nf --input_samplesheet samples.csv -profile singularity

# All tools
nextflow run main.nf --input_bam sample.bam --tools spechla,hlahd,hlala,arcashla,optitype,xhla -profile singularity

# SLURM cluster
nextflow run main.nf --input_samplesheet samples.csv -profile slurm,singularity

# Resume failed run
nextflow run main.nf -resume --input_samplesheet samples.csv -profile singularity
```

### Useful Nextflow Commands

```bash
# Show pipeline help
nextflow run main.nf --help

# View execution history
nextflow log

# Clean work directories
nextflow clean -f

# Clean runs older than 7 days
nextflow clean -f -before 7d

# Generate execution report
nextflow run main.nf -with-report report.html -with-timeline timeline.html -with-trace trace.txt
```

---

## Appendix B: Configuration Templates

### Template: SLURM HPC Configuration

```groovy
// slurm_hpc.config
params {
    container_dir = '/project/shared/hla/containers'
    hlahd_db = '/project/shared/hla/databases/hlahd'
    hlala_graphs = '/project/shared/hla/databases/hlala'
    max_cpus = 16
    max_memory = '128.GB'
    max_time = '72.h'
}

process {
    executor = 'slurm'
    queue = 'normal'
    clusterOptions = '--account=project123'

    withLabel: 'process_high' {
        queue = 'bigmem'
        clusterOptions = '--account=project123 --constraint=bigmem'
    }
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/project/shared/singularity_cache'
    runOptions = '--writable-tmpfs --cleanenv'
}

executor {
    queueSize = 100
    submitRateLimit = '10/1min'
}
```

### Template: PBS HPC Configuration

```groovy
// pbs_hpc.config
params {
    container_dir = '/work/containers/hla'
    max_cpus = 12
    max_memory = '64.GB'
}

process {
    executor = 'pbs'
    queue = 'workq'
    clusterOptions = '-A project123 -l select=1:ncpus=${task.cpus}:mem=${task.memory.toGiga()}gb'
}

singularity {
    enabled = true
    autoMounts = true
}
```

---

*End of Supplementary Documentation*
