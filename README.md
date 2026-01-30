# HLA Typing Pipeline

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A522.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/Singularity-%E2%89%A53.0-blue.svg)](https://sylabs.io/singularity/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive **multi-tool HLA typing pipeline** with weighted consensus voting for reliable allele determination from WGS, WES, and RNA-seq data.

---

## Why This Pipeline?

Human Leukocyte Antigen (HLA) typing is critical for:

| Application | Description |
|-------------|-------------|
| **Transplantation** | Donor-recipient compatibility matching |
| **Pharmacogenomics** | Drug hypersensitivity prediction (e.g., HLA-B*57:01 and abacavir) |
| **Disease Studies** | Autoimmune disease risk assessment |
| **Cancer Immunotherapy** | Neoantigen prediction and vaccine design |

**The Challenge:** Different HLA typing tools often produce inconsistent results.

**Our Solution:** Integrate **6 state-of-the-art tools** and apply **weighted consensus voting** to produce reliable, high-confidence HLA calls.

---

## Pipeline Overview

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

## Key Features

| Feature | Description |
|---------|-------------|
| **6 HLA Tools** | SpecHLA, HLA-HD, HLA*LA, arcasHLA, OptiType, xHLA |
| **Consensus Voting** | Weighted voting based on read confidence |
| **Dual Input** | BAM files or paired-end FASTQ files |
| **Quality Control** | Automatic QC with warnings for low-quality samples |
| **Visual Reports** | Interactive HTML reports with charts and heatmaps |
| **LOH Detection** | Loss of Heterozygosity analysis for tumor samples |
| **HPC Ready** | SLURM cluster and CSC Puhti supercomputer support |
| **Containerized** | Singularity/Docker for full reproducibility |

---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/uoozcan/local_pipeline.git
cd local_pipeline/hla_typing_pipeline
```

### 2. Run with BAM File

```bash
nextflow run main.nf \
    --input_bam /path/to/sample.bam \
    --outdir results \
    --tools spechla,hlahd \
    -profile singularity
```

### 3. Run with FASTQ Files

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
cat results/SampleName/SampleName_consensus.txt
```

**Example Output:**
```
Gene      Allele1      Allele2      Confidence
HLA-A     A*02:01      A*03:01      0.95
HLA-B     B*07:02      B*44:02      0.88
HLA-C     C*07:01      C*05:01      1.00
HLA-DRB1  DRB1*15:01   DRB1*04:01   0.92
```

---

## HLA Tools Comparison

| Tool | Algorithm | Class I | Class II | Best For |
|------|-----------|:-------:|:--------:|----------|
| **SpecHLA** | Assembly + Alignment | ✓ | ✓ | High resolution, phasing |
| **HLA-HD** | Exhaustive alignment | ✓ | ✓ | Accuracy, rare alleles |
| **HLA*LA** | Graph-based | ✓ | ✓ | WGS data |
| **arcasHLA** | Kallisto pseudoalign | ✓ | ✓ | RNA-seq, speed |
| **OptiType** | Integer programming | ✓ | ✗ | Class I accuracy |
| **xHLA** | K-mer based | ✓ | ✓ | Speed, WGS |

### Recommended Tool Combinations

| Data Type | Recommended Tools |
|-----------|-------------------|
| **WGS** | `spechla,hlahd,hlala` |
| **WES** | `spechla,hlahd,arcashla` |
| **RNA-seq** | `arcashla,optitype,spechla` |
| **Quick analysis** | `arcashla,optitype` |

---

## Installation

### Prerequisites

- **Nextflow** ≥22.10.0
- **Singularity** ≥3.0 (or Docker)
- **Samtools** ≥1.10

### Automated Setup (CSC Puhti)

```bash
# Clone repository
git clone https://github.com/uoozcan/local_pipeline.git
cd local_pipeline

# Run setup script
./hla_typing_pipeline/scripts/puhti_full_setup.sh
```

### Manual Setup

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Pull containers
mkdir -p hla_references/containers
singularity pull hla_references/containers/hlahd.sif docker://quay.io/biocontainers/hlahd:1.7.0
singularity pull hla_references/containers/arcashla.sif docker://quay.io/biocontainers/arcas-hla:0.5.0
```

---

## Output Structure

```
results/
├── SampleName/
│   ├── SampleName_consensus.txt      # Final HLA calls
│   ├── SampleName_comparison.txt     # Tool comparison
│   ├── qc/                           # Quality metrics
│   ├── spechla/                      # SpecHLA results
│   ├── hlahd/                        # HLA-HD results
│   ├── visualizations/
│   │   ├── SampleName_confidence.png
│   │   ├── SampleName_agreement.png
│   │   └── SampleName_report.html    # Interactive report
│   └── loh/                          # LOH analysis (if enabled)
├── summary/
│   └── hla_summary_report.html       # Multi-sample summary
└── multiqc/
    └── multiqc_report.html           # Aggregated QC
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [Full Documentation](hla_typing_pipeline/README.md) | Complete pipeline documentation |
| [Puhti Guide](hla_typing_pipeline/docs/PUHTI_GUIDE.md) | CSC Puhti supercomputer guide |
| [Parameters](hla_typing_pipeline/README.md#parameters) | All available parameters |
| [Troubleshooting](hla_typing_pipeline/README.md#troubleshooting) | Common issues and solutions |

---

## Repository Structure

```
local_pipeline/
├── hla_typing_pipeline/      # Main Nextflow pipeline
│   ├── main.nf               # Pipeline entry point
│   ├── nextflow.config       # Configuration
│   ├── modules/              # Tool modules (*.nf)
│   ├── bin/                  # Python/shell scripts
│   ├── conf/                 # Profile configs
│   ├── scripts/              # Setup & submission scripts
│   └── docs/                 # Documentation
├── scripts/                  # Standalone tool runners
│   ├── run_spechla.sh
│   ├── run_hlahd.sh
│   └── ...
└── spechla_local/            # Local SpecHLA installation
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
  <b>Reliable HLA typing for clinical and research applications</b>
</p>
