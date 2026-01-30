# HLA Typing Pipeline

Multi-tool HLA typing pipeline with consensus voting for reliable allele determination from WGS/WES/RNA-seq data.

## Overview

This pipeline integrates multiple HLA typing tools and applies majority voting to determine consensus HLA alleles with confidence scores. It supports:

- **SpecHLA**: High-resolution HLA typing with phasing capability
- **HLA-HD**: Accurate HLA typing using an exhaustive alignment approach
- **HLA\*LA**: Graph-based HLA typing using population reference graphs (BAM only)
- **arcasHLA**: Fast HLA typing from RNA-seq data
- **OptiType**: HLA Class I typing from various data types

## Features

- **Dual input support**: BAM files or paired FASTQ files
- Multi-tool HLA typing with configurable tool selection
- **Quality control with warnings** for insufficient HLA reads
- **Weighted consensus voting** based on read confidence
- Support for 2-field and 4-field resolution
- Configurable minimum tool agreement threshold
- **Comprehensive visualizations**:
  - Per-sample confidence charts and tool agreement heatmaps
  - Multi-sample summary reports with allele frequency plots
  - Interactive HTML reports
- **FastQC and MultiQC integration** for input data quality
- Singularity/Docker container support
- SLURM cluster compatibility

## Requirements

- Nextflow >= 22.10.0
- Singularity >= 3.0 or Docker
- Reference databases (see Setup)

## Quick Start

```bash
# Clone repository
git clone https://github.com/your-org/hla-typing-pipeline.git
cd hla-typing-pipeline

# Run with single BAM file
nextflow run main.nf \
    --input_bam sample.bam \
    --outdir results \
    -profile singularity

# Run with paired FASTQ files
nextflow run main.nf \
    --input_fastq_1 sample_R1.fastq.gz \
    --input_fastq_2 sample_R2.fastq.gz \
    --outdir results \
    -profile singularity

# Run with multiple tools
nextflow run main.nf \
    --input_bam sample.bam \
    --tools hlahd,arcashla \
    --outdir results \
    -profile singularity
```

## Input

### Single BAM File

```bash
nextflow run main.nf --input_bam /path/to/sample.bam
```

### Single Paired FASTQ Files

```bash
nextflow run main.nf \
    --input_fastq_1 /path/to/sample_R1.fastq.gz \
    --input_fastq_2 /path/to/sample_R2.fastq.gz
```

### Multiple Samples (Samplesheet)

#### BAM Samplesheet
Create a CSV file `samples_bam.csv`:
```csv
sample_id,bam_path
Sample1,/path/to/Sample1.bam
Sample2,/path/to/Sample2.bam
```

#### FASTQ Samplesheet
Create a CSV file `samples_fastq.csv`:
```csv
sample_id,fastq_1,fastq_2
Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
Sample2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
```

Run:
```bash
nextflow run main.nf --input_samplesheet samples_bam.csv
# or
nextflow run main.nf --input_samplesheet samples_fastq.csv
```

The pipeline automatically detects the input type from the samplesheet header.

## Parameters

### Input/Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_bam` | - | Path to input BAM file |
| `--input_fastq_1` | - | Path to R1 FASTQ file |
| `--input_fastq_2` | - | Path to R2 FASTQ file |
| `--input_samplesheet` | - | Path to samplesheet CSV |
| `--outdir` | `./results` | Output directory |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reference` | `hg38` | Reference genome (hg38 or hg19) |
| `--tools` | `spechla,hlahd` | HLA typing tools (comma-separated) |
| `--seq_type` | `dna` | Sequence type for OptiType (dna or rna) |
| `--hlala_graph` | `PRG_MHC_GRCh38_withIMGT` | HLA*LA graph to use |
| `--hla_genes` | classical | HLA genes to type |
| `--resolution` | `2-field` | Output resolution (2-field or 4-field) |
| `--min_tools` | `1` | Minimum tools for consensus |

### QC Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_hla_reads` | `1000` | Minimum HLA reads (warning threshold) |
| `--min_read_length` | `50` | Minimum average read length |

### Consensus Weighting Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--expected_reads` | `1000` | Expected reads per allele (for confidence calculation) |
| `--weighting` | `read_confidence` | Weighting method: `equal`, `read_confidence`, or `tool_quality` |

### LOH (Loss of Heterozygosity) Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_loh` | `false` | Enable LOH analysis (requires SpecHLA) |
| `--tumor_purity` | - | Tumor purity estimate (0-1), required for LOH |
| `--tumor_ploidy` | - | Tumor ploidy estimate, required for LOH |
| `--loh_het_cutoff` | `5` | Minimum heterozygous SNPs for LOH call |

### Resource Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `8` | Maximum CPUs per process |
| `--max_memory` | `32.GB` | Maximum memory per process |

### Tool Installation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--spechla_path` | `${projectDir}/../spechla_local` | Path to local SpecHLA installation |
| `--use_local_spechla` | `true` | Use local SpecHLA instead of container |

## Output

```
results/
├── SampleName/
│   ├── qc/
│   │   └── SampleName_qc_report.txt     # QC metrics and warnings
│   ├── fastqc/
│   │   ├── SampleName_fastqc.html       # FastQC HTML report
│   │   └── SampleName_fastqc.zip        # FastQC data archive
│   ├── spechla/
│   │   ├── SampleName_spechla.txt       # Parsed SpecHLA results
│   │   └── SampleName/                   # Full SpecHLA output
│   ├── hlahd/
│   │   ├── SampleName_hlahd.txt         # Parsed HLA-HD results
│   │   └── result/                       # Full HLA-HD output
│   ├── hlala/
│   │   ├── SampleName_hlala.txt         # Parsed HLA*LA results
│   │   └── SampleName/                   # Full HLA*LA output
│   ├── arcashla/
│   │   ├── SampleName_arcashla.txt      # Parsed arcasHLA results
│   │   └── SampleName/                   # Full arcasHLA output
│   ├── visualizations/
│   │   ├── SampleName_confidence.png    # Confidence bar chart
│   │   ├── SampleName_coverage.png      # Read coverage plot
│   │   ├── SampleName_agreement.png     # Tool agreement heatmap
│   │   ├── SampleName_report.html       # Interactive HTML report
│   │   └── SampleName_statistics.json   # Statistics in JSON format
│   ├── loh/                              # LOH analysis (if enabled)
│   │   ├── SampleName_hla_loh.txt       # LOH results per gene
│   │   ├── SampleName_loh_plot.png      # Copy number visualization
│   │   └── SampleName_loh_report.html   # Interactive LOH report
│   ├── SampleName_consensus.txt         # Consensus HLA types
│   └── SampleName_comparison.txt        # Tool comparison matrix
├── summary/
│   ├── hla_summary_report.html          # Multi-sample summary report
│   ├── hla_summary_statistics.tsv       # Summary statistics table
│   ├── allele_frequency.png             # Allele frequency plot
│   └── sample_quality.png               # Sample quality comparison
├── multiqc/
│   ├── multiqc_report.html              # Aggregated QC report
│   └── multiqc_data/                    # MultiQC data files
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

### QC Report Format

The QC report provides statistics and warnings:

```
# HLA Typing QC Report for SampleName
# Generated: 2024-01-01
# Reference: hg38

=== Input Statistics ===
Total reads: 50000000
HLA region reads: 15000
HLA read percentage: 0.03%
Average read length: 150 bp
Average mapping quality: 45.2

=== QC Warnings ===
No warnings - input data appears suitable for HLA typing

=== QC Status ===
Status: true
```

#### Warning Levels

- **WARNING**: Results may be less reliable but analysis continues
  - Low HLA read count (< 1000 reads)
  - Short average read length (< 50 bp)
  - Low mapping quality (< 20)

- **CRITICAL**: Results are unlikely to be reliable
  - Very few HLA reads (< 100 reads)

### Consensus Output Format

```
# HLA Consensus Results for SampleName
# Tools used: hlahd, arcashla
# Resolution: 2-field
# Weighting method: read_confidence
# Expected reads per allele: 1000
# Minimum tools for consensus: 1
#
Gene    Allele1    Allele2    Confidence    Reads1    Reads2
HLA-A   A*03:01    A*02:01    0.95          1250      980
HLA-B   B*35:01    B*40:01    0.88          850       720
HLA-C   C*03:04    C*04:01    1.00          1100      1050
...
```

### Visualization Outputs

Each sample generates interactive visualizations:

- **Confidence Bar Chart**: Shows confidence scores per HLA gene
- **Read Coverage Plot**: Displays read support for each allele
- **Tool Agreement Heatmap**: Shows agreement between different HLA typing tools
- **Interactive HTML Report**: Comprehensive report with all statistics and plots

### Summary Report

For multi-sample runs, a summary report aggregates results across all samples:

- **Allele Frequency Plot**: Distribution of HLA alleles across samples
- **Sample Quality Comparison**: Compare QC metrics across samples
- **Confidence Distribution**: Overall confidence score distribution
- **Summary Statistics Table**: TSV file for downstream analysis

## Tool Support by Input Type

| Tool | BAM | FASTQ | Notes |
|------|-----|-------|-------|
| SpecHLA | ✓ | ✓ | Full support |
| HLA-HD | ✓ | ✓ | Full support |
| HLA*LA | ✓ | ✗ | BAM only (requires alignment) |
| arcasHLA | ✓ | ✓ | Full support |
| OptiType | ✓ | ✓ | Full support |

## Setup

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```

### 2. Prepare Containers

Download or build Singularity containers:

```bash
# SpecHLA container
singularity pull spechla_1.0.7-3.sif docker://your-registry/spechla:1.0.7

# HLA-HD container
singularity pull hlahd.sif docker://your-registry/hlahd:latest

# arcasHLA container
singularity pull arcashla.sif docker://your-registry/arcashla:latest

# OptiType container
singularity pull optitype.sif docker://your-registry/optitype:latest
```

Place containers in `hla_references/containers/`.

### 2b. Local SpecHLA Installation (Alternative)

Instead of using a container, you can use a local SpecHLA installation:

```bash
# Clone SpecHLA
git clone https://github.com/deepomicslab/SpecHLA.git spechla_local
cd spechla_local

# Follow SpecHLA installation instructions
# Make sure all dependencies are installed (bwa, samtools, etc.)
```

Configure in `nextflow.config`:
```groovy
params {
    spechla_path = "/path/to/spechla_local"
    use_local_spechla = true  // Set to false to use container instead
}
```

Or pass on command line:
```bash
nextflow run main.nf \
    --spechla_path /path/to/spechla_local \
    --use_local_spechla true \
    --input_bam sample.bam
```

### 3. Prepare Databases

#### HLA-HD Database
Download from the HLA-HD website and place in `hla_references/databases/hlahd_db/`.

#### arcasHLA Database
arcasHLA downloads its database automatically on first run.

## Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Use Singularity containers |
| `docker` | Use Docker containers |
| `slurm` | Run on SLURM cluster |
| `test` | Test with minimal resources |

### SLURM Configuration

Edit `nextflow.config` to set your account:

```groovy
profiles {
    slurm {
        process.clusterOptions = '--account=your_account'
    }
}
```

## Consensus Voting Algorithm

The pipeline uses a weighted voting algorithm based on read confidence:

1. **Collect alleles** from each tool for each HLA gene, including read counts
2. **Normalize alleles** to specified resolution (2-field or 4-field)
3. **Calculate read confidence** for each allele: `confidence = min(1.0, reads / expected_reads)`
4. **Apply weighted voting** using the selected weighting method:
   - `equal`: Each tool has equal weight (1.0)
   - `read_confidence`: Weight based on read support confidence
   - `tool_quality`: Weight based on historical tool accuracy
5. **Select top 2 alleles** as consensus based on weighted votes
6. **Calculate final confidence** as weighted proportion of agreeing tools

### Weighting Methods

| Method | Description | Best for |
|--------|-------------|----------|
| `equal` | All tools weighted equally (1.0) | Quick analysis, equal trust |
| `read_confidence` | Weight by read support ratio | General use (default) |
| `tool_quality` | Weight by tool accuracy scores | Established benchmarks |

### Confidence Interpretation

The confidence score reflects both tool agreement and read support:

- `1.00`: High agreement with strong read support
- `0.75-0.99`: Good agreement, most tools support the call
- `0.50-0.74`: Moderate agreement, consider reviewing
- `<0.50`: Low confidence, manual review recommended

## HLA Loss of Heterozygosity (LOH) Analysis

The pipeline includes optional LOH detection for tumor samples. This feature identifies HLA allele loss events that may contribute to immune evasion.

### Requirements

- **SpecHLA** must be included in tools (uses frequency data from SpecHLA)
- **Tumor purity** estimate (from tools like ABSOLUTE, ASCAT, or Sequenza)
- **Tumor ploidy** estimate

### Running LOH Analysis

```bash
nextflow run main.nf \
    --input_bam tumor_sample.bam \
    --tools spechla,hlahd \
    --run_loh true \
    --tumor_purity 0.75 \
    --tumor_ploidy 2.1 \
    -profile singularity
```

### LOH Output

The LOH analysis produces:

| File | Description |
|------|-------------|
| `*_hla_loh.txt` | Copy number and LOH status for each HLA gene |
| `*_loh_plot.png` | Visualization of allele copy numbers |
| `*_loh_report.html` | Interactive report with clinical interpretation |

### LOH Output Format

```
Sample  HLA  Allele1    Allele2    CopyRatio  KeptHLA    LostHLA    LOH
Sample1 A    A*02:01    A*03:01    2:0        A*02:01    A*03:01    Y
Sample1 B    B*07:02    B*44:02    1:1        B*07:02    B*44:02    N
```

### Clinical Significance

HLA LOH is a mechanism of immune evasion in tumors:
- **LOH = Y**: One HLA allele is lost or significantly reduced
- May impact response to immune checkpoint inhibitors
- Affects eligibility for personalized cancer vaccines
- Consider for patient stratification in immunotherapy trials

## Standalone Scripts

### Consensus Voting Script

The consensus voting script can be used independently:

```bash
./bin/consensus_voting.py \
    --sample SampleName \
    --tools hlahd,arcashla \
    --files hlahd_result.txt,arcashla_result.txt \
    --resolution 2-field \
    --min-tools 1 \
    --expected-reads 1000 \
    --weighting read_confidence \
    --output-consensus consensus.txt \
    --output-comparison comparison.txt
```

### Visualization Script

Generate visualizations for a single sample:

```bash
./bin/hla_visualize.py \
    --sample SampleName \
    --consensus consensus.txt \
    --comparison comparison.txt \
    --qc-report qc_report.txt \
    --output-dir ./visualizations
```

### Summary Report Script

Generate multi-sample summary report:

```bash
./bin/hla_summary_report.py \
    --consensus-files sample1_consensus.txt sample2_consensus.txt \
    --comparison-files sample1_comparison.txt sample2_comparison.txt \
    --statistics-files sample1_statistics.json sample2_statistics.json \
    --output-dir ./summary
```

## Troubleshooting

### Low HLA Read Warnings

If you receive warnings about low HLA reads:
- For WGS data: Ensure sufficient sequencing depth (>30x recommended)
- For RNA-seq: HLA expression may be low in some tissues
- Consider using HLA-enriched sequencing for better results

### Memory Issues

If HLA-HD fails with memory errors, increase memory in `nextflow.config`:

```groovy
withName: 'HLAHD' {
    memory = '32.GB'
}
```

### Container Permissions

If containers fail with permission errors, try:

```bash
nextflow run main.nf -profile singularity \
    --singularity_runOptions '--writable-tmpfs'
```

### Missing Index Files

The pipeline automatically creates BAM index files if missing. However, for large BAM files, pre-indexing is recommended:

```bash
samtools index sample.bam
```

## Citation

If you use this pipeline, please cite the original tools:

- **SpecHLA**: Liu et al. (2023) SpecHLA: High-resolution HLA typing from sequencing data
- **HLA-HD**: Kawaguchi et al. (2017) HLA-HD: An accurate HLA typing algorithm...
- **HLA\*LA**: Dilthey et al. (2019) HLA*LA: Fast HLA type inference from whole-genome data
- **arcasHLA**: Orenbuch et al. (2020) arcasHLA: high-resolution HLA typing...
- **OptiType**: Szolek et al. (2014) OptiType: precision HLA typing...

## License

MIT License

## Contact

For issues and feature requests, please open an issue on GitHub.
