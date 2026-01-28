# Quick Start Guide

This guide will help you get started with the HLA Typing Pipeline in 5 minutes.

## 🎯 Prerequisites

Before you begin, ensure you have:
- Nextflow installed (≥23.04.0)
- Docker or Singularity installed
- Your sequencing data ready (BAM, CRAM, or FASTQ files)

## 📦 Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/yourusername/hla-typing-pipeline.git
cd hla-typing-pipeline
```

### Step 2: Test Installation

```bash
nextflow -version  # Should show ≥23.04.0
nextflow run main.nf --help  # Should display help message
```

## 🚀 Your First Run

### Example 1: RNA-seq Data (Simplest)

If you have RNA-seq FASTQ files:

```bash
nextflow run main.nf \
    --input /path/to/fastq/files/ \
    --input_type fastq \
    --tools optitype,arcashla \
    --optitype_seq_type rna \
    --outdir results/ \
    -profile docker
```

### Example 2: DNA Exome Data

If you have exome BAM files:

```bash
nextflow run main.nf \
    --input /path/to/bam/files/ \
    --input_type bam \
    --tools spechla \
    --spechla_exon_only 1 \
    --outdir results/ \
    -profile docker
```

⚠️ **Critical**: Always use `--spechla_exon_only 1` for exome data!

### Example 3: WGS with Consensus Calling

For whole genome sequencing with multiple tools:

```bash
nextflow run main.nf \
    --input /path/to/fastq/files/ \
    --input_type fastq \
    --tools optitype,arcashla,spechla \
    --enable_majority_voting \
    --outdir results/ \
    -profile docker -resume
```

## 🏔️ CSC Puhti Specific

If you're on CSC Puhti (Finnish HPC):

### Step 1: Load Modules

```bash
module load java/21
module load biopython-env/3.10.6
module load nextflow/25.10.0
```

### Step 2: Edit SLURM Script

```bash
nano submit_slurm.sh
```

Modify these key parameters:
```bash
INPUT_DIR="/scratch/project_2008084/your_samples"
INPUT_TYPE="fastq"                    # or "bam"
TOOLS="optitype,arcashla"
SLURM_ACCOUNT="project_2008084"       # YOUR project number!
```

For exome data, also set:
```bash
SPECHLA_EXON_ONLY=1
```

### Step 3: Submit Job

```bash
mkdir -p logs
sbatch submit_slurm.sh
```

### Step 4: Monitor Progress

```bash
# Check job status
squeue -u $USER

# Watch log file (replace JOBID with your job ID)
tail -f logs/hla_pipeline_JOBID.out
```

## 📁 Understanding Results

After the pipeline completes, you'll find:

```
results/
├── sample1/
│   ├── optitype/         # OptiType results
│   │   └── sample1_result.tsv
│   ├── arcashla/         # ArcasHLA results
│   │   └── sample1.genotype.json
│   └── spechla/          # SpecHLA results
│       └── sample1_spechla.tsv
├── majority_voting/      # Consensus results (if enabled)
│   └── all_samples.consensus.tsv
└── pipeline_info/        # Pipeline execution reports
    └── execution_report.html
```

## 🔍 Interpreting Results

### OptiType Output (TSV)
```
A1      A2      B1      B2      C1      C2
A*02:01 A*24:02 B*15:01 B*44:03 C*03:04 C*05:01
```

### ArcasHLA Output (JSON)
```json
{
  "A": ["A*02:01", "A*24:02"],
  "B": ["B*15:01", "B*44:03"],
  "C": ["C*03:04", "C*05:01"]
}
```

### SpecHLA Output (TSV)
```
Sample  HLA_A_1  HLA_A_2  HLA_B_1  HLA_B_2  HLA_C_1  HLA_C_2
sample1 A*02:01  A*24:02  B*15:01  B*44:03  C*03:04  C*05:01
```

## ❓ Common Questions

### Q: Which tools should I use?

- **RNA-seq**: OptiType + ArcasHLA
- **DNA exome**: SpecHLA (with `--spechla_exon_only 1`)
- **WGS**: All three + majority voting

### Q: How do I resume a failed run?

Add `-resume` to your command:
```bash
nextflow run main.nf <your_params> -resume
```

### Q: My SpecHLA results are empty!

For exome data, you MUST use:
```bash
--spechla_exon_only 1
```

### Q: How much time/resources do I need?

Typical requirements per sample:
- **RNA-seq (OptiType+ArcasHLA)**: ~5-10 min, 8 CPUs, 16 GB RAM
- **Exome (SpecHLA exon-only)**: ~10-20 min, 8 CPUs, 32 GB RAM
- **WGS (all tools)**: ~30-60 min, 16 CPUs, 64 GB RAM

### Q: Can I use my own containers?

Yes! Edit `nextflow.config` and modify the container paths in the profiles section.

## 🆘 Getting Help

If something goes wrong:

1. **Check the help message**:
   ```bash
   nextflow run main.nf --help
   ```

2. **Look at the logs**:
   ```bash
   cat .nextflow.log
   cat results/sample_name/tool_name/*.log
   ```

3. **Check common issues** in the main README.md

4. **Open an issue** on GitHub with:
   - Your command
   - Error message
   - `.nextflow.log` file

## 🎓 Next Steps

Once you're comfortable with basic usage:

1. Read the full [README.md](README.md) for all parameters
2. Try enabling majority voting for consensus calling
3. Optimize resource usage for your cluster
4. Set up custom container paths if needed

## 📚 Useful Commands

```bash
# Show pipeline version
nextflow run main.nf --version

# Validate pipeline
nextflow run main.nf --help

# Clean work directory
nextflow clean -f

# Resume previous run
nextflow run main.nf <params> -resume

# Generate reports
nextflow run main.nf <params> \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-dag dag.html
```

## 🎉 You're Ready!

You now know enough to start HLA typing! Remember:
- Use the right input type
- Set `--spechla_exon_only 1` for exome data
- Use `-resume` to save time
- Check the logs if something fails

Good luck with your HLA typing! 🧬
