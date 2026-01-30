# HLA Typing Pipeline - Example Scenarios

This directory contains detailed tutorials for different HLA typing scenarios.

## Available Scenarios

### BAM Input Scenarios

| Scenario | Data Type | Document |
|----------|-----------|----------|
| **Scenario 1** | RNA-seq BAM | [SCENARIO_1_BAM_RNASEQ.md](SCENARIO_1_BAM_RNASEQ.md) |
| **Scenario 2** | WES BAM | [SCENARIO_2_BAM_WES.md](SCENARIO_2_BAM_WES.md) |
| **Scenario 3** | WGS BAM | [SCENARIO_3_BAM_WGS.md](SCENARIO_3_BAM_WGS.md) |

### FASTQ Input Scenarios

| Scenario | Data Type | Document |
|----------|-----------|----------|
| **Scenario 4** | RNA-seq FASTQ | [SCENARIO_4_FASTQ_RNASEQ.md](SCENARIO_4_FASTQ_RNASEQ.md) |
| **Scenario 5** | WES/WGS FASTQ | [SCENARIO_5_FASTQ_DNA.md](SCENARIO_5_FASTQ_DNA.md) |

### Quick Reference

| Document | Description |
|----------|-------------|
| [COMPARISON_SUMMARY.md](COMPARISON_SUMMARY.md) | Side-by-side comparison of all scenarios |

---

## Which Scenario Should I Use?

```
                     ┌─────────────────────────┐
                     │   What's your input?    │
                     └───────────┬─────────────┘
                                 │
              ┌──────────────────┴──────────────────┐
              ▼                                      ▼
        ┌───────────┐                          ┌───────────┐
        │    BAM    │                          │   FASTQ   │
        └─────┬─────┘                          └─────┬─────┘
              │                                      │
    ┌─────────┼─────────┐                  ┌─────────┴─────────┐
    ▼         ▼         ▼                  ▼                   ▼
┌───────┐ ┌───────┐ ┌───────┐         ┌───────┐           ┌───────┐
│RNA-seq│ │  WES  │ │  WGS  │         │RNA-seq│           │WES/WGS│
└───┬───┘ └───┬───┘ └───┬───┘         └───┬───┘           └───┬───┘
    │         │         │                 │                   │
    ▼         ▼         ▼                 ▼                   ▼
Scenario 1 Scenario 2 Scenario 3    Scenario 4           Scenario 5
```

---

## Quick Start Commands

### RNA-seq (BAM)
```bash
nextflow run main.nf \
    --input_samplesheet samples.csv \
    --tools arcashla,optitype,spechla \
    --seq_type rna \
    -profile singularity
```

### WES (BAM)
```bash
nextflow run main.nf \
    --input_samplesheet samples.csv \
    --tools spechla,hlahd,arcashla,optitype \
    -profile singularity
```

### WGS (BAM)
```bash
nextflow run main.nf \
    --input_samplesheet samples.csv \
    --tools spechla,hlahd,hlala,xhla \
    -profile singularity
```

### FASTQ (any type)
```bash
nextflow run main.nf \
    --input_samplesheet samples_fastq.csv \
    --tools spechla,hlahd,arcashla,optitype \
    -profile singularity
```

---

## Each Scenario Includes

1. **Overview** - Data type, tools, expected runtime
2. **Sample preparation** - Upload and verification steps
3. **Samplesheet creation** - Template and auto-generation
4. **SLURM script** - Ready-to-use submission script
5. **Monitoring** - How to track job progress
6. **Results interpretation** - Understanding outputs
7. **Troubleshooting** - Common issues and solutions
8. **Clinical/research templates** - Report formats

---

## Need Help?

- **Full documentation**: [../../README.md](../../README.md)
- **Puhti guide**: [../PUHTI_GUIDE.md](../PUHTI_GUIDE.md)
- **Issues**: [GitHub Issues](https://github.com/uoozcan/local_pipeline/issues)
