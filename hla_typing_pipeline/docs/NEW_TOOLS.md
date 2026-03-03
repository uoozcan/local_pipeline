# New Tools: BAMQC and flow-OptiType

This document describes the new tools added to the HLA typing pipeline in the `feature/bamqc-flow-optitype` branch.

## BAMQC - Comprehensive BAM Quality Control

### Overview
BAMQC provides detailed quality control metrics for BAM files, offering insights into:
- Read statistics (total, mapped, paired, duplicates)
- Mapping quality metrics
- Coverage statistics
- Insert size distribution
- Base quality distribution

### Container
- **Docker Image**: `kennethlim206/bamqc:latest`
- **Singularity**: Automatically converted to `.sif` format by `install_puhti.sh`

### Usage
Enable BAMQC by adding the `--run_bamqc` flag:

```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --tools spechla,arcashla \
    --run_bamqc \
    -profile singularity
```

### Outputs
- `<sample>_bamqc_report.txt` - Summary QC report with key metrics
- `<sample>_bamqc/flagstat.txt` - Detailed flagstat output
- `<sample>_bamqc/stats.txt` - Comprehensive BAM statistics
- `<sample>_bamqc/idxstats.txt` - Index statistics per chromosome

### When to Use
- For comprehensive quality assessment beyond standard QC
- When detailed read statistics are required
- For troubleshooting mapping or sequencing issues

---

## flow-OptiType - Direct BAM Processing for HLA Typing

### Overview
flow-OptiType is NMDP Bioinformatics' implementation of OptiType that can process BAM files directly without intermediate FASTQ conversion. This approach is:
- **Faster**: No BAM-to-FASTQ conversion overhead
- **More efficient**: Directly extracts HLA-relevant reads
- **Streamlined**: Single-step HLA typing from aligned data

### Container
- **Docker Image**: `nmdpbioinformatics/flow-optitype:latest`
- **Singularity**: Automatically converted to `.sif` format by `install_puhti.sh`
- **Source**: [GitHub - nmdp-bioinformatics/flow-OptiType](https://github.com/nmdp-bioinformatics/flow-OptiType)

### Usage
Enable flow-OptiType by adding the `--use_flow_optitype` flag when using `optitype`:

```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --tools spechla,arcashla,optitype \
    --use_flow_optitype \
    -profile singularity
```

**Note**: The `optitype` tool name triggers flow-OptiType when `--use_flow_optitype` is set.

### Comparison: Standard OptiType vs flow-OptiType

| Feature | Standard OptiType | flow-OptiType |
|---------|------------------|---------------|
| Input | BAM → FASTQ → OptiType | BAM → OptiType (direct) |
| Speed | Slower (conversion overhead) | Faster |
| Disk Usage | Higher (intermediate FASTQ) | Lower |
| Memory | Standard | Standard |
| Accuracy | High | High (same algorithm) |
| HLA Genes | Class I (A, B, C) | Class I (A, B, C) |

### Outputs
- `<sample>_flow_optitype.txt` - Standard-format HLA typing results
- `<sample>/<sample>_result.tsv` - Original OptiType result file
- Compatible with consensus voting pipeline

### When to Use
- **Recommended**: For all BAM-based workflows
- When processing time is a concern
- When disk space is limited
- For large-scale HLA typing projects

### Technical Details
- Extracts reads from HLA region (chr6:28510120-33480577 for hg38)
- Runs OptiType's core algorithm directly on extracted reads
- Fallback mechanism: If direct BAM processing fails, automatically converts to FASTQ
- Fully compatible with the pipeline's consensus voting system

---

## Installation on CSC Puhti

Both containers are automatically pulled and configured by the installation script:

```bash
bash scripts/install_puhti.sh project_XXXXXXX
```

The script will:
1. Pull `kennethlim206/bamqc:latest` → `bamqc.sif`
2. Pull `nmdpbioinformatics/flow-optitype:latest` → `flow_optitype.sif`
3. Configure Nextflow to use these containers
4. Add process definitions to `conf/user.config`

---

## Example Workflows

### Complete Pipeline with New Tools
```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --tools spechla,arcashla,optitype \
    --run_bamqc \
    --use_flow_optitype \
    --outdir results \
    -profile singularity
```

### BAMQC Only (Quality Control)
```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --tools spechla \
    --run_bamqc \
    -profile singularity
```

### flow-OptiType for Fast Class I Typing
```bash
nextflow run main.nf \
    --input_bam sample.bam \
    --tools optitype \
    --use_flow_optitype \
    -profile singularity
```

---

## References

- **BAMQC**: Docker Hub - [kennethlim206/bamqc](https://hub.docker.com/r/kennethlim206/bamqc)
- **flow-OptiType**:
  - GitHub: [nmdp-bioinformatics/flow-OptiType](https://github.com/nmdp-bioinformatics/flow-OptiType)
  - OptiType Paper: Szolek et al., Bioinformatics, 30(23):3310-6 (2014)
  - Nextflow: Di Tommaso et al., Nat Biotech, 35(4):316-319 (2017)

---

## Troubleshooting

### BAMQC Issues
- **Container not found**: Ensure `bamqc.sif` exists in `containers/` directory
- **Permission errors**: Check that BAM file is readable and indexed

### flow-OptiType Issues
- **No reads extracted**: Verify BAM file uses correct reference (hg38 expected)
- **Fallback to FASTQ**: Automatic, but slower - check HLA region coordinates if persistent
- **Memory errors**: Increase process memory in `nextflow.config` (default: 8GB)

---

*Generated for feature branch: `feature/bamqc-flow-optitype`*
*Date: February 2026*
