# SpecHLA BAM Analysis Implementation Summary

## Overview

This package provides a complete solution for analyzing your BAM file (`alignment-sorted.bam`) with SpecHLA through your Nextflow HLA typing pipeline on CSC Puhti.

**Important Context**: Your data is RNA-seq, which presents specific considerations for SpecHLA analysis. See recommendations below.

---

## What Has Been Created

### 1. Core Scripts

#### `step3_prepare_bam_input.sh`
**Purpose**: Validate and prepare BAM file for SpecHLA analysis

**What it does**:
- ✓ Validates BAM file existence and format
- ✓ Creates/checks BAM index
- ✓ Analyzes reference genome compatibility (CRITICAL)
- ✓ Detects chromosome naming convention (chr6 vs 6)
- ✓ Counts reads in HLA region
- ✓ Assesses data type and coverage
- ✓ Creates sample sheet for Nextflow
- ✓ Generates configuration file
- ✓ Provides specific recommendations

**Run time**: ~5-15 minutes

#### `step4_run_spechla_bam.sh`
**Purpose**: Execute SpecHLA analysis through Nextflow pipeline

**What it does**:
- ✓ Loads required modules (Nextflow, Singularity)
- ✓ Configures pipeline parameters
- ✓ Sets up HLA typing for all specified genes
- ✓ Executes SpecHLA with optimal settings
- ✓ Generates HTML reports and visualizations
- ✓ Handles errors with retry logic

**Run time**: ~2-4 hours for single BAM

### 2. Validation Tools

#### `validate_bam_for_spechla.sh`
**Purpose**: Pre-flight compatibility checker

**Tests performed**:
1. File existence and accessibility
2. BAM index availability
3. Format validation (sorted, indexed)
4. Chromosome naming detection
5. HLA-aware reference check (MOST CRITICAL)
6. Coverage assessment in HLA region
7. Data type identification
8. Read quality metrics

**Output**: Detailed compatibility report with specific recommendations

### 3. Configuration Files

#### `spechla_bam.config`
**Purpose**: Nextflow configuration for BAM analysis

**Includes**:
- Input/output settings
- SpecHLA-specific parameters
- Resource allocation
- Error handling
- Container configuration
- Puhti-specific optimizations

### 4. Documentation

#### `SPECHLA_BAM_ANALYSIS_GUIDE.md`
**Comprehensive 50+ page guide covering**:
- Complete requirements analysis
- Reference genome setup instructions
- Step-by-step workflow
- Pipeline integration details
- Extensive troubleshooting section
- Known limitations
- Best practices
- Comparison with other tools

#### `SPECHLA_BAM_QUICK_REFERENCE.txt`
**Quick reference card with**:
- Essential commands
- Common issues and fixes
- Parameter settings
- Monitoring instructions
- Decision trees

---

## Critical Requirements Analysis

### 1. **HLA-Aware Reference Genome** ⚠️ MOST CRITICAL

**The Issue**: SpecHLA absolutely requires that your BAM was aligned to a reference genome containing:
- HLA gene sequences with alternative contigs (ALT contigs)
- Full HLA haplotype information
- Proper ALT contig annotations

**Check your BAM**:
```bash
samtools view -H alignment-sorted.bam | grep -i HLA
```

**If you see HLA-specific contigs**: ✅ Good to proceed with SpecHLA  
**If you don't see any**: ⚠️ SpecHLA will have LIMITED accuracy

**References that work**:
- ✅ bwakit hs38DH.fa (recommended)
- ✅ GRCh38 full analysis set with HLA ALT contigs
- ❌ Standard GRCh38 without ALT contigs
- ❌ GRCh37/hg19

### 2. **Data Type Compatibility**

**Your data**: RNA-seq (bulk)

| Tool | RNA-seq Suitability | Recommendation |
|------|-------------------|----------------|
| SpecHLA | Limited | Use only with HLA-aware reference |
| OptiType | Excellent | **PRIMARY CHOICE** for RNA-seq |
| ArcasHLA | Excellent | **SECONDARY CHOICE** for RNA-seq |

**Why this matters**: 
- SpecHLA is optimized for WGS, not RNA-seq
- RNA-seq has expression variability that affects SpecHLA's algorithm
- OptiType and ArcasHLA are specifically designed for transcriptome data

### 3. **Coverage Requirements**

Minimum recommended:
- **WGS**: 30x average, 10,000+ reads in HLA region
- **RNA-seq**: Variable (depends on expression)
- **Exome**: Usually insufficient for SpecHLA

Your BAM will be checked during step 3.

---

## Decision Framework

### Should You Use SpecHLA for Your RNA-seq BAM?

```
START
  │
  ├─ Was BAM aligned to HLA-aware reference (with ALT contigs)?
  │  │
  │  ├─ YES → SpecHLA can work, but consider:
  │  │         • OptiType + ArcasHLA are better for RNA-seq
  │  │         • Run all three and use majority voting
  │  │         • SpecHLA as validation tool
  │  │
  │  └─ NO → SpecHLA NOT RECOMMENDED
  │           • Results will be inaccurate
  │           • Use OptiType + ArcasHLA instead
  │           • OR realign to HLA-aware reference
  │
  └─ Recommendation: Use OptiType (RNA mode) + ArcasHLA
                      Add SpecHLA only if you have proper reference
```

---

## Implementation Plan

### Phase 1: Pre-Flight Validation (15 minutes)

```bash
# SSH to Puhti
ssh ozcanumu@puhti.csc.fi

# Navigate to analysis directory
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis

# Upload scripts to scripts/ directory
cd scripts

# Load required modules
module load biokit

# Run validation (IMPORTANT - do this first!)
bash validate_bam_for_spechla.sh ../raw_fastq/alignment-sorted.bam
```

**Review the validation report carefully!** It will tell you:
- Whether SpecHLA is suitable for your data
- What issues need to be resolved
- Specific recommendations for your situation

### Phase 2: Reference Genome Setup (IF NEEDED)

**If validation shows "No HLA ALT contigs"**:

```bash
# Create references directory
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
mkdir -p references
cd references

# Download bwakit (contains HLA-aware reference)
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2
tar -xjf bwakit-0.7.15_x64-linux.tar.bz2

# Copy reference to accessible location
cp bwakit/resource-GRCh38/hs38DH.fa ./
cp bwakit/resource-GRCh38/hs38DH.fa.alt ./

# Index the reference
module load biokit
samtools faidx hs38DH.fa
```

**If your BAM was NOT aligned to this reference**, you would need to realign your data (time-consuming) OR use OptiType/ArcasHLA instead.

### Phase 3: Prepare Input (5-15 minutes)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Make executable
chmod +x step3_prepare_bam_input.sh

# Submit as SLURM job
sbatch step3_prepare_bam_input.sh

# Monitor progress
squeue -u ozcanumu
tail -f ../logs/prepare_bam_*.log
```

**This will create**:
- `pipeline_input_bam/` directory
- `alignment-sorted.bam` symlink
- `samplesheet.csv` for Nextflow
- `spechla_config.txt` with detected settings

### Phase 4: Configure Pipeline (5 minutes)

```bash
# Edit the run script
nano step4_run_spechla_bam.sh

# Find and update these lines (around line 97-99):
REF_GENOME="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/references/hs38DH.fa"

# This MUST match the reference used to create your BAM!
# If you don't have hs38DH.fa, specify your actual reference path
```

**Other parameters to verify**:
```bash
SPECHLA_GENES="A,B,C,DRB1,DQB1,DPB1"  # Which genes to type
SPECHLA_RESOLUTION="2field"            # Resolution level
CHR6_NAME="chr6"                       # Detected in step3, verify
```

### Phase 5: Run Analysis (2-4 hours)

```bash
# Make executable
chmod +x step4_run_spechla_bam.sh

# Submit job
sbatch step4_run_spechla_bam.sh

# Check status
squeue -u ozcanumu

# Monitor log
tail -f ../logs/spechla_bam_*.log

# Check Nextflow progress
cd ../pipeline
tail -f .nextflow.log
```

### Phase 6: Review Results

```bash
# View results
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
cd results_bam_spechla/spechla/alignment-sorted/

# Main result file
cat alignment-sorted.spechla.result.txt

# Example output format:
# Gene    Allele1         Allele2         Confidence
# A       A*02:01         A*03:01         0.95
# B       B*07:02         B*44:03         0.92
# ...

# Copy HTML reports to local machine
scp ozcanumu@puhti.csc.fi:/scratch/project_2008084/ozcanumu/\
hla_rnaseq_analysis/results_bam_spechla/*.html .
```

---

## Alternative Recommendation (Preferred for RNA-seq)

Given that your data is RNA-seq, here's a better approach:

### Use Your Existing FASTQ-Based Pipeline

You already have scripts (`step1_setup_puhti.sh` through `step4_run_pipeline.sh`) that work with FASTQ files. This approach is **better suited for RNA-seq**:

```bash
# This uses the workflow you've already set up
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Run your existing pipeline with optimal tools for RNA-seq
sbatch step4_run_pipeline.sh

# Which uses:
# --tools "optitype,arcashla"  (or add spechla if you want)
# --optitype_seq_type rna      (RNA mode!)
# --enable_majority_voting true
```

**Why this is better**:
1. OptiType RNA mode is designed for RNA-seq
2. ArcasHLA is optimized for transcriptome data
3. No HLA-aware reference requirement
4. Majority voting increases confidence
5. You've already set this up!

---

## Troubleshooting Common Issues

### Issue 1: "No HLA ALT contigs found"

**Impact**: SpecHLA will work but with limited accuracy

**Solutions**:
1. **Best**: Realign your RNA-seq data to hs38DH.fa
2. **Recommended**: Use OptiType + ArcasHLA instead
3. **Acceptable**: Continue with SpecHLA but validate with other tools

### Issue 2: "No reads in HLA region"

**Check**:
```bash
# Try both chromosome naming conventions
samtools view -c alignment-sorted.bam chr6:28000000-34000000
samtools view -c alignment-sorted.bam 6:28000000-34000000
```

**Fix**: Update `CHR6_NAME` parameter in step4 script

### Issue 3: "Reference genome not found"

**Check**:
```bash
ls -lh /path/to/reference.fa
```

**Fix**: Update `REF_GENOME` path in step4 script

### Issue 4: "Low coverage warning"

**For RNA-seq**: This is expected due to variable expression

**Options**:
- Accept the warning and proceed
- Use OptiType/ArcasHLA which handle this better
- Check if HLA genes are expressed in your samples

---

## File Locations Summary

```
/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/
│
├── raw_fastq/
│   └── alignment-sorted.bam                    # Your BAM file
│
├── scripts/
│   ├── step3_prepare_bam_input.sh              # Preparation script
│   ├── step4_run_spechla_bam.sh                # Execution script
│   ├── validate_bam_for_spechla.sh             # Validation script
│   ├── spechla_bam.config                      # Configuration
│   ├── SPECHLA_BAM_ANALYSIS_GUIDE.md           # Full documentation
│   └── SPECHLA_BAM_QUICK_REFERENCE.txt         # Quick reference
│
├── pipeline_input_bam/                          # Created by step3
│   ├── alignment-sorted.bam -> ../raw_fastq/...
│   ├── alignment-sorted.bam.bai
│   ├── samplesheet.csv
│   └── spechla_config.txt
│
├── results_bam_spechla/                         # Created by step4
│   ├── spechla/
│   │   └── alignment-sorted/
│   │       ├── alignment-sorted.spechla.result.txt
│   │       └── alignment-sorted.spechla.detailed.txt
│   ├── spechla_report.html
│   ├── spechla_timeline.html
│   └── spechla_dag.svg
│
├── logs/                                        # Job logs
│   ├── prepare_bam_*.log
│   └── spechla_bam_*.log
│
└── references/                                  # HLA-aware references
    └── hs38DH.fa
```

---

## Next Steps

### Immediate Actions (Required)

1. **Upload scripts** to Puhti:
   ```bash
   scp step3_prepare_bam_input.sh step4_run_spechla_bam.sh \
       validate_bam_for_spechla.sh spechla_bam.config \
       ozcanumu@puhti.csc.fi:/scratch/project_2008084/ozcanumu/\
       hla_rnaseq_analysis/scripts/
   ```

2. **Run validation** (CRITICAL):
   ```bash
   bash validate_bam_for_spechla.sh ../raw_fastq/alignment-sorted.bam
   ```

3. **Review validation report** and make decision based on findings

### Decision Points

**If validation shows HLA-aware reference**:
- ✓ Proceed with SpecHLA
- ✓ Consider also running OptiType/ArcasHLA for comparison
- ✓ Use majority voting for best results

**If validation shows NO HLA-aware reference**:
- Option A: Use OptiType + ArcasHLA (recommended for RNA-seq)
- Option B: Realign data to HLA-aware reference
- Option C: Proceed with SpecHLA but understand limitations

**If validation shows low coverage**:
- For RNA-seq: Expected, proceed but validate results
- For WGS: May need deeper sequencing

---

## Expected Outcomes

### Success Criteria

After completion, you should have:
- ✓ HLA typing results for genes A, B, C, DRB1, DQB1, DPB1
- ✓ Confidence scores for each allele call
- ✓ Detailed coverage information
- ✓ HTML reports with visualizations
- ✓ Timeline and DAG of pipeline execution

### Result Interpretation

**High confidence (>0.9)**: Strong evidence for allele call  
**Medium confidence (0.7-0.9)**: Good evidence, consider validation  
**Low confidence (<0.7)**: Ambiguous, requires validation with other tools

For RNA-seq data, lower confidence is expected due to:
- Variable gene expression levels
- Some alleles may not be expressed
- Splicing complexity

---

## Support and Resources

### Documentation Files
1. `SPECHLA_BAM_ANALYSIS_GUIDE.md` - Complete guide (50+ pages)
2. `SPECHLA_BAM_QUICK_REFERENCE.txt` - Quick reference card
3. This file - Implementation summary

### CSC Support
- Email: servicedesk@csc.fi
- Documentation: https://docs.csc.fi
- Status: https://status.csc.fi

### Tool Documentation
- SpecHLA: https://github.com/deepomicslab/SpecHLA
- OptiType: https://github.com/FRED-2/OptiType
- ArcasHLA: https://github.com/RabadanLab/arcasHLA

---

## Final Recommendations

### For Your RNA-seq Data

**Primary Recommendation**:
```bash
# Use your existing FASTQ pipeline with OptiType + ArcasHLA
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
sbatch step4_run_pipeline.sh  # Already configured for RNA-seq!
```

**Secondary Option** (if you want to try SpecHLA):
```bash
# First validate compatibility
bash validate_bam_for_spechla.sh ../raw_fastq/alignment-sorted.bam

# If compatible, run SpecHLA
sbatch step3_prepare_bam_input.sh
# Configure reference in step4_run_spechla_bam.sh
sbatch step4_run_spechla_bam.sh

# Compare with OptiType/ArcasHLA results
```

**Best Practice**:
- Run all three tools (OptiType, ArcasHLA, SpecHLA)
- Use majority voting to reconcile differences
- Focus on calls agreed upon by multiple tools

---

## Checklist

Before running SpecHLA:

- [ ] Scripts uploaded to Puhti
- [ ] Validation script executed
- [ ] Validation report reviewed
- [ ] Decision made based on validation
- [ ] Reference genome checked/configured
- [ ] BAM file is coordinate-sorted
- [ ] BAM index exists
- [ ] Chromosome naming verified
- [ ] step4 script configured with correct reference
- [ ] Sufficient disk space available
- [ ] Alternative tools considered (OptiType/ArcasHLA)

---

**Document**: Implementation Summary  
**Version**: 1.0  
**Date**: December 2024  
**Author**: Umut Ozcan  
**Project**: project_2008084 (CSC Puhti)
