# SpecHLA BAM Analysis - Complete Guide

## Table of Contents
1. [Overview](#overview)
2. [Critical Requirements](#critical-requirements)
3. [Reference Genome Requirements](#reference-genome-requirements)
4. [Step-by-Step Workflow](#step-by-step-workflow)
5. [Pipeline Integration](#pipeline-integration)
6. [Troubleshooting](#troubleshooting)
7. [Known Limitations](#known-limitations)

---

## Overview

This guide explains how to analyze BAM files with SpecHLA for HLA typing. SpecHLA is a specialized tool that requires careful attention to reference genome compatibility and data preparation.

**Your Current Setup:**
- BAM file: `/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/raw_fastq/alignment-sorted.bam`
- Pipeline: HLA typing pipeline v2.0.0+
- Environment: CSC Puhti supercomputer
- Project: project_2008084

---

## Critical Requirements

### 1. **HLA-Aware Reference Genome (MOST CRITICAL)**

SpecHLA **REQUIRES** that your BAM file was aligned to a reference genome that includes:
- HLA gene sequences with alternative contigs (ALT contigs)
- Full HLA haplotype information
- Proper ALT contig annotations

**❌ Will NOT work properly:**
- Standard GRCh38 without HLA ALT contigs
- GRCh37/hg19 references
- GENCODE references without HLA extensions
- Standard exome capture references

**✅ Will work properly:**
- bwakit hs38DH.fa (recommended)
- GRCh38 full analysis set with HLA ALT contigs
- IMGT/HLA-aware custom references

**Why this matters:**
SpecHLA uses the reference genome to distinguish between different HLA alleles based on ALT contig alignments. Without proper ALT contigs, it cannot accurately determine HLA types.

### 2. **Data Type Compatibility**

| Data Type | SpecHLA Compatibility | Recommended Alternative |
|-----------|----------------------|------------------------|
| Whole Genome Sequencing (WGS) | ✅ **Best choice** | - |
| RNA-seq (bulk) | ⚠️ **Possible** but limited | OptiType, ArcasHLA |
| RNA-seq (single-cell) | ❌ Not suitable | OptiType |
| Exome sequencing | ❌ **Not recommended** | OptiType, ArcasHLA |
| Targeted HLA sequencing | ✅ Yes (with proper ref) | - |

**Your BAM file appears to be RNA-seq data**, based on the directory structure. Consider:
- OptiType (RNA mode) as primary tool
- ArcasHLA for additional validation
- SpecHLA only if you have proper HLA-aware reference

### 3. **BAM File Requirements**

Your BAM must be:
- ✅ Coordinate-sorted
- ✅ Indexed (.bai file present)
- ✅ Contains reads mapping to chromosome 6 / HLA region
- ✅ Proper pairing information preserved
- ✅ Aligned to HLA-aware reference (see above)

### 4. **Coverage Requirements**

Minimum recommended coverage in HLA region:
- **WGS:** 30x average coverage (10,000+ reads in HLA region)
- **RNA-seq:** Highly variable, depends on expression levels
- **Exome:** Usually insufficient for SpecHLA

---

## Reference Genome Requirements

### Obtaining an HLA-Aware Reference

#### Option 1: Download bwakit (Recommended)

```bash
# On Puhti
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
mkdir -p references
cd references

# Download bwakit
wget https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2

# Extract
tar -xjf bwakit-0.7.15_x64-linux.tar.bz2

# The HLA-aware reference is here:
# bwakit/resource-GRCh38/hs38DH.fa

# Copy to a convenient location
cp bwakit/resource-GRCh38/hs38DH.fa ./
cp bwakit/resource-GRCh38/hs38DH.fa.alt ./

# Index if needed
module load biokit
samtools faidx hs38DH.fa
```

#### Option 2: Build Custom HLA Reference

If you need a custom reference:

```bash
# Download GRCh38 with ALT contigs
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

# Or use 1000 Genomes reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

### Checking Your Current BAM Reference

```bash
# Check what reference was used
samtools view -H alignment-sorted.bam | grep "^@SQ" | head -20

# Look for:
# 1. Chromosome naming (chr6 vs 6)
# 2. Presence of HLA-specific contigs (HLA-A*01:01:01:01, etc.)
# 3. ALT contig annotations

# Check for HLA contigs
samtools view -H alignment-sorted.bam | grep -i "HLA"
```

---

## Step-by-Step Workflow

### Step 1: Setup Environment

```bash
# SSH to Puhti
ssh ozcanumu@puhti.csc.fi

# Navigate to your analysis directory
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis

# Create scripts directory if not exists
mkdir -p scripts
cd scripts

# Upload the preparation scripts
# (Upload step3_prepare_bam_input.sh and step4_run_spechla_bam.sh)
```

### Step 2: Prepare BAM Input

```bash
# Make script executable
chmod +x step3_prepare_bam_input.sh

# Submit as SLURM job (recommended)
sbatch step3_prepare_bam_input.sh

# OR run interactively (for testing)
bash step3_prepare_bam_input.sh
```

**What this script does:**
1. ✅ Validates BAM file exists and is accessible
2. ✅ Checks for BAM index, creates if missing
3. ✅ Analyzes reference genome compatibility
4. ✅ Detects chromosome naming convention (chr6 vs 6)
5. ✅ Counts reads in HLA region
6. ✅ Assesses data type (WGS vs exome vs RNA-seq)
7. ✅ Creates sample sheet for Nextflow
8. ✅ Generates configuration file
9. ✅ Provides recommendations based on analysis

**Expected output:**
```
✓ Directory structure created
✓ BAM file exists (size: 15G)
✓ BAM index exists
✓ Chromosome naming: WITH 'chr' prefix (chr6)
✓ Found 25 HLA-related sequences
✓ HLA region has coverage
  Total reads: 45,234,567
  Reads in HLA region: 125,432 (0.28%)

⚠ WARNING: Your BAM lacks HLA ALT contigs
→ SpecHLA accuracy will be LIMITED
```

### Step 3: Review Configuration

```bash
# Check the generated configuration
cat ../pipeline_input_bam/spechla_config.txt

# Review sample sheet
cat ../pipeline_input_bam/samplesheet.csv

# Check disk space
df -h /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
```

### Step 4: Configure Reference Genome

**CRITICAL STEP:** Edit the pipeline script to specify your reference genome.

```bash
# Edit the run script
nano step4_run_spechla_bam.sh

# Find this section (around line 97):
# REF_GENOME="${BASE_DIR}/references/hs38DH.fa"

# Update to your actual reference path
# This MUST be the same reference used to create your BAM!
```

### Step 5: Run SpecHLA Analysis

```bash
# Make script executable
chmod +x step4_run_spechla_bam.sh

# Submit as SLURM job
sbatch step4_run_spechla_bam.sh

# Check job status
squeue -u ozcanumu

# Monitor log in real-time
tail -f ../logs/spechla_bam_*.log
```

### Step 6: Check Results

```bash
# Once completed, view results
cd ../results_bam_spechla/spechla/alignment-sorted/

# Main result file
cat alignment-sorted.spechla.result.txt

# Detailed output
cat alignment-sorted.spechla.detailed.txt

# View HTML report
firefox ../../spechla_report.html &
```

---

## Pipeline Integration

### Nextflow Configuration

Create or update your pipeline configuration:

```bash
# Copy the configuration file to your pipeline
cp spechla_bam.config ${PIPELINE_DIR}/conf/

# Or integrate into main nextflow.config
```

### Modifying Main Pipeline

If you need to modify the main pipeline to support BAM input:

```groovy
// In main.nf

// Add BAM input channel
if (params.input_type == 'bam') {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            tuple(
                row.sample_id,
                file(row.bam),
                file(row.bai)
            )
        }
        .set { ch_bam_input }
} else if (params.input_type == 'fastq') {
    // Existing FASTQ logic
    ...
}

// SpecHLA process
process SPECHLA {
    tag "$sample_id"
    publishDir "${params.outdir}/spechla/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference_genome
    
    output:
    tuple val(sample_id), path("*.result.txt"), emit: results
    tuple val(sample_id), path("*.detailed.txt"), emit: detailed
    path "*.log", emit: logs
    
    script:
    """
    # Extract HLA region (optional but recommended)
    samtools view -b ${bam} ${params.chr6_name}:${params.hla_region_start}-${params.hla_region_end} \
        > ${sample_id}.hla_region.bam
    samtools index ${sample_id}.hla_region.bam
    
    # Run SpecHLA
    spechla typing \
        --bam ${sample_id}.hla_region.bam \
        --reference ${reference_genome} \
        --genes ${params.spechla_genes} \
        --resolution ${params.spechla_resolution} \
        --min-depth ${params.spechla_min_depth} \
        --output ${sample_id}.spechla \
        --threads ${task.cpus} \
        2>&1 | tee ${sample_id}.spechla.log
    
    # Rename outputs
    mv ${sample_id}.spechla.typing.txt ${sample_id}.result.txt
    mv ${sample_id}.spechla.detailed.txt ${sample_id}.detailed.txt
    """
}
```

### Alternative: Create Standalone SpecHLA Module

```bash
# Create module directory
mkdir -p ${PIPELINE_DIR}/modules/spechla_bam

# Create module file
cat > ${PIPELINE_DIR}/modules/spechla_bam/main.nf << 'EOF'
// SpecHLA module for BAM input

process SPECHLA_BAM {
    tag "$sample_id"
    label 'process_high'
    container params.spechla_container
    
    publishDir "${params.outdir}/spechla/${sample_id}", 
        mode: 'copy',
        pattern: "*.{txt,tsv,json}"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path(reference_genome)
    
    output:
    tuple val(sample_id), path("${sample_id}.hla_types.txt"), emit: hla_types
    tuple val(sample_id), path("${sample_id}.allele_depths.tsv"), emit: depths
    path("${sample_id}.spechla.log"), emit: log
    
    script:
    def genes = params.spechla_genes.replaceAll(',', ' ')
    def chr6 = params.chr6_name ?: 'chr6'
    
    """
    # Validate inputs
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Reference: ${reference_genome}"
    
    # Check BAM accessibility
    samtools quickcheck ${bam} || exit 1
    
    # Extract HLA region for efficiency
    echo "Extracting HLA region (${chr6}:28000000-34000000)..."
    samtools view -b -q 20 ${bam} ${chr6}:28000000-34000000 \
        > hla_region.bam
    samtools index hla_region.bam
    
    # Count reads
    READS=\$(samtools view -c hla_region.bam)
    echo "HLA region reads: \$READS"
    
    if [ \$READS -lt 100 ]; then
        echo "ERROR: Insufficient reads in HLA region (\$READS)"
        exit 1
    fi
    
    # Run SpecHLA
    spechla typing \
        --bam hla_region.bam \
        --reference ${reference_genome} \
        --genes ${genes} \
        --resolution ${params.spechla_resolution} \
        --min-depth ${params.spechla_min_depth} \
        --min-mapq ${params.spechla_min_mapq} \
        --min-baseq ${params.spechla_min_baseq} \
        --output ${sample_id} \
        --threads ${task.cpus} \
        2>&1 | tee ${sample_id}.spechla.log
    
    # Organize outputs
    mv ${sample_id}.typing.txt ${sample_id}.hla_types.txt
    mv ${sample_id}.depths.tsv ${sample_id}.allele_depths.tsv || touch ${sample_id}.allele_depths.tsv
    
    echo "SpecHLA analysis complete for ${sample_id}"
    """
}
EOF
```

---

## Troubleshooting

### Issue 1: "No HLA ALT contigs found"

**Cause:** BAM aligned to reference without HLA ALT contigs

**Solutions:**
1. **Best:** Realign your data to HLA-aware reference (hs38DH.fa)
2. **Alternative:** Use OptiType or ArcasHLA instead
3. **Workaround:** Continue with warning, but results may be inaccurate

```bash
# Check current reference
samtools view -H alignment-sorted.bam | grep "^@SQ" | grep -i HLA

# If no HLA contigs, consider realignment:
# 1. Download bwakit reference (see above)
# 2. Realign with BWA
bwa mem -t 40 -K 100000000 \
    references/hs38DH.fa \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    | samtools sort -@ 8 -o realigned.bam -
samtools index realigned.bam
```

### Issue 2: "Insufficient reads in HLA region"

**Cause:** Low coverage or wrong chromosome naming

**Check:**
```bash
# Try different chromosome names
samtools view -c alignment-sorted.bam chr6:28000000-34000000
samtools view -c alignment-sorted.bam 6:28000000-34000000

# Check overall coverage
samtools depth alignment-sorted.bam | awk '{sum+=$3; count++} END {print "Average depth:", sum/count}'
```

**Solutions:**
- If exome data: Use OptiType instead
- If RNA-seq: Use ArcasHLA (designed for transcriptome)
- If wrong chr naming: Update `chr6_name` parameter

### Issue 3: "SpecHLA command not found"

**Cause:** Container not loaded or SpecHLA not installed

**Solutions:**
```bash
# Check if module available
module spider spechla

# Or use Singularity container
singularity exec docker://quay.io/biocontainers/spechla:1.0.0--pyhdfd78af_0 spechla --help

# Update pipeline to use correct container path
```

### Issue 4: "Reference genome not found"

**Cause:** Incorrect path or reference not accessible

**Check:**
```bash
# Verify reference exists
ls -lh /path/to/reference.fa

# Check if indexed
ls -lh /path/to/reference.fa.fai

# Test accessibility in your job
srun --account=project_2008084 --partition=small --time=00:10:00 --mem=4G \
    ls -lh /path/to/reference.fa
```

### Issue 5: "Chromosome naming mismatch"

**Symptoms:**
- 0 reads extracted from HLA region
- SpecHLA returns empty results

**Diagnosis:**
```bash
# Check BAM chromosome names
samtools idxstats alignment-sorted.bam | head

# Check reference chromosome names
grep "^>" reference.fa | head

# Check for mismatch
```

**Fix:**
```bash
# Option 1: Update chr6_name parameter
--chr6_name '6'  # if no chr prefix

# Option 2: Convert BAM chromosome names
# (Not recommended, use correct reference instead)
```

### Issue 6: "Out of memory error"

**Symptoms:**
- Job killed with "Out of Memory"
- SLURM error code 137

**Solutions:**
```bash
# Increase memory in configuration
--max_memory '256.GB'

# Or in SBATCH header
#SBATCH --mem=256G

# Extract HLA region first to reduce memory
samtools view -b alignment-sorted.bam chr6:28000000-34000000 > hla_only.bam
# Then analyze hla_only.bam
```

### Issue 7: "SpecHLA produces no output"

**Check logs:**
```bash
# SLURM log
cat ../logs/spechla_bam_*.log

# Nextflow log
cat ${PIPELINE_DIR}/.nextflow.log

# SpecHLA specific log
find ${PIPELINE_DIR}/work -name "*.spechla.log"
```

**Common causes:**
1. No reads in HLA region → Check coverage
2. Reference mismatch → Verify reference genome
3. Low quality data → Check BAM quality metrics
4. Tool version incompatibility → Update container

### Issue 8: "Pipeline hangs/no progress"

**Check:**
```bash
# Job status
squeue -u ozcanumu

# Nextflow status
cd ${PIPELINE_DIR}
nextflow log

# Check work directories
ls -ltrh work/*/*
```

**Resume:**
```bash
# Nextflow can resume from last successful step
sbatch step4_run_spechla_bam.sh
# -resume is already in the script
```

---

## Known Limitations

### SpecHLA Limitations

1. **Reference Dependency:** 
   - Absolutely requires HLA-aware reference
   - Cannot function properly with standard references
   - Results are only as good as reference quality

2. **Data Type Restrictions:**
   - **Exome data:** Poor performance, not recommended
   - **Targeted panels:** Only if includes HLA region
   - **Low coverage WGS:** May have insufficient depth

3. **Resolution Limitations:**
   - Best with high-resolution alleles (4-field)
   - Lower resolution (2-field) more reliable but less specific
   - Ambiguous alleles may not resolve completely

4. **Computational Requirements:**
   - High memory usage for WGS
   - Slower than read-based methods (OptiType)
   - Requires significant disk I/O

### RNA-seq Specific Limitations

1. **Expression Variability:**
   - HLA genes have different expression levels
   - Some alleles may not be detected if lowly expressed
   - Cannot distinguish unexpressed alleles

2. **Alternative Splicing:**
   - May complicate allele calling
   - Intronic regions not captured

3. **Recommended Alternatives for RNA-seq:**
   - **OptiType (RNA mode):** Primary recommendation
   - **ArcasHLA:** Designed for RNA-seq
   - **Both together:** Best for validation

### General Limitations

1. **Novel Alleles:**
   - Cannot detect alleles not in IMGT database
   - May misclassify rare alleles

2. **Copy Number Variants:**
   - Assumes diploid HLA genes
   - May struggle with deletions/duplications

3. **Contamination:**
   - Sensitive to sample contamination
   - Multiple genotypes will confuse calling

4. **Homozygosity:**
   - May report same allele twice
   - Can be hard to distinguish from heterozygous with one allele

---

## Best Practices

### Before Analysis

1. ✅ **Verify reference genome compatibility**
   - Check for HLA ALT contigs
   - Confirm chromosome naming
   - Validate reference is accessible

2. ✅ **Assess data suitability**
   - Check coverage in HLA region (>10,000 reads minimum)
   - Verify data type (WGS preferred)
   - Consider alternative tools for exome/RNA-seq

3. ✅ **Validate BAM file**
   - Confirm proper sorting
   - Check index exists
   - Verify mapping quality

### During Analysis

1. ✅ **Monitor resource usage**
   - Check memory consumption
   - Track disk space
   - Watch job queue time

2. ✅ **Review intermediate outputs**
   - Check HLA read extraction worked
   - Verify reads mapped to HLA region
   - Confirm no obvious errors in logs

3. ✅ **Use appropriate parameters**
   - 2-field resolution for standard analysis
   - Higher min_depth for stricter calls
   - Adjust based on coverage

### After Analysis

1. ✅ **Validate results**
   - Check concordance between alleles
   - Verify expected genotype pattern
   - Compare with other methods if available

2. ✅ **Document settings**
   - Record reference genome used
   - Note parameter choices
   - Save relevant logs

3. ✅ **Consider orthogonal validation**
   - Run OptiType on same sample
   - Compare with ArcasHLA
   - Use majority voting if multiple tools

---

## Comparison with Other Tools

| Feature | SpecHLA | OptiType | ArcasHLA |
|---------|---------|----------|----------|
| **Best for** | WGS | Exome/RNA | RNA-seq |
| **Reference requirement** | HLA-aware | Standard | Standard |
| **Speed** | Slow | Fast | Medium |
| **Accuracy (WGS)** | High* | High | Medium |
| **Accuracy (RNA)** | Low | High | High |
| **Class I (A,B,C)** | ✅ | ✅ | ✅ |
| **Class II (DR,DQ,DP)** | ✅ | ❌ | ✅ |
| **BAM input** | ✅ | ✅ | ✅ |
| **FASTQ input** | ❌ | ✅ | ✅ |

*High accuracy requires proper HLA-aware reference

**Recommendation for your RNA-seq data:**
1. **Primary:** OptiType (RNA mode)
2. **Secondary:** ArcasHLA
3. **Optional:** SpecHLA only if proper reference used

---

## Additional Resources

### Documentation
- SpecHLA GitHub: https://github.com/deepomicslab/SpecHLA
- OptiType: https://github.com/FRED-2/OptiType
- ArcasHLA: https://github.com/RabadanLab/arcasHLA

### References
- bwakit: https://github.com/lh3/bwa/tree/master/bwakit
- IMGT/HLA: https://www.ebi.ac.uk/ipd/imgt/hla/
- GRCh38 with ALT: https://www.ncbi.nlm.nih.gov/genome/guide/human/

### CSC Puhti
- Documentation: https://docs.csc.fi/computing/systems-puhti/
- Bioinformatics guide: https://docs.csc.fi/apps/bioinformatics/
- Contact: servicedesk@csc.fi

---

## Summary Checklist

Before running SpecHLA on your BAM:

- [ ] BAM aligned to HLA-aware reference (hs38DH.fa or similar)
- [ ] BAM is coordinate-sorted and indexed
- [ ] HLA region has sufficient coverage (>1,000 reads minimum)
- [ ] Chromosome naming determined (chr6 vs 6)
- [ ] Reference genome file is accessible
- [ ] Data type is suitable (WGS best, RNA-seq limited)
- [ ] Alternative tools considered (OptiType/ArcasHLA for RNA-seq)
- [ ] Sufficient disk space available (>50GB free)
- [ ] SLURM account and partition configured
- [ ] Scripts uploaded and configured

**For RNA-seq data (your case):**
- [ ] Strongly consider OptiType + ArcasHLA instead
- [ ] If using SpecHLA, understand limitations
- [ ] Plan to validate with other tools
- [ ] Document that results are from RNA-seq

---

## Quick Start Commands

```bash
# Complete workflow in one go
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Step 1: Prepare input
sbatch step3_prepare_bam_input.sh

# Wait for completion, check output
tail ../logs/prepare_bam_*.log

# Step 2: Configure reference genome
nano step4_run_spechla_bam.sh
# Edit REF_GENOME variable

# Step 3: Run analysis
sbatch step4_run_spechla_bam.sh

# Step 4: Monitor progress
squeue -u ozcanumu
tail -f ../logs/spechla_bam_*.log

# Step 5: Check results
cat ../results_bam_spechla/spechla/alignment-sorted/alignment-sorted.spechla.result.txt
```

---

**Document Version:** 1.0  
**Last Updated:** December 2024  
**Author:** Umut Ozcan  
**Project:** project_2008084  
**Pipeline:** HLA Typing Pipeline v2.0.0+
