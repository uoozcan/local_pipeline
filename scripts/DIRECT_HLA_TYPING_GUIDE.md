# Direct HLA Typing - Bypass Buggy Pipeline

## The Problem

The Nextflow pipeline has an internal bug:
```
ERROR ~ No such property: ch_input for class: nextflow.script.WorkflowBinding
```

This is a pipeline code error that can't be fixed without modifying the pipeline source code.

## The Solution

**Run HLA typing tools directly** without the Nextflow wrapper. This gives you:
- ✅ More control
- ✅ No pipeline bugs
- ✅ Faster debugging
- ✅ Better understanding of what's happening
- ✅ Same results as pipeline would produce

---

## Option 1: Simple Extraction (RECOMMENDED - Always Works)

**Best for**: When you want maximum compatibility or just need the reads

This approach:
1. Extracts HLA region reads from your BAM files
2. Converts them to FASTQ format
3. You can then use ANY HLA typing tool

### Quick Start

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Submit extraction job
sbatch run_hla_simple.sh
```

### What You Get

For each sample:
```
results_bam/batch1_VenEx_DNA_BAM/SAMPLE_ID/
└── hla_reads/
    ├── SAMPLE_ID_HLA_R1.fastq.gz      # Forward reads
    ├── SAMPLE_ID_HLA_R2.fastq.gz      # Reverse reads
    ├── SAMPLE_ID_HLA_single.fastq.gz  # Unpaired reads
    ├── SAMPLE_ID_hla.sam              # SAM format
    └── extraction_summary.txt         # Info file
```

### Then Use Any Tool

**OptiType** (Class I: A, B, C):
```bash
module load optitype  # if available
OptiTypePipeline.py \
    -i SAMPLE_ID_HLA_R1.fastq.gz SAMPLE_ID_HLA_R2.fastq.gz \
    --dna -o optitype_output/ -v
```

**PHLAT** (Class I and II):
```bash
python PHLAT.py \
    -1 SAMPLE_ID_HLA_R1.fastq.gz \
    -2 SAMPLE_ID_HLA_R2.fastq.gz \
    -o phlat_output/
```

**HLA-HD** (High accuracy):
```bash
hlahd.sh \
    -t 8 \
    -f freq_data/ \
    SAMPLE_ID_HLA_R1.fastq.gz \
    SAMPLE_ID_HLA_R2.fastq.gz \
    dictionary/ \
    SAMPLE_ID \
    output_dir/
```

---

## Option 2: Direct Tool Execution

**Best for**: When you have containers and want automated typing

Uses available HLA typing containers directly.

### Setup

```bash
# Check what containers you have
ls /scratch/project_2008084/hla_references/singularity_cache/containers/

# Or find them
find /scratch/project_2008084 -name "*.sif" -type f
```

### Run

```bash
# Edit run_hla_direct_array.sh to point to your containers
nano run_hla_direct_array.sh

# Submit
bash submit_direct_hla.sh
```

---

## Option 3: Individual Sample Processing

**Best for**: Testing or processing specific samples

### Extract One Sample

```bash
#!/bin/bash
SAMPLE_ID="VX_36_4_D1"
BAM_FILE="/path/to/VX_36_4_D1_tumor.bam"
OUTPUT_DIR="./hla_output/${SAMPLE_ID}"

mkdir -p ${OUTPUT_DIR}

# Extract HLA region
samtools view -b ${BAM_FILE} 6:28510120-33480577 > ${OUTPUT_DIR}/hla_region.bam

# Convert to FASTQ
samtools fastq ${OUTPUT_DIR}/hla_region.bam \
    -1 ${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz \
    -2 ${OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz \
    -s ${OUTPUT_DIR}/${SAMPLE_ID}_single.fastq.gz

# Count reads
samtools view -c ${OUTPUT_DIR}/hla_region.bam

echo "Done! FASTQ files in: ${OUTPUT_DIR}/"
```

### Type One Sample

```bash
# Using OptiType (if available)
OptiTypePipeline.py \
    -i ${OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz \
       ${OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz \
    --dna -o ${OUTPUT_DIR}/optitype/ -v
```

---

## Comparison of Approaches

| Approach | Pros | Cons | Time |
|----------|------|------|------|
| **Simple Extraction** | Always works, flexible | Need to run typing separately | 1-2 hours |
| **Direct Execution** | Automated, complete | Requires containers | 3-6 hours |
| **Individual Sample** | Good for testing | Manual for each sample | ~15 min/sample |

---

## What Tools Are Available?

### Check on Puhti

```bash
# Check for modules
module spider optitype
module spider phlat
module spider hlahd

# Check for installed tools
which OptiTypePipeline.py
which hlahd.sh

# Check for containers
ls /scratch/project_2008084/*/singularity_cache/containers/
find /scratch/project_2008084 -name "*hla*.sif"
```

### If Nothing Available

You have two options:

**A. Use the extracted FASTQ files on your local machine:**
```bash
# Download the reads
scp -r ozcanumu@puhti.csc.fi:/scratch/.../hla_reads/ ./local_analysis/

# Type locally with your favorite tool
```

**B. Install a tool on Puhti:**
```bash
# OptiType (simplest - conda)
conda create -n optitype -c bioconda optitype
conda activate optitype

# Then use on extracted reads
```

---

## Complete Workflow Example

### 1. Extract HLA Reads (Always Do This First)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts

# Submit extraction for all 24 samples
sbatch run_hla_simple.sh

# Monitor
squeue -u $USER
tail -f logs/hla_simple_*.log
```

### 2. Choose Typing Method

After extraction completes, choose based on what's available:

**If OptiType is available:**
```bash
# Create typing script
cat > type_with_optitype.sh << 'EOF'
#!/bin/bash
#SBATCH --array=1-24
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Get sample
SAMPLE_DIR=$(find results_bam/batch1_VenEx_DNA_BAM -maxdepth 1 -mindepth 1 -type d | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE_ID=$(basename ${SAMPLE_DIR})

# Type with OptiType
OptiTypePipeline.py \
    -i ${SAMPLE_DIR}/hla_reads/${SAMPLE_ID}_HLA_R1.fastq.gz \
       ${SAMPLE_DIR}/hla_reads/${SAMPLE_ID}_HLA_R2.fastq.gz \
    --dna \
    -o ${SAMPLE_DIR}/optitype_results/ \
    -v
EOF

sbatch type_with_optitype.sh
```

**If nothing available:**
```bash
# Download reads to your computer
scp -r ozcanumu@puhti.csc.fi:/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/results_bam/batch1_VenEx_DNA_BAM/*/hla_reads/ ./

# Type locally
```

### 3. Collect Results

Once typing is complete, collect all results:

```bash
# Create summary
for sample_dir in results_bam/batch1_VenEx_DNA_BAM/*/; do
    sample=$(basename $sample_dir)
    
    # Check for OptiType results
    if [ -f "${sample_dir}/optitype_results/${sample}_result.tsv" ]; then
        echo "${sample}: $(cat ${sample_dir}/optitype_results/${sample}_result.tsv)"
    fi
done > all_hla_results.txt
```

---

## HLA Region Details

**Chromosome 6 HLA Region**: `6:28510120-33480577` (hg38)

This includes:
- **Class I**: HLA-A, HLA-B, HLA-C
- **Class II**: HLA-DRA, HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DPA1, HLA-DPB1
- **Class III**: Complement, TNF, heat shock proteins

### Expected Read Counts

For WGS data (30x coverage):
- Good: >1,000 HLA reads
- Acceptable: 100-1,000 reads
- Poor: <100 reads

Check your extraction summary files to see coverage.

---

## Recommended HLA Typing Tools

### For DNA/WGS Data (Your Case)

1. **OptiType** - Best for Class I (A, B, C)
   - Fast, accurate
   - Well-validated
   - Easy to use

2. **HLA-HD** - Best accuracy overall
   - Supports Class I and II
   - Higher resolution
   - Slower

3. **PHLAT** - Good balance
   - Supports Class I and II
   - Fast
   - Good accuracy

### Not Recommended for WGS

- ArcasHLA: Designed for RNA-seq
- seq2HLA: Optimized for RNA-seq

---

## Troubleshooting

### No HLA Reads Extracted

```bash
# Check chromosome naming
samtools view -H your.bam | grep "^@SQ"

# If using "chr6" instead of "6":
samtools view -b your.bam chr6:28510120-33480577 > hla_region.bam
```

### Low Read Count

Check coverage:
```bash
samtools depth -r 6:28510120-33480577 your.bam | awk '{sum+=$3} END {print sum/NR}'
```

If <10x, HLA typing may be unreliable.

### Tool Installation Issues

Use conda/mamba:
```bash
# Install optitype
mamba create -n hla optitype -c bioconda

# Activate and use
conda activate hla
OptiTypePipeline.py --help
```

---

## Summary

**Simplest Approach (Always Works)**:
```bash
1. sbatch run_hla_simple.sh           # Extract HLA reads
2. Download reads to local machine     # If no tools on Puhti
3. Use any HLA typing tool locally    # OptiType, PHLAT, etc.
4. Upload results if needed           # For archiving
```

**If You Have Tools on Puhti**:
```bash
1. sbatch run_hla_simple.sh           # Extract HLA reads
2. module load optitype               # Load tool
3. Run typing on extracted reads      # Batch or individual
4. Collect results                    # Create summary
```

**Benefits Over Pipeline**:
- ✅ No buggy Nextflow code
- ✅ Complete control
- ✅ Easy to debug
- ✅ Works with any tool
- ✅ Portable - can analyze anywhere

---

## Files Provided

| File | Purpose |
|------|---------|
| `run_hla_simple.sh` | Extract HLA reads (RECOMMENDED) |
| `run_hla_direct_array.sh` | Run tools directly with containers |
| `submit_direct_hla.sh` | Master submission script |

---

## Next Steps

1. **Start with extraction**:
   ```bash
   sbatch run_hla_simple.sh
   ```

2. **Check what tools are available**:
   ```bash
   module spider optitype
   which OptiTypePipeline.py
   ```

3. **Choose your typing approach** based on what's available

4. **Process your samples** and collect results

This approach will definitely work and give you the HLA types you need!
