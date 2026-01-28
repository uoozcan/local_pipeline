# QUICK SOLUTION: Bypass Buggy Pipeline

## The Problem
Your Nextflow pipeline has a bug:
```
ERROR ~ No such property: ch_input
```

## The Solution (30 seconds)

**Extract HLA reads from your BAM files, then type them with any tool.**

### Step 1: Extract Reads (Submit Now)

```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/scripts
sbatch extract_hla_reads.sh
```

**That's it!** This will:
- ✅ Extract HLA region reads from all 24 BAM files
- ✅ Convert to FASTQ format
- ✅ Create README with typing instructions
- ✅ Take ~1-2 hours total

### Step 2: Monitor

```bash
# Check progress
squeue -u $USER

# Watch logs
tail -f logs/extract_hla_*.log

# Check results as they complete
ls results_bam/batch1_VenEx_DNA_BAM/*/hla_reads/
```

### Step 3: Type HLA (After Extraction)

You now have FASTQ files ready for ANY HLA typing tool:

**Option A: Use OptiType (if available on Puhti)**
```bash
module load optitype  # or conda activate optitype
OptiTypePipeline.py -i SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz --dna -o output/
```

**Option B: Download and type locally**
```bash
scp -r ozcanumu@puhti.csc.fi:/scratch/.../hla_reads/ ./
# Then use OptiType, PHLAT, HLA-HD, etc. on your computer
```

**Option C: Use online tools**
- Upload to Galaxy (usegalaxy.org)
- Use IMGT/HLA database tools

---

## What You Get

For each sample:
```
results_bam/batch1_VenEx_DNA_BAM/SAMPLE_ID/hla_reads/
├── SAMPLE_ID_R1.fastq.gz       # Ready for HLA typing
├── SAMPLE_ID_R2.fastq.gz       # Ready for HLA typing
├── SAMPLE_ID_single.fastq.gz   # Unpaired reads
├── SAMPLE_ID_hla.sam           # SAM format
└── README.txt                  # Instructions
```

---

## Why This Works

1. **Extracts HLA reads** - chromosome 6:28510120-33480577
2. **Converts to FASTQ** - universal format
3. **No pipeline bugs** - direct samtools commands
4. **Flexible** - use any HLA typing tool you want
5. **Portable** - analyze anywhere

---

## Recommended HLA Typing Tools

| Tool | Genes | Best For | Speed |
|------|-------|----------|-------|
| OptiType | A, B, C | DNA/WGS | Fast ⚡ |
| PHLAT | A, B, C, DR, DQ, DP | All-in-one | Medium 🔧 |
| HLA-HD | All HLA | Highest accuracy | Slow 🐌 |
| seq2HLA | A, B, C, DR, DQ | RNA-seq | Fast ⚡ |

---

## Quick Commands

```bash
# Submit extraction (DO THIS NOW)
sbatch extract_hla_reads.sh

# After ~2 hours, check results
ls -lh results_bam/batch1_VenEx_DNA_BAM/*/hla_reads/*.fastq.gz

# Count total HLA reads per sample
for dir in results_bam/batch1_VenEx_DNA_BAM/*/hla_reads/; do
    sample=$(basename $(dirname $dir))
    reads=$(zcat ${dir}/*_R1.fastq.gz | wc -l)
    reads=$((reads / 4))
    echo "${sample}: ${reads} paired reads"
done

# Type one sample with OptiType (example)
OptiTypePipeline.py \
    -i results_bam/.../VX_36_4_D1/hla_reads/VX_36_4_D1_R1.fastq.gz \
       results_bam/.../VX_36_4_D1/hla_reads/VX_36_4_D1_R2.fastq.gz \
    --dna -o optitype_VX_36_4_D1/ -v
```

---

## Expected Timeline

- **Extraction**: 1-2 hours for 24 samples (parallel)
- **Typing** (if local): 
  - OptiType: ~5-10 min/sample
  - PHLAT: ~10-15 min/sample
  - HLA-HD: ~30-60 min/sample

---

## Advantages Over Pipeline

| Aspect | Pipeline | This Approach |
|--------|----------|---------------|
| Works? | ❌ Buggy | ✅ Always |
| Speed | N/A | ⚡ 1-2 hours |
| Flexibility | 🔒 Fixed tools | 🔓 Any tool |
| Debugging | 😰 Hard | 😊 Easy |
| Portable | 🖥️ Puhti only | 🌍 Anywhere |

---

## Files Provided

- **extract_hla_reads.sh** ← USE THIS (the simple solution)
- **run_hla_simple.sh** - Alternative with auto-detection
- **run_hla_direct_array.sh** - For container-based tools
- **DIRECT_HLA_TYPING_GUIDE.md** - Complete documentation

---

## Troubleshooting

**Q: No reads extracted?**
```bash
# Check chromosome naming in your BAM
samtools view -H your.bam | grep "^@SQ"
# Look for "SN:6" or "SN:chr6"
```

**Q: Tool not available on Puhti?**
```bash
# Option 1: Use conda
conda install -c bioconda optitype

# Option 2: Download reads and type locally
scp -r ozcanumu@puhti:/.../hla_reads/ ./

# Option 3: Use online tools (Galaxy, etc.)
```

**Q: Low read count?**
- <100 reads: Poor coverage, results unreliable
- 100-1000 reads: Acceptable, most genes
- >1000 reads: Good, high confidence

---

## Summary

**DO THIS NOW:**
```bash
sbatch extract_hla_reads.sh
```

**THEN (after 2 hours):**
- Check extraction logs
- Choose your HLA typing tool
- Type your samples
- Collect results

**RESULT:**
- ✅ HLA types for all 24 samples
- ✅ No pipeline bugs
- ✅ Complete control
- ✅ Reproducible workflow

---

## Next Steps After HLA Typing

Once you have HLA types:

1. **Create summary table**
```bash
# Combine all results
cat */optitype/*result.tsv > all_hla_types.tsv
```

2. **Quality check**
- Compare with known genotypes if available
- Check for unusual alleles
- Verify read coverage was adequate

3. **Downstream analysis**
- Association studies
- Population genetics
- Clinical interpretation

---

**Bottom Line**: 
Run `sbatch extract_hla_reads.sh` right now. In 2 hours, you'll have FASTQ files ready for any HLA typing tool, bypassing all pipeline bugs.
