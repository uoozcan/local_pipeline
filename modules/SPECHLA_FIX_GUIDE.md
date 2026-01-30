# SpecHLA Module Syntax Error - Diagnosis and Fix

## Error Summary

**Error Type:** Nextflow Module Compilation Error  
**Location:** `/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/modules/spechla.nf`  
**Line:** 7, Column 21  
**Message:** `Unexpected input: '{' @ line 7, column 21`

## Root Cause

The error "Unexpected input: '{'" at a process declaration indicates a **Nextflow DSL2 syntax error**. Common causes:

1. **Missing closing brace** from a previous process or workflow block
2. **Incorrect process directive syntax** before the opening brace
3. **Malformed tag, label, or publishDir directives**
4. **Missing or extra commas** in directive lists
5. **Improper use of quotes** in directive values

The error specifically occurs at:
```nextflow
process SPECHLA_BAM {
                    ^
                    Line 7, Column 21
```

This suggests the parser encounters an unexpected opening brace, likely because:
- Something before line 7 is malformed
- The process declaration itself is missing required elements
- There's an unclosed statement from previous lines

## Fix Strategy

### Quick Fix (Automated)

1. **Upload and run the fix script:**
   ```bash
   cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis
   
   # Upload fix_spechla_syntax.sh to this directory
   
   chmod +x fix_spechla_syntax.sh
   ./fix_spechla_syntax.sh
   ```

2. **Verify the fix:**
   ```bash
   nextflow run main.nf -resume --help
   ```

### Manual Fix

If you prefer to fix manually or want to understand what's wrong:

1. **View the problematic file:**
   ```bash
   cat -n modules/spechla.nf | head -20
   ```

2. **Check for common issues:**
   
   **Issue A: Missing closing brace from previous section**
   ```nextflow
   // WRONG - previous section missing closing brace
   process SOME_PREVIOUS_PROCESS {
       // ... process content
   // Missing closing brace here!
   
   process SPECHLA_BAM {  // ERROR: unexpected '{'
   ```
   
   **Fix:** Add missing closing brace before SPECHLA_BAM
   
   **Issue B: Incorrect directive syntax**
   ```nextflow
   // WRONG - malformed directives
   process SPECHLA_BAM
   tag "$sample_id"  // Missing: at beginning or inside { }
   {
   ```
   
   **Fix:** Move directives inside the process block:
   ```nextflow
   process SPECHLA_BAM {
       tag "$sample_id"
       label 'process_high'
   ```
   
   **Issue C: Missing DSL2 declaration**
   ```nextflow
   // WRONG - missing at top of file
   process SPECHLA_BAM {
   ```
   
   **Fix:** Not applicable here (DSL2 is set in main.nf)

3. **Compare with working module structure:**
   ```bash
   # View a working module for reference
   cat modules/optitype.nf | head -30
   cat modules/arcashla.nf | head -30
   ```

## Correct SpecHLA Module Structure

The corrected `spechla.nf` should follow this structure:

```nextflow
/*
========================================================================================
    SpecHLA Module - HLA Typing from DNA/RNA Sequencing Data
========================================================================================
*/

process SPECHLA_BAM {
    tag "$sample_id"
    label 'process_high'
    
    publishDir "${params.outdir}/${sample_id}/spechla", mode: 'copy'
    
    container "${params.singularity_cache_dir}/spechla.sif"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_spechla"), emit: results
    path "${sample_id}_spechla/logs/*.log", emit: logs, optional: true
    
    script:
    """
    # Process script here
    """
}

process SPECHLA_FASTQ {
    // Similar structure
}
```

### Key Elements (All Required)

1. **Process declaration:** `process PROCESS_NAME {`
2. **Tag directive:** `tag "$sample_id"` (for job naming)
3. **Label directive:** `label 'process_high'` (for resource allocation)
4. **publishDir directive:** Where to save results
5. **Container directive:** Singularity/Docker container path
6. **Input block:** `input:` with tuple/path declarations
7. **Output block:** `output:` with emit channels
8. **Script block:** `script:` or `script:` with triple-quoted shell commands
9. **Closing brace:** `}` to close the process

## Verification Steps

After applying the fix:

### 1. Check Nextflow Syntax
```bash
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis

# Test pipeline parsing (should show no compilation errors)
nextflow inspect main.nf
```

Expected output: Pipeline DAG structure (no errors)

### 2. Dry Run Test
```bash
# Test without actually running jobs
nextflow run main.nf \
    --input /path/to/test/sample.bam \
    --input_type bam \
    --tools optitype,arcashla,spechla \
    -profile singularity \
    -dry-run
```

### 3. Resume Your Batch Job
```bash
# Resume the failed pipeline run
sbatch scripts/step3_run_pipeline_bam.sh
```

## Troubleshooting

### If Error Persists

1. **Check all module files for syntax:**
   ```bash
   for module in modules/*.nf; do
       echo "Checking: $module"
       nextflow inspect main.nf 2>&1 | grep -i error || echo "✓ OK"
   done
   ```

2. **Look for unclosed braces:**
   ```bash
   # Count braces in each module
   for module in modules/*.nf; do
       OPEN=$(grep -o '{' "$module" | wc -l)
       CLOSE=$(grep -o '}' "$module" | wc -l)
       echo "$module: Open=$OPEN, Close=$CLOSE"
       if [ $OPEN -ne $CLOSE ]; then
           echo "  ⚠ MISMATCH!"
       fi
   done
   ```

3. **View Nextflow's detailed error log:**
   ```bash
   cat .nextflow.log | tail -50
   ```

### Other Common Module Issues

**Issue: Container not found**
```
ERROR: Container image not found
```
Fix: Verify container path in nextflow.config:
```bash
ls -lh /scratch/project_2008084/hla_references/singularity_cache/containers/spechla.sif
```

**Issue: Reference genome not found**
```
ERROR: Reference file not found
```
Fix: Update reference path in config_bam_batch.sh:
```bash
export REFERENCE_GENOME="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/references/hs38DH.fa"
```

## Additional Resources

### Files Provided

1. **fix_spechla_syntax.sh** - Automated fix script
2. **spechla_corrected.nf** - Complete corrected module
3. **diagnose_spechla.sh** - Diagnostic script for analysis

### Running the Fix

```bash
# On CSC Puhti
cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis

# If files were uploaded to a different location, copy them:
# cp /path/to/uploaded/fix_spechla_syntax.sh .
# cp /path/to/uploaded/spechla_corrected.nf .

# Run the fix
chmod +x fix_spechla_syntax.sh
./fix_spechla_syntax.sh

# Verify
nextflow inspect main.nf

# Resume pipeline
sbatch scripts/step3_run_pipeline_bam.sh
```

## Summary

**Problem:** Syntax error in `spechla.nf` at process declaration  
**Cause:** Malformed Nextflow DSL2 process definition  
**Solution:** Replace with properly structured module file  
**Verification:** Run `nextflow inspect main.nf` - should show no errors  
**Next Step:** Resume pipeline with `sbatch scripts/step3_run_pipeline_bam.sh`

The corrected module follows proper Nextflow DSL2 syntax with:
- Correct process declaration
- Proper directive formatting  
- Complete input/output specifications
- Robust error handling in script blocks
- Compatibility with OptiType, ArcasHLA, and SpecHLA workflow

---

**Need Help?** Check the Nextflow log: `.nextflow.log`  
**Still Issues?** Verify other modules: `modules/optitype.nf`, `modules/arcashla.nf`
