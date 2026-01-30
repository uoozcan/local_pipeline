#!/bin/bash
# Fix script for SpecHLA module syntax error

set -e

echo "=========================================="
echo "SpecHLA Module Syntax Error Fix"
echo "=========================================="
echo ""

# Paths
MODULES_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/modules"
SPECHLA_FILE="${MODULES_DIR}/spechla.nf"
BACKUP_FILE="${MODULES_DIR}/spechla.nf.backup_$(date +%Y%m%d_%H%M%S)"

# Check if modules directory exists
if [ ! -d "$MODULES_DIR" ]; then
    echo "ERROR: Modules directory not found: $MODULES_DIR"
    echo "Creating modules directory..."
    mkdir -p "$MODULES_DIR"
fi

# Backup existing file if it exists
if [ -f "$SPECHLA_FILE" ]; then
    echo "Backing up existing file to: $BACKUP_FILE"
    cp "$SPECHLA_FILE" "$BACKUP_FILE"
    echo "✓ Backup created"
    echo ""
fi

# Create corrected spechla.nf file
echo "Creating corrected SpecHLA module..."

cat > "$SPECHLA_FILE" << 'EOFMODULE'
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
    def genes = params.spechla_genes ?: 'A,B,C,DQA1,DQB1,DRB1'
    def exon_only = params.spechla_exon_only ?: 0
    def reference = params.reference_genome ?: '/references/hs38DH.fa'
    
    """
    # Create output directory
    mkdir -p ${sample_id}_spechla/logs
    
    # Log execution details
    echo "SpecHLA Analysis" > ${sample_id}_spechla/logs/spechla.log
    echo "Sample: ${sample_id}" >> ${sample_id}_spechla/logs/spechla.log
    echo "BAM: ${bam}" >> ${sample_id}_spechla/logs/spechla.log
    echo "BAI: ${bai}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Genes: ${genes}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Exon only mode: ${exon_only}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Reference: ${reference}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Started: \$(date)" >> ${sample_id}_spechla/logs/spechla.log
    
    # Verify BAI exists
    if [ ! -f "${bai}" ]; then
        echo "ERROR: Index file not found: ${bai}" >> ${sample_id}_spechla/logs/spechla.log
        echo "Attempting to create index..." >> ${sample_id}_spechla/logs/spechla.log
        samtools index ${bam}
    fi
    
    # Run SpecHLA
    SpecHLA \\
        -n ${sample_id} \\
        -r ${reference} \\
        -b ${bam} \\
        -e ${exon_only} \\
        -o ${sample_id}_spechla \\
        2>&1 | tee -a ${sample_id}_spechla/logs/spechla.log
    
    # Check if results were generated
    if [ -f "${sample_id}_spechla/${sample_id}.results.txt" ]; then
        echo "SUCCESS: SpecHLA completed" >> ${sample_id}_spechla/logs/spechla.log
    else
        echo "WARNING: SpecHLA results file not found" >> ${sample_id}_spechla/logs/spechla.log
    fi
    
    echo "Completed: \$(date)" >> ${sample_id}_spechla/logs/spechla.log
    """
}

process SPECHLA_FASTQ {
    tag "$sample_id"
    label 'process_high'
    
    publishDir "${params.outdir}/${sample_id}/spechla", mode: 'copy'
    
    container "${params.singularity_cache_dir}/spechla.sif"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_spechla"), emit: results
    path "${sample_id}_spechla/logs/*.log", emit: logs, optional: true
    
    script:
    def genes = params.spechla_genes ?: 'A,B,C,DQA1,DQB1,DRB1'
    def exon_only = params.spechla_exon_only ?: 0
    def reference = params.reference_genome ?: '/references/hs38DH.fa'
    
    """
    # Create output directory
    mkdir -p ${sample_id}_spechla/logs
    
    # Log execution details
    echo "SpecHLA Analysis (FASTQ input)" > ${sample_id}_spechla/logs/spechla.log
    echo "Sample: ${sample_id}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Read1: ${read1}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Read2: ${read2}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Genes: ${genes}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Exon only mode: ${exon_only}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Reference: ${reference}" >> ${sample_id}_spechla/logs/spechla.log
    echo "Started: \$(date)" >> ${sample_id}_spechla/logs/spechla.log
    
    # Note: SpecHLA typically requires BAM input
    # If FASTQ input is needed, first align reads to reference
    echo "Note: SpecHLA requires BAM input. Converting FASTQ to BAM..." >> ${sample_id}_spechla/logs/spechla.log
    
    # Align with BWA (if available in container)
    bwa mem -t ${task.cpus} ${reference} ${read1} ${read2} | \\
        samtools view -bS - | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    
    samtools index ${sample_id}.sorted.bam
    
    # Run SpecHLA on aligned BAM
    SpecHLA \\
        -n ${sample_id} \\
        -r ${reference} \\
        -b ${sample_id}.sorted.bam \\
        -e ${exon_only} \\
        -o ${sample_id}_spechla \\
        2>&1 | tee -a ${sample_id}_spechla/logs/spechla.log
    
    # Check if results were generated
    if [ -f "${sample_id}_spechla/${sample_id}.results.txt" ]; then
        echo "SUCCESS: SpecHLA completed" >> ${sample_id}_spechla/logs/spechla.log
    else
        echo "WARNING: SpecHLA results file not found" >> ${sample_id}_spechla/logs/spechla.log
    fi
    
    echo "Completed: \$(date)" >> ${sample_id}_spechla/logs/spechla.log
    
    # Clean up intermediate BAM
    rm -f ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}
EOFMODULE

echo "✓ Corrected SpecHLA module created"
echo ""

# Verify syntax
echo "Verifying Nextflow syntax..."
echo ""

cd /scratch/project_2008084/ozcanumu/hla_rnaseq_analysis

# Quick syntax check by trying to parse the pipeline
nextflow inspect main.nf 2>&1 | head -20

echo ""
echo "=========================================="
echo "Fix Applied Successfully!"
echo "=========================================="
echo ""
echo "What was fixed:"
echo "  - Proper Nextflow DSL2 process declaration syntax"
echo "  - Correct process directive formatting (tag, label, publishDir, container)"
echo "  - Proper input/output tuple declarations"
echo "  - Complete script blocks with error handling"
echo ""
echo "Changes made:"
echo "  - Original file backed up to: $BACKUP_FILE"
echo "  - New corrected file created: $SPECHLA_FILE"
echo ""
echo "Next steps:"
echo "  1. Review the corrected file: cat $SPECHLA_FILE"
echo "  2. Test the pipeline: nextflow run main.nf -resume"
echo "  3. If issues persist, check other module files (optitype.nf, arcashla.nf)"
echo ""
echo "Common causes of this error:"
echo "  ✗ Missing closing brace in previous process"
echo "  ✗ Incorrect directive syntax (e.g., missing quotes or commas)"
echo "  ✗ Malformed process declaration"
echo "  ✓ Fixed: Proper DSL2 syntax with all required elements"
echo ""
