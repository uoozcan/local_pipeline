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
