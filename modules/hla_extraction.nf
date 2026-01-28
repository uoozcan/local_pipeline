// modules/hla_extraction.nf
// Extract HLA-specific reads from FASTQ files

process EXTRACT_HLA_READS {
    tag "$sample_id"
    label 'process_medium'
    
    container "${params.container_cache}/samtools_latest.sif"
    
    publishDir "${params.outdir}/hla_extracted", mode: 'copy', pattern: "*.log"
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_HLA_R1.fastq.gz"), path("${sample_id}_HLA_R2.fastq.gz"), emit: hla_fastq
    path "${sample_id}_extraction.log", emit: log
    
    script:
    // HLA region coordinates for hg38
    def hla_chr6_start = 28477797
    def hla_chr6_end = 33448354
    
    """
    #!/bin/bash
    set -eo pipefail
    
    LOG="${sample_id}_extraction.log"
    
    echo "==================================================" > \$LOG
    echo "HLA Read Extraction - ${sample_id}" >> \$LOG
    echo "==================================================" >> \$LOG
    echo "Started: \$(date)" >> \$LOG
    echo "Input R1: ${r1}" >> \$LOG
    echo "Input R2: ${r2}" >> \$LOG
    echo "" >> \$LOG
    
    # Count input reads
    INPUT_READS=\$(zcat ${r1} | wc -l | awk '{print \$1/4}')
    echo "Input reads: \$INPUT_READS" >> \$LOG
    
    # Method 1: Use seqtk to extract reads mapping to HLA keywords
    # This is faster than BAM conversion for RNA-seq
    echo "Extracting HLA-related reads by sequence similarity..." >> \$LOG
    
    # Create a simple HLA seed pattern file
    cat > hla_seeds.txt << 'EOF'
GCTCCCACTCCATGAGGTATTTC
GTCGCAGCCATAGGGTCTCAG
CACTCCATGAGGTATTTCTAC
GGAGGCGGTCATGACCAGTAC
EOF
    
    # Use grep to find reads containing HLA-like sequences
    # More permissive for RNA-seq data
    zcat ${r1} | paste - - - - | grep -f hla_seeds.txt | tr "\\t" "\\n" | gzip > temp_R1.fastq.gz || true
    zcat ${r2} | paste - - - - | grep -f hla_seeds.txt | tr "\\t" "\\n" | gzip > temp_R2.fastq.gz || true
    
    # If too few reads, use full file (for RNA-seq with HLA transcripts)
    EXTRACTED_READS=\$(zcat temp_R1.fastq.gz 2>/dev/null | wc -l | awk '{print \$1/4}' || echo "0")
    
    if [ \$EXTRACTED_READS -lt 100 ]; then
        echo "WARNING: Extracted only \$EXTRACTED_READS reads. Using full input." >> \$LOG
        cp ${r1} ${sample_id}_HLA_R1.fastq.gz
        cp ${r2} ${sample_id}_HLA_R2.fastq.gz
    else
        echo "Extracted \$EXTRACTED_READS HLA-enriched reads" >> \$LOG
        mv temp_R1.fastq.gz ${sample_id}_HLA_R1.fastq.gz
        mv temp_R2.fastq.gz ${sample_id}_HLA_R2.fastq.gz
    fi
    
    # Calculate extraction efficiency
    OUTPUT_READS=\$(zcat ${sample_id}_HLA_R1.fastq.gz | wc -l | awk '{print \$1/4}')
    PERCENT=\$(echo "scale=2; \$OUTPUT_READS / \$INPUT_READS * 100" | bc)
    
    echo "Output reads: \$OUTPUT_READS" >> \$LOG
    echo "Extraction efficiency: \${PERCENT}%" >> \$LOG
    echo "==================================================" >> \$LOG
    echo "Completed: \$(date)" >> \$LOG
    echo "==================================================" >> \$LOG
    
    # Cleanup
    rm -f hla_seeds.txt temp_*.fastq.gz
    """
}


process FASTQ_TO_BAM_HLA_REGION {
    tag "$sample_id"
    label 'process_medium'
    
    container "${params.container_cache}/samtools_latest.sif"
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_chr6_HLA.bam"), emit: hla_bam
    
    script:
    """
    #!/bin/bash
    set -eo pipefail
    
    # Index reference if needed
    if [ ! -f ${reference_fasta}.fai ]; then
        samtools faidx ${reference_fasta}
    fi
    
    # Extract chr6 HLA region from reference
    samtools faidx ${reference_fasta} chr6:28477797-33448354 > hla_region.fa || \\
        samtools faidx ${reference_fasta} 6:28477797-33448354 > hla_region.fa
    
    # Index HLA region
    bwa index hla_region.fa
    
    # Align reads to HLA region only
    bwa mem -t ${task.cpus} hla_region.fa ${r1} ${r2} | \\
        samtools view -b -F 4 | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}_chr6_HLA.bam
    
    # Index BAM
    samtools index ${sample_id}_chr6_HLA.bam
    
    # Report alignment statistics
    samtools flagstat ${sample_id}_chr6_HLA.bam
    """
}
