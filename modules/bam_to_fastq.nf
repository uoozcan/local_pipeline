// modules/bam_to_fastq.nf
// Optimized BAM processing with improved error handling

process EXTRACT_HLA_REGION {
    tag "$sample_id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${sample_id}/hla_extraction", mode: 'copy', pattern: "*.log"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.hla.bam"), emit: hla_bam
    path "${sample_id}.extract.log", emit: log
    
    script:
    def hla_start = params.hla_region_start ?: 28000000
    def hla_end = params.hla_region_end ?: 34000000
    """
    echo "==================================================" > ${sample_id}.extract.log
    echo "HLA Region Extraction - ${sample_id}" >> ${sample_id}.extract.log
    echo "==================================================" >> ${sample_id}.extract.log
    echo "Started: \$(date -Iseconds)" >> ${sample_id}.extract.log
    echo "BAM file: ${bam}" >> ${sample_id}.extract.log
    echo "" >> ${sample_id}.extract.log
    
    # Ensure BAM is indexed
    if [ ! -f ${bam}.bai ] && [ ! -f ${bam.baseName}.bai ]; then
        echo "Creating BAM index..." >> ${sample_id}.extract.log
        samtools index -@ ${task.cpus} ${bam}
    fi
    
    # Detect chromosome 6 naming convention (chr6 vs 6)
    CHR=\$(samtools view -H ${bam} | awk '
        \$1=="@SQ" {
            for(i=1; i<=NF; i++) {
                if(\$i ~ /^SN:/) {
                    chr = substr(\$i, 4)
                    if(chr == "6") { print "6"; exit }
                    if(chr == "chr6") { chr6_found = "chr6" }
                }
            }
        }
        END { if(chr6_found) print chr6_found }
    ')
    
    if [ -z "\$CHR" ]; then
        echo "ERROR: Could not detect chromosome 6 naming" >> ${sample_id}.extract.log
        exit 1
    fi
    
    echo "Chromosome 6 detected as: \$CHR" >> ${sample_id}.extract.log
    HLA_REGION="\${CHR}:${hla_start}-${hla_end}"
    echo "HLA region: \$HLA_REGION" >> ${sample_id}.extract.log
    echo "CPUs: ${task.cpus}" >> ${sample_id}.extract.log
    echo "Memory: ${task.memory}" >> ${sample_id}.extract.log
    echo "" >> ${sample_id}.extract.log
    
    # Extract HLA region
    echo "Extracting HLA region..." >> ${sample_id}.extract.log
    samtools view -b -@ ${task.cpus} ${bam} \$HLA_REGION > ${sample_id}.hla.bam
    
    # Index the extracted BAM
    samtools index -@ ${task.cpus} ${sample_id}.hla.bam
    
    # Report statistics
    echo "" >> ${sample_id}.extract.log
    echo "Statistics:" >> ${sample_id}.extract.log
    TOTAL_READS=\$(samtools view -c -@ ${task.cpus} ${bam})
    HLA_READS=\$(samtools view -c -@ ${task.cpus} ${sample_id}.hla.bam)
    
    echo "  Total reads in BAM: \$TOTAL_READS" >> ${sample_id}.extract.log
    echo "  HLA region reads: \$HLA_READS" >> ${sample_id}.extract.log
    
    if [ \$TOTAL_READS -gt 0 ]; then
        PERCENTAGE=\$(echo "scale=2; \$HLA_READS * 100 / \$TOTAL_READS" | bc)
        echo "  Percentage: \${PERCENTAGE}%" >> ${sample_id}.extract.log
    fi
    
    BAM_SIZE=\$(du -h ${sample_id}.hla.bam | cut -f1)
    echo "  HLA BAM size: \$BAM_SIZE" >> ${sample_id}.extract.log
    
    echo "" >> ${sample_id}.extract.log
    echo "Completed: \$(date -Iseconds)" >> ${sample_id}.extract.log
    echo "SUCCESS" >> ${sample_id}.extract.log
    echo "==================================================" >> ${sample_id}.extract.log
    """
}

process BAM_TO_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${sample_id}/fastq", mode: 'copy', enabled: params.save_intermediate
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    path "${sample_id}_bam2fastq.log", emit: log
    
    script:
    """
    echo "==================================================" > ${sample_id}_bam2fastq.log
    echo "BAM to FASTQ Conversion - ${sample_id}" >> ${sample_id}_bam2fastq.log
    echo "==================================================" >> ${sample_id}_bam2fastq.log
    echo "Started: \$(date -Iseconds)" >> ${sample_id}_bam2fastq.log
    echo "BAM file: ${bam}" >> ${sample_id}_bam2fastq.log
    echo "CPUs: ${task.cpus}" >> ${sample_id}_bam2fastq.log
    echo "Memory: ${task.memory}" >> ${sample_id}_bam2fastq.log
    echo "" >> ${sample_id}_bam2fastq.log
    
    # Check for pigz, fall back to gzip
    if command -v pigz &> /dev/null; then
        GZBIN="pigz"
        echo "Using pigz for compression" >> ${sample_id}_bam2fastq.log
    else
        GZBIN="gzip"
        echo "Using gzip for compression" >> ${sample_id}_bam2fastq.log
    fi
    
    echo "Converting BAM to FASTQ..." >> ${sample_id}_bam2fastq.log
    
    # Optimized pipeline: collate -> fastq -> compress in parallel
    set -euo pipefail
    samtools collate -@ ${task.cpus} -O ${bam} ${sample_id}_tmp \\
    | samtools fastq -@ ${task.cpus} -N -F 0x900 \\
        -1 >(\$GZBIN -c > ${sample_id}_R1.fastq.gz) \\
        -2 >(\$GZBIN -c > ${sample_id}_R2.fastq.gz) \\
        -0 /dev/null -s /dev/null -
    
    # Report statistics
    echo "" >> ${sample_id}_bam2fastq.log
    echo "Statistics:" >> ${sample_id}_bam2fastq.log
    READ_PAIRS=\$(zcat ${sample_id}_R1.fastq.gz | wc -l | awk '{print int(\$1/4)}')
    R1_SIZE=\$(du -h ${sample_id}_R1.fastq.gz | cut -f1)
    R2_SIZE=\$(du -h ${sample_id}_R2.fastq.gz | cut -f1)
    
    echo "  Read pairs: \$READ_PAIRS" >> ${sample_id}_bam2fastq.log
    echo "  R1 size: \$R1_SIZE" >> ${sample_id}_bam2fastq.log
    echo "  R2 size: \$R2_SIZE" >> ${sample_id}_bam2fastq.log
    
    echo "" >> ${sample_id}_bam2fastq.log
    echo "Completed: \$(date -Iseconds)" >> ${sample_id}_bam2fastq.log
    echo "SUCCESS" >> ${sample_id}_bam2fastq.log
    echo "==================================================" >> ${sample_id}_bam2fastq.log
    """
}

process EXTRACT_HLA_AND_CONVERT {
    tag "$sample_id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*.log"
    publishDir "${params.outdir}/${sample_id}/fastq", mode: 'copy', pattern: "*.fastq.gz", enabled: params.save_intermediate
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    path "${sample_id}_hla_conversion.log", emit: log
    
    script:
    def hla_start = params.hla_region_start ?: 28000000
    def hla_end = params.hla_region_end ?: 34000000
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "==================================================" > ${sample_id}_hla_conversion.log
    echo "HLA Region Extraction & FASTQ Conversion" >> ${sample_id}_hla_conversion.log
    echo "Sample: ${sample_id}" >> ${sample_id}_hla_conversion.log
    echo "==================================================" >> ${sample_id}_hla_conversion.log
    echo "Started: \$(date -Iseconds)" >> ${sample_id}_hla_conversion.log
    echo "BAM file: ${bam}" >> ${sample_id}_hla_conversion.log
    echo "" >> ${sample_id}_hla_conversion.log
    
    # Ensure BAM is indexed
    if [ ! -f ${bam}.bai ] && [ ! -f ${bam.baseName}.bai ]; then
        echo "Creating BAM index..." >> ${sample_id}_hla_conversion.log
        samtools index -@ ${task.cpus} ${bam} 2>&1 | tee -a ${sample_id}_hla_conversion.log
    fi
    
    # Detect chromosome 6 naming convention
    echo "Detecting chromosome 6 naming..." >> ${sample_id}_hla_conversion.log
    CHR=\$(samtools view -H ${bam} 2>&1 | awk '
        \$1=="@SQ" {
            for(i=1; i<=NF; i++) {
                if(\$i ~ /^SN:/) {
                    chr = substr(\$i, 4)
                    if(chr == "6") { print "6"; exit }
                    if(chr == "chr6") { chr6_found = "chr6" }
                }
            }
        }
        END { if(chr6_found) print chr6_found }
    ')
    
    if [ -z "\$CHR" ]; then
        echo "ERROR: Could not detect chromosome 6 naming" >> ${sample_id}_hla_conversion.log
        echo "Available chromosomes:" >> ${sample_id}_hla_conversion.log
        samtools view -H ${bam} | grep "^@SQ" | head -20 >> ${sample_id}_hla_conversion.log
        exit 1
    fi
    
    echo "Configuration:" >> ${sample_id}_hla_conversion.log
    echo "  Chromosome 6 detected as: \$CHR" >> ${sample_id}_hla_conversion.log
    HLA_REGION="\${CHR}:${hla_start}-${hla_end}"
    echo "  HLA region: \$HLA_REGION" >> ${sample_id}_hla_conversion.log
    echo "  CPUs: ${task.cpus}" >> ${sample_id}_hla_conversion.log
    echo "  Memory: ${task.memory}" >> ${sample_id}_hla_conversion.log
    echo "" >> ${sample_id}_hla_conversion.log
    
    # Check for compression tools
    if command -v pigz &> /dev/null; then
        GZBIN="pigz"
        GZCMD="pigz -c"
        echo "Using pigz (fast compression)" >> ${sample_id}_hla_conversion.log
    else
        GZBIN="gzip"
        GZCMD="gzip -c"
        echo "Using gzip for compression" >> ${sample_id}_hla_conversion.log
    fi
    echo "" >> ${sample_id}_hla_conversion.log
    
    # Optimized single-pass pipeline
    echo "Extracting HLA region and converting to FASTQ..." >> ${sample_id}_hla_conversion.log
    
    # Create named pipes for better process handling
    mkfifo ${sample_id}_R1.pipe ${sample_id}_R2.pipe
    
    # Start compression in background
    \$GZCMD < ${sample_id}_R1.pipe > ${sample_id}_R1.fastq.gz &
    PID_R1=\$!
    \$GZCMD < ${sample_id}_R2.pipe > ${sample_id}_R2.fastq.gz &
    PID_R2=\$!
    
    # Run samtools pipeline
    samtools view -@ ${task.cpus} -b ${bam} "\$HLA_REGION" 2>> ${sample_id}_hla_conversion.log \\
    | samtools collate -@ ${task.cpus} -O - ${sample_id}_tmp 2>> ${sample_id}_hla_conversion.log \\
    | samtools fastq -@ ${task.cpus} -N -F 0x900 \\
        -1 ${sample_id}_R1.pipe \\
        -2 ${sample_id}_R2.pipe \\
        -0 /dev/null -s /dev/null - 2>> ${sample_id}_hla_conversion.log
    
    # Wait for compression to finish
    wait \$PID_R1
    wait \$PID_R2
    
    # Clean up pipes
    rm -f ${sample_id}_R1.pipe ${sample_id}_R2.pipe
    
    # Verify output files exist and are not empty
    if [ ! -s ${sample_id}_R1.fastq.gz ] || [ ! -s ${sample_id}_R2.fastq.gz ]; then
        echo "ERROR: Output FASTQ files are missing or empty" >> ${sample_id}_hla_conversion.log
        exit 1
    fi
    
    # Report statistics
    echo "" >> ${sample_id}_hla_conversion.log
    echo "Statistics:" >> ${sample_id}_hla_conversion.log
    
    TOTAL_READS=\$(samtools view -c -@ ${task.cpus} ${bam} 2>/dev/null || echo "0")
    HLA_READ_PAIRS=\$(zcat ${sample_id}_R1.fastq.gz | wc -l | awk '{print int(\$1/4)}')
    R1_SIZE=\$(du -h ${sample_id}_R1.fastq.gz | cut -f1)
    R2_SIZE=\$(du -h ${sample_id}_R2.fastq.gz | cut -f1)
    
    echo "  Total reads in BAM: \$TOTAL_READS" >> ${sample_id}_hla_conversion.log
    echo "  HLA region read pairs: \$HLA_READ_PAIRS" >> ${sample_id}_hla_conversion.log
    
    if [ \$TOTAL_READS -gt 0 ]; then
        PERCENTAGE=\$(echo "scale=2; \$HLA_READ_PAIRS * 200 / \$TOTAL_READS" | bc 2>/dev/null || echo "N/A")
        echo "  Percentage in HLA region: \${PERCENTAGE}%" >> ${sample_id}_hla_conversion.log
    fi
    
    echo "  R1 file size: \$R1_SIZE" >> ${sample_id}_hla_conversion.log
    echo "  R2 file size: \$R2_SIZE" >> ${sample_id}_hla_conversion.log
    
    echo "" >> ${sample_id}_hla_conversion.log
    echo "Completed: \$(date -Iseconds)" >> ${sample_id}_hla_conversion.log
    echo "SUCCESS" >> ${sample_id}_hla_conversion.log
    echo "==================================================" >> ${sample_id}_hla_conversion.log
    """
}
