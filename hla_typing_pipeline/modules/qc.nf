/*
 * QC Module
 * Quality control and validation for HLA typing inputs
 * Provides warnings for insufficient reads or poor quality data
 */

process QC_BAM {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference
    val min_hla_reads
    val min_read_length

    output:
    tuple val(sample_id), path(bam), emit: validated_bam
    tuple val(sample_id), path("${sample_id}_qc_report.txt"), emit: qc_report
    tuple val(sample_id), env(QC_PASS), emit: qc_status

    script:
    def ref = reference == 'hg19' ? 'hg19' : 'hg38'
    """
    # Initialize QC status
    QC_PASS="true"
    WARNINGS=""

    # Create QC report
    echo "# HLA Typing QC Report for ${sample_id}" > ${sample_id}_qc_report.txt
    echo "# Generated: \$(date)" >> ${sample_id}_qc_report.txt
    echo "# Reference: ${ref}" >> ${sample_id}_qc_report.txt
    echo "" >> ${sample_id}_qc_report.txt

    # Check for BAM index
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index for QC..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Determine chromosome naming convention
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -oP "SN:\\K[^\\t]*" | grep -o "^chr" || echo "")

    # Define HLA region based on reference
    if [ "${ref}" == "hg38" ]; then
        if [ -n "\$CHR_PREFIX" ]; then
            HLA_REGION="chr6:28510120-33480577"
        else
            HLA_REGION="6:28510120-33480577"
        fi
    else
        if [ -n "\$CHR_PREFIX" ]; then
            HLA_REGION="chr6:28477797-33448354"
        else
            HLA_REGION="6:28477797-33448354"
        fi
    fi

    echo "=== Input Statistics ===" >> ${sample_id}_qc_report.txt

    # Total reads in BAM
    TOTAL_READS=\$(samtools view -c ${bam})
    echo "Total reads: \$TOTAL_READS" >> ${sample_id}_qc_report.txt

    # HLA region reads
    HLA_READS=\$(samtools view -c ${bam} \$HLA_REGION 2>/dev/null || echo "0")
    echo "HLA region reads: \$HLA_READS" >> ${sample_id}_qc_report.txt
    echo "HLA region: \$HLA_REGION" >> ${sample_id}_qc_report.txt

    # Calculate HLA read percentage
    if [ "\$TOTAL_READS" -gt 0 ]; then
        HLA_PCT=\$(echo "scale=4; \$HLA_READS * 100 / \$TOTAL_READS" | bc)
        echo "HLA read percentage: \${HLA_PCT}%" >> ${sample_id}_qc_report.txt
    fi

    # Average read length (sample first 10000 reads)
    AVG_LENGTH=\$(samtools view ${bam} | head -10000 | awk '{sum+=length(\$10); count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')
    echo "Average read length: \$AVG_LENGTH bp" >> ${sample_id}_qc_report.txt

    # Mapping quality statistics
    AVG_MAPQ=\$(samtools view ${bam} | head -10000 | awk '{sum+=\$5; count++} END {if(count>0) printf "%.1f", sum/count; else print "0"}')
    echo "Average mapping quality: \$AVG_MAPQ" >> ${sample_id}_qc_report.txt

    # Paired-end check
    PAIRED_READS=\$(samtools view -c -f 1 ${bam})
    if [ "\$PAIRED_READS" -gt 0 ]; then
        echo "Read type: Paired-end" >> ${sample_id}_qc_report.txt

        # Proper pairs
        PROPER_PAIRS=\$(samtools view -c -f 2 ${bam})
        PROPER_PCT=\$(echo "scale=2; \$PROPER_PAIRS * 100 / \$PAIRED_READS" | bc)
        echo "Properly paired: \${PROPER_PCT}%" >> ${sample_id}_qc_report.txt
    else
        echo "Read type: Single-end" >> ${sample_id}_qc_report.txt
    fi

    echo "" >> ${sample_id}_qc_report.txt
    echo "=== QC Warnings ===" >> ${sample_id}_qc_report.txt

    # Check minimum HLA reads
    if [ "\$HLA_READS" -lt ${min_hla_reads} ]; then
        echo "WARNING: Low HLA read count (\$HLA_READS < ${min_hla_reads})" >> ${sample_id}_qc_report.txt
        echo "  - HLA typing results may be unreliable" >> ${sample_id}_qc_report.txt
        echo "  - Consider deeper sequencing or enrichment" >> ${sample_id}_qc_report.txt
        WARNINGS="\${WARNINGS}LOW_HLA_READS;"
        QC_PASS="warning"
    fi

    # Check read length
    if [ "\$AVG_LENGTH" -lt ${min_read_length} ]; then
        echo "WARNING: Short average read length (\$AVG_LENGTH < ${min_read_length} bp)" >> ${sample_id}_qc_report.txt
        echo "  - Short reads may reduce typing accuracy" >> ${sample_id}_qc_report.txt
        echo "  - Some tools work better with longer reads" >> ${sample_id}_qc_report.txt
        WARNINGS="\${WARNINGS}SHORT_READS;"
        QC_PASS="warning"
    fi

    # Check mapping quality
    if [ "\$(echo "\$AVG_MAPQ < 20" | bc)" -eq 1 ]; then
        echo "WARNING: Low average mapping quality (\$AVG_MAPQ < 20)" >> ${sample_id}_qc_report.txt
        echo "  - Poor mapping may affect HLA typing accuracy" >> ${sample_id}_qc_report.txt
        WARNINGS="\${WARNINGS}LOW_MAPQ;"
        QC_PASS="warning"
    fi

    # Check if HLA reads are very low (critical)
    if [ "\$HLA_READS" -lt 100 ]; then
        echo "CRITICAL: Very few HLA reads (\$HLA_READS < 100)" >> ${sample_id}_qc_report.txt
        echo "  - HLA typing is unlikely to produce reliable results" >> ${sample_id}_qc_report.txt
        echo "  - Check if sample is from expected organism/tissue" >> ${sample_id}_qc_report.txt
        QC_PASS="critical"
    fi

    if [ -z "\$WARNINGS" ]; then
        echo "No warnings - input data appears suitable for HLA typing" >> ${sample_id}_qc_report.txt
    fi

    echo "" >> ${sample_id}_qc_report.txt
    echo "=== QC Status ===" >> ${sample_id}_qc_report.txt
    echo "Status: \$QC_PASS" >> ${sample_id}_qc_report.txt

    # Print summary to stdout for Nextflow log
    echo ""
    echo "=============================================="
    echo "QC Summary for ${sample_id}"
    echo "=============================================="
    echo "HLA reads: \$HLA_READS (minimum: ${min_hla_reads})"
    echo "Avg read length: \$AVG_LENGTH bp (minimum: ${min_read_length})"
    echo "Avg mapping quality: \$AVG_MAPQ"
    echo "Status: \$QC_PASS"
    if [ "\$QC_PASS" != "true" ]; then
        echo ""
        echo "*** WARNINGS DETECTED - Review ${sample_id}_qc_report.txt ***"
    fi
    echo "=============================================="
    """
}

process QC_FASTQ {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val min_hla_reads
    val min_read_length

    output:
    tuple val(sample_id), path(fastq1), path(fastq2), emit: validated_fastq
    tuple val(sample_id), path("${sample_id}_qc_report.txt"), emit: qc_report
    tuple val(sample_id), env(QC_PASS), emit: qc_status

    script:
    """
    # Initialize QC status
    QC_PASS="true"
    WARNINGS=""

    # Create QC report
    echo "# HLA Typing QC Report for ${sample_id}" > ${sample_id}_qc_report.txt
    echo "# Generated: \$(date)" >> ${sample_id}_qc_report.txt
    echo "# Input type: Paired FASTQ" >> ${sample_id}_qc_report.txt
    echo "" >> ${sample_id}_qc_report.txt

    echo "=== Input Statistics ===" >> ${sample_id}_qc_report.txt
    echo "FASTQ R1: ${fastq1}" >> ${sample_id}_qc_report.txt
    echo "FASTQ R2: ${fastq2}" >> ${sample_id}_qc_report.txt

    # Count reads in FASTQ files
    if [[ "${fastq1}" == *.gz ]]; then
        R1_READS=\$(zcat ${fastq1} | awk 'NR%4==1' | wc -l)
        R2_READS=\$(zcat ${fastq2} | awk 'NR%4==1' | wc -l)

        # Sample read lengths from first 10000 reads
        AVG_LENGTH_R1=\$(zcat ${fastq1} | awk 'NR%4==2' | head -10000 | awk '{sum+=length(\$0); count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')
        AVG_LENGTH_R2=\$(zcat ${fastq2} | awk 'NR%4==2' | head -10000 | awk '{sum+=length(\$0); count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')

        # Sample quality scores
        AVG_QUAL_R1=\$(zcat ${fastq1} | awk 'NR%4==0' | head -10000 | awk '{for(i=1;i<=length(\$0);i++) sum+=ord(substr(\$0,i,1))-33; count+=length(\$0)} END {if(count>0) printf "%.1f", sum/count; else print "0"}' 2>/dev/null || echo "N/A")
    else
        R1_READS=\$(awk 'NR%4==1' ${fastq1} | wc -l)
        R2_READS=\$(awk 'NR%4==1' ${fastq2} | wc -l)

        AVG_LENGTH_R1=\$(awk 'NR%4==2' ${fastq1} | head -10000 | awk '{sum+=length(\$0); count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')
        AVG_LENGTH_R2=\$(awk 'NR%4==2' ${fastq2} | head -10000 | awk '{sum+=length(\$0); count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')
    fi

    TOTAL_READS=\$((R1_READS + R2_READS))
    AVG_LENGTH=\$(echo "scale=0; (\$AVG_LENGTH_R1 + \$AVG_LENGTH_R2) / 2" | bc)

    echo "R1 reads: \$R1_READS" >> ${sample_id}_qc_report.txt
    echo "R2 reads: \$R2_READS" >> ${sample_id}_qc_report.txt
    echo "Total reads: \$TOTAL_READS" >> ${sample_id}_qc_report.txt
    echo "Average read length R1: \$AVG_LENGTH_R1 bp" >> ${sample_id}_qc_report.txt
    echo "Average read length R2: \$AVG_LENGTH_R2 bp" >> ${sample_id}_qc_report.txt
    echo "Average read length: \$AVG_LENGTH bp" >> ${sample_id}_qc_report.txt

    # Check read count balance
    if [ \$R1_READS -ne \$R2_READS ]; then
        DIFF=\$(echo "scale=2; (\$R1_READS - \$R2_READS) * 100 / \$R1_READS" | bc | tr -d '-')
        echo "R1/R2 difference: \${DIFF}%" >> ${sample_id}_qc_report.txt
    fi

    echo "" >> ${sample_id}_qc_report.txt
    echo "=== QC Warnings ===" >> ${sample_id}_qc_report.txt

    # Note: For pre-extracted FASTQ, we assume reads are HLA-relevant
    # The actual HLA read count depends on whether these are whole-genome or pre-extracted

    # Check total read count (assuming these might be HLA-extracted reads)
    if [ "\$TOTAL_READS" -lt ${min_hla_reads} ]; then
        echo "WARNING: Low total read count (\$TOTAL_READS < ${min_hla_reads})" >> ${sample_id}_qc_report.txt
        echo "  - If these are HLA-extracted reads, typing may be unreliable" >> ${sample_id}_qc_report.txt
        echo "  - If whole-genome FASTQ, this is critically low" >> ${sample_id}_qc_report.txt
        WARNINGS="\${WARNINGS}LOW_READS;"
        QC_PASS="warning"
    fi

    # Check read length
    if [ "\$AVG_LENGTH" -lt ${min_read_length} ]; then
        echo "WARNING: Short average read length (\$AVG_LENGTH < ${min_read_length} bp)" >> ${sample_id}_qc_report.txt
        echo "  - Short reads may reduce typing accuracy" >> ${sample_id}_qc_report.txt
        WARNINGS="\${WARNINGS}SHORT_READS;"
        QC_PASS="warning"
    fi

    # Check R1/R2 balance
    if [ \$R1_READS -ne \$R2_READS ]; then
        DIFF_ABS=\$((R1_READS - R2_READS))
        if [ \${DIFF_ABS#-} -gt \$((R1_READS / 10)) ]; then
            echo "WARNING: Unbalanced R1/R2 read counts" >> ${sample_id}_qc_report.txt
            echo "  - R1: \$R1_READS, R2: \$R2_READS" >> ${sample_id}_qc_report.txt
            echo "  - This may indicate truncated or corrupted files" >> ${sample_id}_qc_report.txt
            WARNINGS="\${WARNINGS}UNBALANCED_PAIRS;"
            QC_PASS="warning"
        fi
    fi

    # Check for very low reads (critical)
    if [ "\$TOTAL_READS" -lt 100 ]; then
        echo "CRITICAL: Very few reads (\$TOTAL_READS < 100)" >> ${sample_id}_qc_report.txt
        echo "  - HLA typing is unlikely to produce reliable results" >> ${sample_id}_qc_report.txt
        QC_PASS="critical"
    fi

    if [ -z "\$WARNINGS" ]; then
        echo "No warnings - input data appears suitable for HLA typing" >> ${sample_id}_qc_report.txt
    fi

    echo "" >> ${sample_id}_qc_report.txt
    echo "=== QC Status ===" >> ${sample_id}_qc_report.txt
    echo "Status: \$QC_PASS" >> ${sample_id}_qc_report.txt

    # Print summary to stdout for Nextflow log
    echo ""
    echo "=============================================="
    echo "QC Summary for ${sample_id}"
    echo "=============================================="
    echo "Total reads: \$TOTAL_READS (minimum recommended: ${min_hla_reads})"
    echo "Avg read length: \$AVG_LENGTH bp (minimum: ${min_read_length})"
    echo "Status: \$QC_PASS"
    if [ "\$QC_PASS" != "true" ]; then
        echo ""
        echo "*** WARNINGS DETECTED - Review ${sample_id}_qc_report.txt ***"
    fi
    echo "=============================================="
    """
}
