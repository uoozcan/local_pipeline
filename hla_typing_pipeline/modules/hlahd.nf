/*
 * HLA-HD Module
 * Accurate HLA typing from WGS/WES data
 * Supports both BAM and FASTQ inputs
 */

process HLAHD {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlahd", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference
    val hla_genes

    output:
    tuple val(sample_id), path("${sample_id}_hlahd.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/${sample_id}/result/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def ref_version = reference == 'hg19' ? '37' : '38'
    """
    # Create working directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Extract HLA reads from BAM
    echo "[Step 1] Extracting HLA reads..."
    samtools view -b -h ${bam} chr6:28000000-34000000 > hla_region.bam 2>/dev/null || \
    samtools view -b -h ${bam} 6:28000000-34000000 > hla_region.bam

    # Add unmapped reads
    samtools view -b -f 4 ${bam} > unmapped.bam
    samtools merge -f merged.bam hla_region.bam unmapped.bam
    samtools sort -n -@ ${task.cpus} merged.bam -o sorted.bam

    # Convert to FASTQ
    samtools fastq -@ ${task.cpus} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam

    # Run HLA-HD
    echo "[Step 2] Running HLA-HD..."
    hlahd.sh -t ${task.cpus} -m 100 -c 0.95 -f ${params.hlahd_db}/freq_data \
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results with per-gene read counts
    # Use cut to extract exactly 3 fields (avoids trailing fields for HLA-E etc.)
    echo "[Step 3] Parsing results with read counts...]"

    RESULT_DIR="${sample_id}/${sample_id}/result"
    if [ -f "\${RESULT_DIR}/${sample_id}_final.result.txt" ]; then
        while IFS='' read -r LINE; do
            GENE=\$(printf '%s\n' "\$LINE" | cut -f1)
            ALLELE1=\$(printf '%s\n' "\$LINE" | cut -f2)
            ALLELE2=\$(printf '%s\n' "\$LINE" | cut -f3)
            GENE_SHORT="\${GENE#HLA-}"
            READ_FILE="\${RESULT_DIR}/${sample_id}_\${GENE_SHORT}.read.txt"
            READS=0
            if [ -f "\${READ_FILE}" ]; then
                RAW=\$(awk 'NR==1{print \$2}' "\${READ_FILE}" 2>/dev/null)
                [[ "\${RAW}" =~ ^[0-9]+\$ ]] && READS=\${RAW}
            fi
            printf '%s\t%s\t%s\t%s\t%s\n' "\${GENE}" "\${ALLELE1}" "\${ALLELE2}" "\${READS}" "\${READS}"
        done < "\${RESULT_DIR}/${sample_id}_final.result.txt" > ${sample_id}_hlahd.txt
    else
        echo "# HLA-HD results for ${sample_id}" > ${sample_id}_hlahd.txt
        echo "# No results generated" >> ${sample_id}_hlahd.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: \$(hlahd.sh 2>&1 | grep -i version | head -1 || echo "1.4.0")
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * HLA-HD from paired FASTQ files
 */
process HLAHD_FASTQ {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlahd", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val hla_genes

    output:
    tuple val(sample_id), path("${sample_id}_hlahd.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/${sample_id}/result/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create working directory
    mkdir -p ${sample_id}

    # Prepare FASTQ files (decompress if needed)
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > R1.fastq
        zcat ${fastq2} > R2.fastq
    else
        ln -s ${fastq1} R1.fastq
        ln -s ${fastq2} R2.fastq
    fi

    # Run HLA-HD
    echo "[Running HLA-HD from FASTQ...]"
    hlahd.sh -t ${task.cpus} -m 100 -c 0.95 -f ${params.hlahd_db}/freq_data \
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split.txt ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results with per-gene read counts
    echo "[Parsing results with read counts...]"
    RESULT_DIR="${sample_id}/${sample_id}/result"
    if [ -f "\${RESULT_DIR}/${sample_id}_final.result.txt" ]; then
        while IFS='' read -r LINE; do
            GENE=\$(printf '%s\n' "\$LINE" | cut -f1)
            ALLELE1=\$(printf '%s\n' "\$LINE" | cut -f2)
            ALLELE2=\$(printf '%s\n' "\$LINE" | cut -f3)
            GENE_SHORT="\${GENE#HLA-}"
            READ_FILE="\${RESULT_DIR}/${sample_id}_\${GENE_SHORT}.read.txt"
            READS=0
            if [ -f "\${READ_FILE}" ]; then
                RAW=\$(awk 'NR==1{print \$2}' "\${READ_FILE}" 2>/dev/null)
                [[ "\${RAW}" =~ ^[0-9]+\$ ]] && READS=\${RAW}
            fi
            printf '%s\t%s\t%s\t%s\t%s\n' "\${GENE}" "\${ALLELE1}" "\${ALLELE2}" "\${READS}" "\${READS}"
        done < "\${RESULT_DIR}/${sample_id}_final.result.txt" > ${sample_id}_hlahd.txt
    else
        echo "# HLA-HD results for ${sample_id}" > ${sample_id}_hlahd.txt
        echo "# No results generated" >> ${sample_id}_hlahd.txt
    fi

    # Cleanup
    rm -f R1.fastq R2.fastq

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: \$(hlahd.sh 2>&1 | grep -i version | head -1 || echo "1.4.0")
    END_VERSIONS
    """
}
