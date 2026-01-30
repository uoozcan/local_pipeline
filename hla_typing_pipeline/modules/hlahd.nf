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
    tuple val(sample_id), path("${sample_id}/result/*"), emit: full_results, optional: true
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
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results
    echo "[Step 3] Parsing results..."
    mkdir -p ${sample_id}/result
    if [ -f "${sample_id}/result/${sample_id}_final.result.txt" ]; then
        cp ${sample_id}/result/${sample_id}_final.result.txt ${sample_id}_hlahd.txt
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
    tuple val(sample_id), path("${sample_id}/result/*"), emit: full_results, optional: true
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
        R1.fastq R2.fastq ${params.hlahd_db}/HLA_gene.split ${params.hlahd_db}/dictionary \
        ${sample_id} ${sample_id}

    # Parse results
    echo "[Parsing results...]"
    mkdir -p ${sample_id}/result
    if [ -f "${sample_id}/result/${sample_id}_final.result.txt" ]; then
        cp ${sample_id}/result/${sample_id}_final.result.txt ${sample_id}_hlahd.txt
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
