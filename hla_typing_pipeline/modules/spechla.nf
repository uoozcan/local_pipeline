/*
 * SpecHLA Module
 * High-resolution HLA typing from WGS/WES/RNA-seq data
 * Supports both BAM and FASTQ inputs
 * Uses Singularity container (default) - spechla_1.0.7-3.sif
 * Container includes all dependencies (SpecHap, databases, etc.)
 */

process SPECHLA {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/spechla", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_spechla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def ref = reference == 'hg19' ? 'hg19' : 'hg38'
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Determine chromosome naming convention
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")

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

    # Step 1: Extract HLA reads
    echo "[Step 1] Extracting HLA reads from \$HLA_REGION..."
    samtools view -@ ${task.cpus} -b ${bam} \$HLA_REGION > ${sample_id}/hla_extract.bam
    samtools index ${sample_id}/hla_extract.bam

    # Step 2: Convert to FASTQ
    echo "[Step 2] Converting to FASTQ..."
    samtools sort -n -@ ${task.cpus} ${sample_id}/hla_extract.bam -o ${sample_id}/namesort.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}/R1.fastq.gz \
        -2 ${sample_id}/R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        ${sample_id}/namesort.bam

    # Step 3: Run SpecHLA
    echo "[Step 3] Running SpecHLA..."
    SPECHLA_PATH="${params.spechla_path}"
    cd ${sample_id}
    bash \${SPECHLA_PATH}/script/whole/SpecHLA.sh \
        -n ${sample_id} \
        -1 R1.fastq.gz \
        -2 R2.fastq.gz \
        -o . \
        -j ${task.cpus} \
        -u 1
    cd ..

    # Step 4: Parse results
    echo "[Step 4] Parsing results..."
    if [ -f "${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    elif [ -f "${sample_id}/${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    else
        echo "# SpecHLA results for ${sample_id}" > ${sample_id}_spechla.txt
        echo "# No results generated" >> ${sample_id}_spechla.txt
    fi

    # Cleanup intermediate files
    rm -f ${sample_id}/hla_extract.bam* ${sample_id}/namesort.bam

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spechla: "1.0.7"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * SpecHLA from paired FASTQ files
 */
process SPECHLA_FASTQ {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/spechla", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_spechla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Link or copy FASTQ files to working directory
    if [[ "${fastq1}" == *.gz ]]; then
        ln -s ${fastq1} ${sample_id}/R1.fastq.gz
        ln -s ${fastq2} ${sample_id}/R2.fastq.gz
    else
        gzip -c ${fastq1} > ${sample_id}/R1.fastq.gz
        gzip -c ${fastq2} > ${sample_id}/R2.fastq.gz
    fi

    # Run SpecHLA
    echo "[Running SpecHLA from FASTQ...]"
    SPECHLA_PATH="${params.spechla_path}"
    cd ${sample_id}
    bash \${SPECHLA_PATH}/script/whole/SpecHLA.sh \
        -n ${sample_id} \
        -1 R1.fastq.gz \
        -2 R2.fastq.gz \
        -o . \
        -j ${task.cpus} \
        -u 1
    cd ..

    # Parse results
    echo "[Parsing results...]"
    if [ -f "${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    elif [ -f "${sample_id}/${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    else
        echo "# SpecHLA results for ${sample_id}" > ${sample_id}_spechla.txt
        echo "# No results generated" >> ${sample_id}_spechla.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spechla: "1.0.7"
    END_VERSIONS
    """
}
