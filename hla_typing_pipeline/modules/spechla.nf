/*
 * SpecHLA Module
 * High-resolution HLA typing from WGS/WES/RNA-seq data
 * Supports both BAM and FASTQ inputs
 *
 * Container mode (default): Uses spechla_with_spechap.sif with compiled SpecHap
 * Local mode: Set params.use_local_spechla = true and params.spechla_path
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
    def use_local = params.use_local_spechla ?: false
    """
    # Set up SpecHLA environment
    SPECHLA_PATH="${params.spechla_path}"

    if [ "${use_local}" = "true" ]; then
        # Local installation mode - may need custom library paths
        export PATH="\${SPECHLA_PATH}/spechla_env/bin:\${SPECHLA_PATH}/bin:\${PATH}"
        # Add local library path if set (for htslib compatibility)
        if [ -n "${params.local_lib_path ?: ''}" ]; then
            export LD_LIBRARY_PATH="${params.local_lib_path}:\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
        else
            export LD_LIBRARY_PATH="\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
        fi
    else
        # Container mode - libraries are properly installed in /usr/local
        export PATH="\${SPECHLA_PATH}/spechla_env/bin:\${SPECHLA_PATH}/bin:/usr/local/bin:\${PATH}"
        export LD_LIBRARY_PATH="/usr/local/lib:\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
    fi

    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index ${bam}
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
    samtools view -b ${bam} \$HLA_REGION > ${sample_id}/hla_extract.bam
    samtools index ${sample_id}/hla_extract.bam

    # Step 2: Convert to FASTQ
    echo "[Step 2] Converting to FASTQ..."
    samtools sort -n ${sample_id}/hla_extract.bam -o ${sample_id}/namesort.bam
    # samtools 1.3.1 doesn't auto-compress, output to uncompressed then gzip
    samtools fastq \
        -1 ${sample_id}/R1.fastq \
        -2 ${sample_id}/R2.fastq \
        -0 /dev/null -s /dev/null \
        ${sample_id}/namesort.bam
    gzip ${sample_id}/R1.fastq
    gzip ${sample_id}/R2.fastq

    # Step 3: Run SpecHLA
    echo "[Step 3] Running SpecHLA..."
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
    def use_local = params.use_local_spechla ?: false
    """
    # Set up SpecHLA environment
    SPECHLA_PATH="${params.spechla_path}"

    if [ "${use_local}" = "true" ]; then
        # Local installation mode - may need custom library paths
        export PATH="\${SPECHLA_PATH}/spechla_env/bin:\${SPECHLA_PATH}/bin:\${PATH}"
        # Add local library path if set (for htslib compatibility)
        if [ -n "${params.local_lib_path ?: ''}" ]; then
            export LD_LIBRARY_PATH="${params.local_lib_path}:\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
        else
            export LD_LIBRARY_PATH="\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
        fi
    else
        # Container mode - libraries are properly installed in /usr/local
        export PATH="\${SPECHLA_PATH}/spechla_env/bin:\${SPECHLA_PATH}/bin:/usr/local/bin:\${PATH}"
        export LD_LIBRARY_PATH="/usr/local/lib:\${SPECHLA_PATH}/spechla_env/lib:\${LD_LIBRARY_PATH:-}"
    fi

    # Create output directory
    mkdir -p ${sample_id}

    # Link FASTQ files using absolute paths (relative symlinks break after 'cd ${sample_id}')
    if [[ "${fastq1}" == *.gz ]]; then
        ln -s "\$(realpath ${fastq1})" ${sample_id}/R1.fastq.gz
        ln -s "\$(realpath ${fastq2})" ${sample_id}/R2.fastq.gz
    else
        gzip -c ${fastq1} > ${sample_id}/R1.fastq.gz
        gzip -c ${fastq2} > ${sample_id}/R2.fastq.gz
    fi

    # Run SpecHLA
    echo "[Running SpecHLA from FASTQ...]"
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
