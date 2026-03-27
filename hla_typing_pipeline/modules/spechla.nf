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
    tuple val(sample_id), path(hla_bam)
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
    if [ ! -f "${hla_bam}.bai" ] && [ ! -f "${hla_bam.baseName}.bai" ]; then
        echo "Creating HLA BAM index..."
        samtools index ${hla_bam}
    fi

    # Step 1: Convert pre-extracted HLA BAM to FASTQ for the local SpecHLA install
    echo "[Step 1] Preparing FASTQs from pipeline-generated HLA BAM..."
    samtools sort -n ${hla_bam} -o ${sample_id}/namesort.bam
    # samtools 1.3.1 doesn't auto-compress, output to uncompressed then gzip
    samtools fastq \
        -1 ${sample_id}/R1.fastq \
        -2 ${sample_id}/R2.fastq \
        -0 /dev/null -s /dev/null \
        ${sample_id}/namesort.bam
    gzip ${sample_id}/R1.fastq
    gzip ${sample_id}/R2.fastq

    # Step 2: Run SpecHLA
    echo "[Step 2] Running SpecHLA from pipeline-generated HLA BAM..."
    cd ${sample_id}
    bash \${SPECHLA_PATH}/script/whole/SpecHLA.sh \
        -n ${sample_id} \
        -1 R1.fastq.gz \
        -2 R2.fastq.gz \
        -o . \
        -j ${task.cpus} \
        -u ${ref == 'hg19' ? 0 : 1}
    cd ..

    # Step 3: Parse results
    echo "[Step 3] Parsing results..."
    if [ -f "${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    elif [ -f "${sample_id}/${sample_id}/hla.result.txt" ]; then
        cp ${sample_id}/${sample_id}/hla.result.txt ${sample_id}_spechla.txt
    else
        echo "# SpecHLA results for ${sample_id}" > ${sample_id}_spechla.txt
        echo "# No results generated" >> ${sample_id}_spechla.txt
    fi

    # Cleanup intermediate files
    rm -f ${sample_id}/namesort.bam

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
