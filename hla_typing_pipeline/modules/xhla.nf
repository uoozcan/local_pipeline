/*
 * xHLA Module
 * HLA typing using Human Longevity Inc.'s xHLA algorithm
 * Supports BAM input only (requires aligned reads)
 * Uses k-mer based approach for HLA typing
 */

process XHLA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/xhla", mode: 'copy'

    container "${params.container_dir}/xhla.sif"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_xhla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Determine chromosome naming convention (chr6 vs 6)
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")

    # xHLA expects chr6 coordinates, may need to adjust for non-chr BAMs
    echo "Chromosome prefix detected: '\${CHR_PREFIX}'"

    # Run xHLA
    echo "[Running xHLA...]"
    python /opt/bin/run.py \
        --sample_id ${sample_id} \
        --input_bam_path ${bam} \
        --output_path ${sample_id}/ \
        --full

    # Parse results to standard format
    echo "[Parsing xHLA results...]"
    parse_xhla_results.py \
        --sample ${sample_id} \
        --input ${sample_id}/report-${sample_id}-hla.json \
        --output ${sample_id}_xhla.txt

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xhla: "1.0"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * xHLA from FASTQ files
 * Note: xHLA requires BAM input, so we need to align FASTQ first
 */
process XHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/xhla", mode: 'copy'

    container "${params.container_dir}/xhla.sif"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_xhla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # xHLA requires aligned BAM, so we need to align to HLA reference first
    echo "[Aligning FASTQ to HLA reference...]"

    # Use the xHLA internal reference for alignment
    # Extract reference from container data
    HLA_REF="/opt/data/chr6/hla-chr6.fa"

    # Check if bwa index exists, if not use alternative alignment
    if [ -f "\${HLA_REF}.bwt" ]; then
        # Align with BWA
        bwa mem -t ${task.cpus} \${HLA_REF} ${fastq1} ${fastq2} | \
            samtools sort -@ ${task.cpus} -o ${sample_id}/aligned.bam -
        samtools index ${sample_id}/aligned.bam
    else
        # Use diamond for protein-level alignment (xHLA's internal method)
        # Create a minimal BAM from FASTQ for xHLA processing
        echo "Creating pseudo-aligned BAM for xHLA..."

        # Combine FASTQs and create unaligned BAM
        samtools import -@ ${task.cpus} \
            -1 ${fastq1} -2 ${fastq2} \
            -o ${sample_id}/unaligned.bam

        # For xHLA, we need chromosome 6 aligned reads
        # Use minimap2 if available, otherwise skip alignment
        if command -v minimap2 &> /dev/null; then
            minimap2 -ax sr -t ${task.cpus} \${HLA_REF} ${fastq1} ${fastq2} | \
                samtools sort -@ ${task.cpus} -o ${sample_id}/aligned.bam -
            samtools index ${sample_id}/aligned.bam
        else
            echo "WARNING: Cannot align FASTQ for xHLA without proper index"
            echo "# xHLA results for ${sample_id}" > ${sample_id}_xhla.txt
            echo "# ERROR: FASTQ input not supported without alignment" >> ${sample_id}_xhla.txt
            echo "Gene\tAllele1\tAllele2" >> ${sample_id}_xhla.txt
            exit 0
        fi
    fi

    # Run xHLA on aligned BAM
    echo "[Running xHLA...]"
    python /opt/bin/run.py \
        --sample_id ${sample_id} \
        --input_bam_path ${sample_id}/aligned.bam \
        --output_path ${sample_id}/ \
        --full \
        --delete

    # Parse results to standard format
    echo "[Parsing xHLA results...]"
    if [ -f "${sample_id}/report-${sample_id}-hla.json" ]; then
        parse_xhla_results.py \
            --sample ${sample_id} \
            --input ${sample_id}/report-${sample_id}-hla.json \
            --output ${sample_id}_xhla.txt
    else
        echo "# xHLA results for ${sample_id}" > ${sample_id}_xhla.txt
        echo "# No results generated" >> ${sample_id}_xhla.txt
        echo "Gene\tAllele1\tAllele2" >> ${sample_id}_xhla.txt
    fi

    # Cleanup alignment files
    rm -f ${sample_id}/aligned.bam* ${sample_id}/unaligned.bam

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xhla: "1.0"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
