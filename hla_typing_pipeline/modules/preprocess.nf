/*
 * Preprocessing Module
 * Extract HLA reads and convert to FASTQ for downstream HLA typing
 * Also handles CRAM → BAM conversion (requires reference FASTA)
 * Uses basetools container with samtools
 */

/*
 * Validate BAM file
 */
process CHECK_BAM {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), emit: validated_bam
    path "bam_info.txt", emit: info

    script:
    """
    # Validate BAM file
    echo "Sample: ${sample_id}" > bam_info.txt
    echo "BAM: ${bam}" >> bam_info.txt

    # Check if BAM is sorted
    if samtools view -H ${bam} | grep -q "SO:coordinate"; then
        echo "Sorted: yes" >> bam_info.txt
    else
        echo "Sorted: unknown (may not be coordinate sorted)" >> bam_info.txt
    fi

    # Get read count estimate
    READ_COUNT=\$(samtools view -c ${bam} | head -1 || echo "unknown")
    echo "Reads: \${READ_COUNT}" >> bam_info.txt

    # Check for index
    if [ -f "${bam}.bai" ] || [ -f "${bam.baseName}.bai" ]; then
        echo "Index: found" >> bam_info.txt
    else
        echo "Index: not found (will be created by HLA tools)" >> bam_info.txt
    fi

    cat bam_info.txt
    """
}

/*
 * Extract HLA reads from BAM and convert to FASTQ
 * This process runs with the basetools container to ensure samtools is available
 */
process EXTRACT_HLA_READS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/preprocessing", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_hla_R1.fastq.gz"), path("${sample_id}_hla_R2.fastq.gz"), emit: fastq
    tuple val(sample_id), path("${sample_id}_extraction.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def ref = reference == 'hg19' ? 'hg19' : 'hg38'
    """
    # Create log file
    echo "HLA Read Extraction Log for ${sample_id}" > ${sample_id}_extraction.log
    echo "Reference: ${ref}" >> ${sample_id}_extraction.log
    echo "Date: \$(date)" >> ${sample_id}_extraction.log
    echo "" >> ${sample_id}_extraction.log

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..." | tee -a ${sample_id}_extraction.log
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Determine chromosome naming convention
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    echo "Chromosome prefix: '\${CHR_PREFIX:-none}'" >> ${sample_id}_extraction.log

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
    echo "HLA region: \$HLA_REGION" >> ${sample_id}_extraction.log

    # Step 1: Extract HLA reads
    echo "Extracting HLA reads from \$HLA_REGION..." | tee -a ${sample_id}_extraction.log
    samtools view -@ ${task.cpus} -b ${bam} \$HLA_REGION > hla_extract.bam

    # Count extracted reads
    HLA_READS=\$(samtools view -c hla_extract.bam)
    echo "Extracted HLA reads: \$HLA_READS" >> ${sample_id}_extraction.log

    # Step 2: Name sort
    echo "Name sorting..." | tee -a ${sample_id}_extraction.log
    samtools sort -n -@ ${task.cpus} -m 2G hla_extract.bam -o namesort.bam

    # Step 3: Convert to FASTQ
    echo "Converting to FASTQ..." | tee -a ${sample_id}_extraction.log
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}_hla_R1.fastq.gz \
        -2 ${sample_id}_hla_R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        -c 6 \
        namesort.bam

    # Count output reads
    R1_READS=\$(zcat ${sample_id}_hla_R1.fastq.gz | awk 'NR%4==1' | wc -l)
    echo "Output R1 reads: \$R1_READS" >> ${sample_id}_extraction.log
    echo "Output R2 reads: \$R1_READS" >> ${sample_id}_extraction.log

    # Cleanup
    rm -f hla_extract.bam namesort.bam

    echo "Extraction complete." >> ${sample_id}_extraction.log

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * Extract HLA reads for HLA-HD (includes unmapped reads)
 */
process EXTRACT_HLA_READS_HLAHD {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_hlahd_R1.fastq"), path("${sample_id}_hlahd_R2.fastq"), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    # Check for BAM index
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Extract HLA region reads (expanded region for HLA-HD)
    samtools view -b -h ${bam} chr6:28000000-34000000 > hla_region.bam 2>/dev/null || \
    samtools view -b -h ${bam} 6:28000000-34000000 > hla_region.bam

    # Extract unmapped reads
    samtools view -b -f 4 ${bam} > unmapped.bam

    # Merge and sort
    samtools merge -f merged.bam hla_region.bam unmapped.bam
    samtools sort -n -@ ${task.cpus} merged.bam -o sorted.bam

    # Convert to FASTQ (uncompressed for HLA-HD)
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}_hlahd_R1.fastq \
        -2 ${sample_id}_hlahd_R2.fastq \
        -0 /dev/null -s /dev/null \
        sorted.bam

    # Cleanup
    rm -f hla_region.bam unmapped.bam merged.bam sorted.bam

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * Convert CRAM to coordinate-sorted BAM
 *
 * CRAM files are compressed against a reference genome and require the same
 * reference FASTA for decoding. The output BAM is indexed and ready for all
 * downstream HLA typing tools.
 *
 * Usage in pipeline (set reference_fasta in nextflow.config or run.config):
 *   CRAM_TO_BAM(ch_cram, params.reference_fasta)
 *
 * Samplesheet: use 'cram_path' column instead of 'bam_path'
 * CLI:         --input_cram sample.cram --reference_fasta /path/to/ref.fa
 */
process CRAM_TO_BAM {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/cram_to_bam", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(cram)
    path(reference_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path("${sample_id}_cram_to_bam.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    echo "CRAM to BAM conversion for ${sample_id}" > ${sample_id}_cram_to_bam.log
    echo "CRAM: ${cram}" >> ${sample_id}_cram_to_bam.log
    echo "Reference: ${reference_fasta}" >> ${sample_id}_cram_to_bam.log
    echo "Date: \$(date)" >> ${sample_id}_cram_to_bam.log

    # Create CRAM index if missing (.crai)
    if [ ! -f "${cram}.crai" ] && [ ! -f "${cram.baseName}.crai" ]; then
        echo "Creating CRAM index..." | tee -a ${sample_id}_cram_to_bam.log
        samtools index ${cram}
    fi

    # Decode CRAM → coordinate-sorted BAM using the reference
    samtools view -@ ${task.cpus} -b -T ${reference_fasta} -o ${sample_id}_unsorted.bam ${cram} \
        2>>${sample_id}_cram_to_bam.log \
        || { echo "ERROR: samtools view (CRAM decode) failed" | tee -a ${sample_id}_cram_to_bam.log; exit 1; }

    # Sort by coordinate (CRAM may not be coordinate-sorted or index may differ)
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam ${sample_id}_unsorted.bam \
        2>>${sample_id}_cram_to_bam.log

    # Index the output BAM
    samtools index ${sample_id}.bam

    # Summary
    READS=\$(samtools view -c ${sample_id}.bam)
    echo "Output BAM reads: \${READS}" | tee -a ${sample_id}_cram_to_bam.log

    rm -f ${sample_id}_unsorted.bam

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
	END_VERSIONS
    """
}

/*
 * Simple BAM to FASTQ conversion (full BAM, no region extraction)
 */
process BAM_TO_FASTQ {
    tag "$sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    # Sort by name and convert to FASTQ
    samtools sort -n -@ ${task.cpus} ${bam} -o namesort.bam
    samtools fastq -@ ${task.cpus} \
        -1 ${sample_id}_R1.fastq.gz \
        -2 ${sample_id}_R2.fastq.gz \
        -0 /dev/null -s /dev/null \
        -c 6 \
        namesort.bam

    rm -f namesort.bam

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
