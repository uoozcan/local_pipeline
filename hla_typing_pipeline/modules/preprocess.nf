/*
 * Preprocessing Module
 * Prepare BAM files for HLA typing (indexing, validation)
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

    # Get read count estimate from first 10000 reads
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
