/*
 * HLAscan Module
 * HLA typing using SyntekaBio's HLAscan algorithm
 * Supports BAM and FASTQ input
 *
 * IMPORTANT LIMITATIONS:
 *   - HLAscan v2.1 container only supports GRCh37/hg19 (-v 37)
 *   - GRCh38/hg38 support is broken in v2.1 (binary exits with "wrong HLA gene")
 *   - Input BAM must be a full genome BAM (not HLA-region-only subset)
 *   - HLAscan internally extracts reads from chr6 coordinates; HLA-region BAMs
 *     will have all reads mixed across all HLA genes, causing gene-resolution failure
 *
 * For hg38 BAMs, the process outputs NA values and exits cleanly.
 * For hg19 BAMs with full genome coverage, all 8 HLA genes are typed per gene.
 */

process HLASCAN {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/hlascan", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_hlascan.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Detect reference genome version from chr6 length
    # hg19 chr6: 171,115,067 bp; hg38 chr6: 170,805,979 bp
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -o "SN:[^\t]*" | cut -d: -f2 | grep -o "^chr" || echo "")
    CHR6_LEN=\$(samtools view -H ${bam} | grep "^@SQ.*SN:\${CHR_PREFIX}6\b" | grep -o "LN:[0-9]*" | cut -d: -f2)

    if [ "\${CHR6_LEN}" -gt "171000000" ]; then
        REF_VERSION="hg19"
        HLASCAN_V="37"
    else
        REF_VERSION="hg38"
        HLASCAN_V="38"
    fi
    echo "Detected reference: \${REF_VERSION} (chr6 length: \${CHR6_LEN})"

    if [ "\${REF_VERSION}" == "hg19" ]; then
        echo "[Running HLAscan on hg19 BAM (-v 37)...]"

        for GENE in HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQA1 HLA-DQB1 HLA-DPA1 HLA-DPB1; do
            echo "[Typing \${GENE}...]"
            hlascan \
                -b ${bam} \
                -d /db \
                -v \${HLASCAN_V} \
                -g \${GENE} \
                -t ${task.cpus} \
                > ${sample_id}/\${GENE}.txt 2>&1 || true
        done

        # Parse all per-gene outputs into standard format
        parse_hlascan_results.py \
            --sample ${sample_id} \
            --input-dir ${sample_id}/ \
            --output ${sample_id}_hlascan.txt
    else
        echo "[HLAscan: hg38 not supported by v2.1 container — writing NA output]"
        echo "# HLAscan results for ${sample_id}" > ${sample_id}_hlascan.txt
        echo "# HLAscan typing skipped - hg38 not supported by HLAscan v2.1 container" >> ${sample_id}_hlascan.txt
        printf 'Gene\\tAllele1\\tAllele2\\tReads1\\tReads2\\n' >> ${sample_id}_hlascan.txt
        for GENE in A B C DRB1 DQA1 DQB1 DPA1 DPB1; do
            printf '%s\\tNA\\tNA\\tNA\\tNA\\n' "\${GENE}" >> ${sample_id}_hlascan.txt
        done
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlascan: "2.1"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * HLAscan from FASTQ files
 * Note: -l / -r FASTQ flags feed raw reads to HLAscan's internal aligner.
 * HLAscan maps ALL input reads to HLA reference; with HLA-enriched reads (from
 * all genes), gene-resolution may fail ("wrong HLA gene" error). Best results
 * are achieved with reads from an HLA-region BAM using the BAM mode above.
 */
process HLASCAN_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/hlascan", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_hlascan.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${sample_id}

    echo "[Running HLAscan on FASTQ...]"
    echo "[Note: HLAscan FASTQ mode works best with hg19 DB and uncompressed FASTQ]"

    # Decompress if needed (hlascan v2.1 may not support gzipped input)
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > ${sample_id}/R1.fq
        zcat ${fastq2} > ${sample_id}/R2.fq
        FQ1="${sample_id}/R1.fq"
        FQ2="${sample_id}/R2.fq"
    else
        FQ1="${fastq1}"
        FQ2="${fastq2}"
    fi

    HLASCAN_SUCCESS=false
    for GENE in HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQA1 HLA-DQB1 HLA-DPA1 HLA-DPB1; do
        echo "[Typing \${GENE}...]"
        hlascan \
            -l \${FQ1} \
            -r \${FQ2} \
            -d ${params.hlascan_db} \
            -g \${GENE} \
            -t ${task.cpus} \
            > ${sample_id}/\${GENE}.txt 2>&1 || true

        # Check if result looks valid (not "wrong HLA gene")
        if grep -q "\\[Type " "${sample_id}/\${GENE}.txt" 2>/dev/null; then
            HLASCAN_SUCCESS=true
        fi
    done

    # Cleanup uncompressed FASTQ to save space
    rm -f ${sample_id}/R1.fq ${sample_id}/R2.fq

    if [ "\${HLASCAN_SUCCESS}" == "true" ]; then
        parse_hlascan_results.py \
            --sample ${sample_id} \
            --input-dir ${sample_id}/ \
            --output ${sample_id}_hlascan.txt
    else
        echo "# HLAscan results for ${sample_id}" > ${sample_id}_hlascan.txt
        echo "# HLAscan FASTQ typing failed - hg38 or gene-resolution error" >> ${sample_id}_hlascan.txt
        printf 'Gene\\tAllele1\\tAllele2\\tReads1\\tReads2\\n' >> ${sample_id}_hlascan.txt
        for GENE in A B C DRB1 DQA1 DQB1 DPA1 DPB1; do
            printf '%s\\tNA\\tNA\\tNA\\tNA\\n' "\${GENE}" >> ${sample_id}_hlascan.txt
        done
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlascan: "2.1"
    END_VERSIONS
    """
}
