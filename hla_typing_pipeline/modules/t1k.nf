/*
 * T1K Module
 * HLA typing from WGS/RNA-seq FASTQ data using the T1K assembler-based approach.
 * Covers Class I + II: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1
 * Reference: https://github.com/mourisl/T1K
 *
 * Container: quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0
 * Runtime:   ~5 min WGS 30x HLA-enriched reads, ~4 GB RAM
 */

/*
 * T1K Long-Read Module
 * HLA typing from PacBio HiFi or Oxford Nanopore single-file reads.
 * Uses T1K's platform-specific presets for long-read alignment.
 *
 * Platform: 'hifi'  → --preset HiFi  (PacBio CCS/HiFi)
 *           'ont'   → --preset ONT   (Oxford Nanopore)
 *
 * Input: single-file BAM or FASTQ (not paired-end)
 * Output: standard pipeline TSV identical to T1K_FASTQ
 */

process T1K_LONGREADS {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/t1k", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val platform  // 'hifi' or 'ont'

    output:
    tuple val(sample_id), path("${sample_id}_t1k.txt"),  emit: results
    tuple val(sample_id), path("t1k_out/*"),              emit: full_results, optional: true
    path "versions.yml",                                  emit: versions

    script:
    def preset = platform == 'hifi' ? 'HiFi' : 'ONT'
    """
    mkdir -p t1k_out

    # Locate bundled HLA reference inside the biocontainers image
    T1K_SHARE=\$(find /usr/local/share -maxdepth 1 -name "t1k-*" -type d 2>/dev/null | head -1 || true)
    if [ -z "\${T1K_SHARE}" ]; then
        T1K_SHARE="/usr/local/share/t1k"
    fi
    HLA_FA="\${T1K_SHARE}/hlaidx/hla_dna.fa"
    GENE_LIST="\${T1K_SHARE}/hlaidx/h_gene_list.txt"

    # Run T1K with single-file input (-u) and platform preset
    if [ -f "\${HLA_FA}" ] && [ -f "\${GENE_LIST}" ]; then
        run-t1k \\
            -u ${reads} \\
            -f "\${HLA_FA}" \\
            --genotype-list "\${GENE_LIST}" \\
            --preset ${preset} \\
            -t ${task.cpus} \\
            -o t1k_out/${sample_id} \\
            2>t1k_out/t1k_stderr.log
    else
        # Fallback: --preset hla uses T1K's built-in HLA reference
        run-t1k \\
            -u ${reads} \\
            --preset hla \\
            -t ${task.cpus} \\
            -o t1k_out/${sample_id} \\
            2>t1k_out/t1k_stderr.log
    fi

    # Parse genotype TSV → standard output format
    GENOTYPE_FILE=\$(ls t1k_out/${sample_id}_genotype.tsv t1k_out/${sample_id}.genotype.tsv 2>/dev/null | head -1 || true)
    if [ -n "\${GENOTYPE_FILE}" ] && [ -f "\${GENOTYPE_FILE}" ]; then
        python3 ${projectDir}/bin/parse_t1k_results.py \\
            "\${GENOTYPE_FILE}" \\
            "${sample_id}_t1k.txt"
    else
        echo "# T1K long-read results for ${sample_id} (platform: ${preset})" > ${sample_id}_t1k.txt
        echo "# No results generated" >> ${sample_id}_t1k.txt
        cat t1k_out/t1k_stderr.log >> ${sample_id}_t1k.txt || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "1.0.9")
        platform: "${preset}"
    END_VERSIONS
    """
}

// ─────────────────────────────────────────────────────────────────────────────

process T1K_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/t1k", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_t1k.txt"),  emit: results
    tuple val(sample_id), path("t1k_out/*"),              emit: full_results, optional: true
    path "versions.yml",                                  emit: versions

    script:
    """
    mkdir -p t1k_out

    # Locate bundled HLA reference inside the biocontainers image.
    # T1K ships its IMGT/HLA index at /usr/local/share/t1k-*/hlaidx/
    T1K_SHARE=\$(find /usr/local/share -maxdepth 1 -name "t1k-*" -type d 2>/dev/null | head -1)
    if [ -z "\${T1K_SHARE}" ]; then
        T1K_SHARE="/usr/local/share/t1k"
    fi
    HLA_FA="\${T1K_SHARE}/hlaidx/hla_dna.fa"
    GENE_LIST="\${T1K_SHARE}/hlaidx/h_gene_list.txt"

    # Decompress FASTQs if needed (T1K can handle .gz natively but some versions prefer plain)
    if [[ "${fastq1}" == *.gz ]]; then
        F1="${fastq1}"
        F2="${fastq2}"
    else
        F1="${fastq1}"
        F2="${fastq2}"
    fi

    # Run T1K — use bundled reference if found, otherwise --preset hla
    if [ -f "\${HLA_FA}" ] && [ -f "\${GENE_LIST}" ]; then
        run-t1k \\
            -1 "\${F1}" \\
            -2 "\${F2}" \\
            -f "\${HLA_FA}" \\
            --genotype-list "\${GENE_LIST}" \\
            -t ${task.cpus} \\
            -o t1k_out/${sample_id} \\
            2>t1k_out/t1k_stderr.log
    else
        run-t1k \\
            --preset hla \\
            -1 "\${F1}" \\
            -2 "\${F2}" \\
            -t ${task.cpus} \\
            -o t1k_out/${sample_id} \\
            2>t1k_out/t1k_stderr.log
    fi

    # Parse genotype TSV → standard output format
    GENOTYPE_FILE=\$(ls t1k_out/${sample_id}_genotype.tsv t1k_out/${sample_id}.genotype.tsv 2>/dev/null | head -1 || true)
    if [ -n "\${GENOTYPE_FILE}" ] && [ -f "\${GENOTYPE_FILE}" ]; then
        python3 ${projectDir}/bin/parse_t1k_results.py \\
            "\${GENOTYPE_FILE}" \\
            "${sample_id}_t1k.txt"
    else
        echo "# T1K results for ${sample_id}" > ${sample_id}_t1k.txt
        echo "# No results generated" >> ${sample_id}_t1k.txt
        cat t1k_out/t1k_stderr.log >> ${sample_id}_t1k.txt || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "1.0.9")
    END_VERSIONS
    """
}
