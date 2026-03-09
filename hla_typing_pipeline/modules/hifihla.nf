/*
 * HiFi-HLA Module
 * HLA typing at up to 4-field resolution from PacBio HiFi (CCS) reads.
 *
 * Container : quay.io/pacbio/hifihla:latest
 * Input     : HiFi BAM (aligned to GRCh38, or unaligned CCS reads)
 * Output    : 4-field allele calls (e.g. A*03:01:01:01)
 *
 * Reference : https://github.com/PacificBiosciences/hifihla
 *
 * Key distinguishing feature: long-read length enables true 4-field
 * resolution across all classical HLA loci without additional phasing steps.
 */

process HIFIHLA {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hifihla", mode: 'copy'
    errorStrategy 'ignore'   // gracefully skip if container unavailable

    input:
    tuple val(sample_id), path(hifi_bam)

    output:
    tuple val(sample_id), path("${sample_id}_hifihla.txt"),  emit: results
    tuple val(sample_id), path("hifihla_out/"),              emit: full_results, optional: true
    path "versions.yml",                                     emit: versions

    script:
    """
    mkdir -p hifihla_out

    # ── Run HiFi-HLA ────────────────────────────────────────────────────────
    # Supports two invocation styles depending on the installed version.
    # v1.x: hifihla star-call-from-reads  (reads may be BAM or FASTQ)
    # v0.x: hifihla call-star-alleles     (older API)
    #
    # The tool auto-detects whether input is aligned or unaligned.
    if hifihla star-call-from-reads --help >/dev/null 2>&1; then
        hifihla star-call-from-reads \\
            --reads ${hifi_bam} \\
            --output-dir hifihla_out/ \\
            2>hifihla_out/hifihla_stderr.log
    else
        hifihla call-star-alleles \\
            --bam ${hifi_bam} \\
            --output-dir hifihla_out/ \\
            2>hifihla_out/hifihla_stderr.log
    fi

    # ── Parse output → standard pipeline TSV ──────────────────────────────
    # HiFi-HLA writes one or more *_calls.tsv files in output-dir.
    # Columns: gene  allele1  confidence1  allele2  confidence2
    CALLS_FILE=\$(ls hifihla_out/*_calls.tsv hifihla_out/*hla_calls.tsv 2>/dev/null | head -1 || true)

    if [ -n "\${CALLS_FILE}" ] && [ -f "\${CALLS_FILE}" ]; then
        python3 ${projectDir}/bin/parse_hifihla_results.py \\
            "\${CALLS_FILE}" \\
            "${sample_id}_hifihla.txt"
    else
        # No calls file produced — write a placeholder so the pipeline
        # continues gracefully (CONSENSUS will simply receive 0 alleles
        # from this tool for this sample).
        echo "# HiFi-HLA results for ${sample_id}" > ${sample_id}_hifihla.txt
        echo "# No calls file produced — check hifihla_out/hifihla_stderr.log" \\
            >> ${sample_id}_hifihla.txt
        cat hifihla_out/hifihla_stderr.log >> ${sample_id}_hifihla.txt || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifihla: \$(hifihla --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.?[0-9]*' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
