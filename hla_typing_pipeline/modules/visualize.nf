/*
 * Visualization Module
 * Generate plots and statistics for HLA typing results
 */

process HLA_VISUALIZE {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/visualizations", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_file), path(comparison_file)
    tuple val(sample_id), path(qc_report)

    output:
    tuple val(sample_id), path("*.png"), emit: plots
    tuple val(sample_id), path("*.html"), emit: html_report
    tuple val(sample_id), path("${sample_id}_statistics.json"), emit: statistics
    path "versions.yml", emit: versions

    script:
    """
    hla_visualize.py \\
        --sample ${sample_id} \\
        --consensus ${consensus_file} \\
        --comparison ${comparison_file} \\
        --qc-report ${qc_report} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)" 2>/dev/null || echo "N/A")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """
}

process HLA_SUMMARY_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(consensus_files)
    path(comparison_files)
    path(statistics_files)

    output:
    path "hla_summary_report.html", emit: html_report
    path "hla_summary_statistics.tsv", emit: statistics
    path "*.png", emit: plots
    path "versions.yml", emit: versions

    script:
    """
    hla_summary_report.py \\
        --consensus-files ${consensus_files} \\
        --comparison-files ${comparison_files} \\
        --statistics-files ${statistics_files} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
