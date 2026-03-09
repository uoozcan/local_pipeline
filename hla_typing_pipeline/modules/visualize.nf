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

process HLA_PIPELINE_METRICS {
    tag "pipeline_metrics"
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(trace_file)

    output:
    path("metrics_*.png"),                    emit: plots,  optional: true
    path("execution_metrics_report.html"),    emit: report, optional: true

    script:
    """
    hla_pipeline_metrics.py \\
        --trace ${trace_file} \\
        --outdir . \\
        --system-ram ${task.memory ? (task.memory.toGiga() as int) : 14}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)" 2>/dev/null || echo "N/A")
        plotly: \$(python3 -c "import plotly; print(plotly.__version__)" 2>/dev/null || echo "N/A")
    END_VERSIONS
    """
}
