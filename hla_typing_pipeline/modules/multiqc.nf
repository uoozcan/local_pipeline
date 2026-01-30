/*
 * MultiQC Module
 * Aggregate quality control reports
 */

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('fastqc/*')
    path('hla_reports/*')
    path(multiqc_config)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "versions.yml", emit: versions

    script:
    def config_arg = multiqc_config ? "--config ${multiqc_config}" : ""
    """
    multiqc . ${config_arg} \\
        --title "HLA Typing Pipeline Report" \\
        --comment "Multi-tool HLA typing with consensus voting"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
