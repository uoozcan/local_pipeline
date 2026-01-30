/*
 * HLA Loss of Heterozygosity (LOH) Module
 * Detect HLA allele loss using SpecHLA's copy number analysis
 *
 * Requirements:
 * - SpecHLA typing results (with frequency files)
 * - Tumor purity estimate (from ABSOLUTE, ASCAT, etc.)
 * - Tumor ploidy estimate
 */

process HLA_LOH {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/loh", mode: 'copy'

    input:
    tuple val(sample_id), path(spechla_dir)
    val purity
    val ploidy

    output:
    tuple val(sample_id), path("${sample_id}_hla_loh.txt"), emit: loh_results
    tuple val(sample_id), path("${sample_id}_loh_summary.txt"), emit: loh_summary
    path "versions.yml", emit: versions

    when:
    purity != null && ploidy != null

    script:
    def het_cutoff = params.loh_het_cutoff ?: 5
    """
    # Create frequency file list
    ls ${spechla_dir}/*_freq.txt > freq.list 2>/dev/null || touch freq.list

    # Check if frequency files exist
    if [ ! -s freq.list ]; then
        echo "WARNING: No frequency files found for LOH analysis" >&2
        echo "Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH" > ${sample_id}_hla_loh.txt
        echo "No frequency files available for LOH analysis" > ${sample_id}_loh_summary.txt
    else
        # Find the typing result file
        TYPING_FILE=\$(find ${spechla_dir} -name "hla.result.txt" -o -name "*_spechla.txt" | head -1)

        if [ -z "\$TYPING_FILE" ]; then
            echo "WARNING: No typing result file found" >&2
            echo "Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH" > ${sample_id}_hla_loh.txt
            echo "No typing result file available" > ${sample_id}_loh_summary.txt
        else
            # Run HLA LOH calculation
            hla_loh.py \\
                --sample ${sample_id} \\
                --purity ${purity} \\
                --ploidy ${ploidy} \\
                --freq-list freq.list \\
                --typing-file \$TYPING_FILE \\
                --het-cutoff ${het_cutoff} \\
                --output ${sample_id}_hla_loh.txt \\
                --summary ${sample_id}_loh_summary.txt
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        hla_loh: "1.0.0"
    END_VERSIONS
    """
}

process HLA_LOH_VISUALIZE {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/loh", mode: 'copy'

    input:
    tuple val(sample_id), path(loh_results)

    output:
    tuple val(sample_id), path("${sample_id}_loh_plot.png"), emit: plot, optional: true
    tuple val(sample_id), path("${sample_id}_loh_report.html"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    hla_loh_visualize.py \\
        --sample ${sample_id} \\
        --loh-results ${loh_results} \\
        --output-plot ${sample_id}_loh_plot.png \\
        --output-report ${sample_id}_loh_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

process HLA_LOH_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(loh_files)

    output:
    path "hla_loh_summary.tsv", emit: summary
    path "hla_loh_summary.html", emit: report
    path "*.png", emit: plots, optional: true
    path "versions.yml", emit: versions

    script:
    """
    hla_loh_summary.py \\
        --loh-files ${loh_files} \\
        --output-tsv hla_loh_summary.tsv \\
        --output-html hla_loh_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
