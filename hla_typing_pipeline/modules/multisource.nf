/*
 * Multi-source Consensus Module
 * Integrates HLA typing results across multiple sequencing data sources
 * (WGS, WES, RNAseq, targeted) for the same patient.
 *
 * Inputs per patient:
 *   - Per-source consensus files (output of CONSENSUS process)
 *   - Parallel lists of sample_ids and seq_types
 *
 * Outputs per patient:
 *   - {patient_id}_integrated_consensus.txt  — weighted cross-source allele calls
 *   - {patient_id}_source_comparison.txt     — per-locus per-source breakdown + discordance flag
 *   - {patient_id}_integrated_report.html    — self-contained HTML summary
 *
 * Version: 1.3.0
 */

process MULTISOURCE_CONSENSUS {
    tag "$patient_id"
    label 'process_low'
    publishDir "${params.outdir}/${patient_id}/integrated", mode: 'copy'

    input:
    tuple val(patient_id), val(sample_ids), val(seq_types), path(consensus_files)
    val source_weights
    val resolution

    output:
    tuple val(patient_id), path("${patient_id}_integrated_consensus.txt"),  emit: integrated_consensus
    tuple val(patient_id), path("${patient_id}_source_comparison.txt"),     emit: source_comparison
    tuple val(patient_id), path("${patient_id}_integrated_report.html"),    emit: html_report
    path "versions.yml",                                                     emit: versions

    script:
    // Convert lists to comma-separated strings safely (handles both single and multi-element input)
    def samples_str  = sample_ids  instanceof List ? sample_ids.join(',')  : sample_ids
    def seqtype_str  = seq_types   instanceof List ? seq_types.join(',')   : seq_types
    def files_str    = consensus_files instanceof List
                         ? consensus_files.collect { it.toString() }.join(',')
                         : consensus_files.toString()
    """
    multisource_consensus.py \\
        --patient       "${patient_id}" \\
        --samples       "${samples_str}" \\
        --seq-types     "${seqtype_str}" \\
        --files         "${files_str}" \\
        --source-weights '${source_weights}' \\
        --resolution    ${resolution} \\
        --output-consensus  ${patient_id}_integrated_consensus.txt \\
        --output-comparison ${patient_id}_source_comparison.txt \\
        --output-report     ${patient_id}_integrated_report.html

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        multisource_consensus: "1.0.0"
    END_VERSIONS
    """
}
