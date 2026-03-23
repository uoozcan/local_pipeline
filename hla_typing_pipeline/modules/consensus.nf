/*
 * Consensus Module
 * Weighted voting across HLA typing tools based on read confidence.
 *
 * Supported output formats (params.output_format):
 *   text       - standard TSV consensus + comparison files (default)
 *   gl_string  - text + GL String Consortium file (<sample>_gl_string.txt)
 *   hml        - text + gl_string + HML v1.0.1 XML  (<sample>.hml.xml)
 *   all        - all of the above
 */

process CONSENSUS {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), val(tools), path(result_files)
    val resolution
    val min_tools

    output:
    tuple val(sample_id), path("${sample_id}_consensus.txt"),              emit: consensus
    tuple val(sample_id), path("${sample_id}_comparison.txt"),             emit: comparison
    tuple val(sample_id), path("${sample_id}_gl_string.txt"), optional: true, emit: gl_string
    tuple val(sample_id), path("${sample_id}.hml.xml"),       optional: true, emit: hml
    path "versions.yml",                                                   emit: versions

    script:
    def tools_str      = tools.join(',')
    def files_str      = result_files.join(',')
    def expected_reads = params.expected_reads ?: 1000
    def weighting      = params.weighting ?: 'read_confidence'
    def output_format  = params.output_format ?: 'text'

    // calibrated weighting: auto-select weights file by seq_type if not explicitly set
    def resolvedWeightsFile = params.weights_file
    if (params.weighting == 'calibrated' && !params.weights_file) {
        def dtMap = [dna: 'wgs', wes: 'wes', rna: 'rna',
                     longreads_hifi: 'hifi', longreads_ont: 'ont']
        def dtLabel = dtMap[params.seq_type] ?: 'wgs'
        resolvedWeightsFile = "${projectDir}/conf/tool_weights_${dtLabel}.json"
    }
    def weights_arg = (params.weighting == 'calibrated' && resolvedWeightsFile)
        ? "--weights-file ${resolvedWeightsFile} --data-type ${params.seq_type ?: 'dna'}"
        : ''

    // gl_string output path arg (passed to consensus_voting.py)
    def do_gl  = output_format in ['gl_string', 'hml', 'all']
    def gl_arg = do_gl ? "--output-gl-string ${sample_id}_gl_string.txt" : ''

    // HML generation args
    def do_hml         = output_format in ['hml', 'all']
    def allele_db_ver  = params.allele_db_version ?: '3.57.0'
    def center_id      = params.hml_center_id ?: 'HLA-PIPELINE'

    """
    # ── Consensus voting ──────────────────────────────────────────────────
    consensus_voting.py \\
        --sample ${sample_id} \\
        --tools ${tools_str} \\
        --files ${files_str} \\
        --resolution ${resolution} \\
        --min-tools ${min_tools} \\
        --expected-reads ${expected_reads} \\
        --weighting ${weighting} \\
        ${weights_arg} \\
        --output-format ${output_format} \\
        ${gl_arg} \\
        --output-consensus ${sample_id}_consensus.txt \\
        --output-comparison ${sample_id}_comparison.txt

    # ── HML v1.0.1 XML ────────────────────────────────────────────────────
    if ${do_hml}; then
        generate_hml.py \\
            --sample ${sample_id} \\
            --consensus ${sample_id}_consensus.txt \\
            --output ${sample_id}.hml.xml \\
            --allele-db IMGT/HLA \\
            --allele-db-version ${allele_db_ver} \\
            --center-id ${center_id}
    fi

    # ── Version manifest ──────────────────────────────────────────────────
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        consensus_voting: "2.1.0"
        generate_hml: "1.0.0"
    END_VERSIONS
    """
}
