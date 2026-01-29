/*
 * Consensus Module
 * Weighted voting across HLA typing tools based on read confidence
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
    tuple val(sample_id), path("${sample_id}_consensus.txt"), emit: consensus
    tuple val(sample_id), path("${sample_id}_comparison.txt"), emit: comparison
    path "versions.yml", emit: versions

    script:
    def tools_str = tools.join(',')
    def files_str = result_files.join(',')
    def expected_reads = params.expected_reads ?: 1000
    def weighting = params.weighting ?: 'read_confidence'
    """
    # Run weighted consensus voting
    consensus_voting.py \
        --sample ${sample_id} \
        --tools ${tools_str} \
        --files ${files_str} \
        --resolution ${resolution} \
        --min-tools ${min_tools} \
        --expected-reads ${expected_reads} \
        --weighting ${weighting} \
        --output-consensus ${sample_id}_consensus.txt \
        --output-comparison ${sample_id}_comparison.txt

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
        consensus_voting: "2.0.0"
    END_VERSIONS
    """
}
