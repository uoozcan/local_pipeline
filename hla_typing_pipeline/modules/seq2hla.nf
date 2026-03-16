/*
 * seq2HLA Module
 * RNA-seq HLA typing, Class I and II, 4-field resolution with per-allele confidence
 * Complements arcasHLA (Kallisto-based) with a read-mapping approach
 * FASTQ input only — best suited for RNA-seq data (params.seq_type == 'rna')
 * Reference is bundled inside the container; no external database required
 */

process SEQ2HLA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/seq2hla", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_seq2hla.txt"), emit: results
    // seq2HLA v2.3 actual output filenames (note: -ClassI-class, not -ClassI):
    path("${sample_id}.-ClassI-class.HLAgenotype4digits"),    emit: class1,    optional: true
    path("${sample_id}.-ClassII.HLAgenotype4digits"),          emit: class2,    optional: true
    path("${sample_id}.-ClassI-nonclass.HLAgenotype4digits"),  emit: nonclass1, optional: true
    path "versions.yml", emit: versions

    script:
    """
    echo "[seq2HLA] Running on ${sample_id}..."

    # seq2HLA writes output files prefixed by the -r value.
    # With -r ${sample_id}. the tool produces:
    #   ${sample_id}.-ClassI-class.HLAgenotype4digits    (A, B, C)
    #   ${sample_id}.-ClassI-nonclass.HLAgenotype4digits (E, F, G, ... — not used)
    #   ${sample_id}.-ClassII.HLAgenotype4digits          (DRB1, DQA1, DQB1, DPA1, DPB1)
    seq2HLA \
        -r ${sample_id}. \
        -p ${task.cpus} \
        -1 ${fastq1} \
        -2 ${fastq2} || true

    # Parse classical Class I + Class II result files → standard pipeline TSV
    # Uses the Python 2/3 compatible parse_seq2hla_results.py script
    # (seq2HLA container has Python 2.7 only — script uses 'python' shebang)
    python ${projectDir}/bin/parse_seq2hla_results.py \
        --sample ${sample_id} \
        --prefix "${sample_id}." \
        --output ${sample_id}_seq2hla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seq2hla: "2.3"
    END_VERSIONS
    """
}
