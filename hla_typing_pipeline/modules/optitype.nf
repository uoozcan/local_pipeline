/*
 * OptiType Module
 * HLA Class I typing from WGS/WES/RNA-seq data
 * BAM-mode execution is routed through shared preprocessing and FASTQ-mode execution.
 */

process OPTITYPE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_optitype.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    echo "OPTITYPE BAM mode is deprecated in favor of shared preprocessing + OPTITYPE_FASTQ." >&2
    exit 1
    """
}

/*
 * OptiType from paired FASTQ files
 */
process OPTITYPE_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)
    val seq_type  // 'dna' or 'rna'

    output:
    tuple val(sample_id), path("${sample_id}_optitype.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def type_flag = seq_type == 'rna' ? '--rna' : '--dna'
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Prepare FASTQ files
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} > R1.fastq
        zcat ${fastq2} > R2.fastq
    else
        ln -s ${fastq1} R1.fastq
        ln -s ${fastq2} R2.fastq
    fi

    # Run OptiType
    echo "[Running OptiType from FASTQ...]"
    OptiTypePipeline.py \
        -i R1.fastq R2.fastq \
        ${type_flag} \
        -v \
        -o ${sample_id} \
        -p ${sample_id}

    # Parse results
    echo "[Parsing results...]"
    result_tsv=\$(find ${sample_id} -name "*_result.tsv" | head -1)
    if [ -n "\$result_tsv" ] && [ -f "\$result_tsv" ]; then
        python3 <<EOF
import csv
import sys

with open('\${result_tsv}', 'r') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    row = next(reader)

with open('${sample_id}_optitype.txt', 'w') as out:
    out.write("# OptiType results for ${sample_id}\\n")
    out.write("Gene\\tAllele1\\tAllele2\\n")
    out.write("HLA-A\\t{}\\t{}\\n".format(row.get('A1', '-'), row.get('A2', '-')))
    out.write("HLA-B\\t{}\\t{}\\n".format(row.get('B1', '-'), row.get('B2', '-')))
    out.write("HLA-C\\t{}\\t{}\\n".format(row.get('C1', '-'), row.get('C2', '-')))
EOF
    else
        echo "# OptiType results for ${sample_id}" > ${sample_id}_optitype.txt
        echo "# No results generated" >> ${sample_id}_optitype.txt
    fi

    # Cleanup
    rm -f R1.fastq R2.fastq

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(OptiTypePipeline.py --version 2>&1 | head -1 || echo "1.3.5")
    END_VERSIONS
    """
}
