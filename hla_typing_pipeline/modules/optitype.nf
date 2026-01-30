/*
 * OptiType Module
 * HLA Class I typing from WGS/WES/RNA-seq data
 * Supports both BAM and FASTQ inputs
 */

process OPTITYPE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_optitype.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Extract reads and convert to FASTQ
    echo "[Step 1] Extracting reads..."
    samtools sort -n -@ ${task.cpus} ${bam} -o sorted.bam
    samtools fastq -@ ${task.cpus} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam

    # Run OptiType
    echo "[Step 2] Running OptiType..."
    OptiTypePipeline.py \
        -i R1.fastq R2.fastq \
        --dna \
        -v \
        -o ${sample_id} \
        -p ${sample_id}

    # Parse results
    echo "[Step 3] Parsing results..."
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
    out.write(f"HLA-A\\t{row.get('A1', '-')}\\t{row.get('A2', '-')}\\n")
    out.write(f"HLA-B\\t{row.get('B1', '-')}\\t{row.get('B2', '-')}\\n")
    out.write(f"HLA-C\\t{row.get('C1', '-')}\\t{row.get('C2', '-')}\\n")
EOF
    else
        echo "# OptiType results for ${sample_id}" > ${sample_id}_optitype.txt
        echo "# No results generated" >> ${sample_id}_optitype.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(OptiTypePipeline.py --version 2>&1 | head -1 || echo "1.3.5")
    END_VERSIONS
    """
}

/*
 * OptiType from paired FASTQ files
 */
process OPTITYPE_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'

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
    out.write(f"HLA-A\\t{row.get('A1', '-')}\\t{row.get('A2', '-')}\\n")
    out.write(f"HLA-B\\t{row.get('B1', '-')}\\t{row.get('B2', '-')}\\n")
    out.write(f"HLA-C\\t{row.get('C1', '-')}\\t{row.get('C2', '-')}\\n")
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
