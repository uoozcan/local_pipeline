/*
 * arcasHLA Module
 * HLA typing from RNA-seq data
 * Supports both BAM and FASTQ inputs
 */

process ARCASHLA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/arcashla", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_arcashla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Extract HLA reads
    echo "[Step 1] Extracting HLA reads with arcasHLA..."
    arcasHLA extract ${bam} \
        -o ${sample_id} \
        -t ${task.cpus} \
        -v

    # Genotype
    echo "[Step 2] Running HLA genotyping..."
    arcasHLA genotype ${sample_id}/*.extracted.fq.gz \
        -o ${sample_id} \
        -t ${task.cpus} \
        -v

    # Parse results
    echo "[Step 3] Parsing results..."
    if [ -f "${sample_id}/${sample_id}.genotype.json" ]; then
        python3 <<EOF
import json
import sys

with open('${sample_id}/${sample_id}.genotype.json', 'r') as f:
    data = json.load(f)

with open('${sample_id}_arcashla.txt', 'w') as out:
    out.write("# arcasHLA results for ${sample_id}\\n")
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene, alleles in sorted(data.items()):
        if isinstance(alleles, list):
            a1 = alleles[0] if len(alleles) > 0 else '-'
            a2 = alleles[1] if len(alleles) > 1 else '-'
            out.write(f"{gene}\\t{a1}\\t{a2}\\n")
EOF
    else
        echo "# arcasHLA results for ${sample_id}" > ${sample_id}_arcashla.txt
        echo "# No results generated" >> ${sample_id}_arcashla.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: \$(arcasHLA --version 2>&1 | head -1 || echo "0.5.0")
    END_VERSIONS
    """
}

/*
 * arcasHLA from paired FASTQ files
 */
process ARCASHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/arcashla", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_arcashla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Prepare FASTQ files
    if [[ "${fastq1}" == *.gz ]]; then
        ln -s ${fastq1} ${sample_id}/${sample_id}.extracted.1.fq.gz
        ln -s ${fastq2} ${sample_id}/${sample_id}.extracted.2.fq.gz
    else
        gzip -c ${fastq1} > ${sample_id}/${sample_id}.extracted.1.fq.gz
        gzip -c ${fastq2} > ${sample_id}/${sample_id}.extracted.2.fq.gz
    fi

    # Run arcasHLA genotype directly on FASTQ
    echo "[Running arcasHLA genotype from FASTQ...]"
    arcasHLA genotype \
        ${sample_id}/${sample_id}.extracted.1.fq.gz \
        ${sample_id}/${sample_id}.extracted.2.fq.gz \
        -o ${sample_id} \
        -t ${task.cpus} \
        -v

    # Parse results
    echo "[Parsing results...]"
    # Find the genotype JSON file
    GENOTYPE_FILE=\$(find ${sample_id} -name "*.genotype.json" | head -1)
    if [ -n "\$GENOTYPE_FILE" ] && [ -f "\$GENOTYPE_FILE" ]; then
        python3 <<EOF
import json
import sys

with open('\$GENOTYPE_FILE', 'r') as f:
    data = json.load(f)

with open('${sample_id}_arcashla.txt', 'w') as out:
    out.write("# arcasHLA results for ${sample_id}\\n")
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene, alleles in sorted(data.items()):
        if isinstance(alleles, list):
            a1 = alleles[0] if len(alleles) > 0 else '-'
            a2 = alleles[1] if len(alleles) > 1 else '-'
            out.write(f"{gene}\\t{a1}\\t{a2}\\n")
EOF
    else
        echo "# arcasHLA results for ${sample_id}" > ${sample_id}_arcashla.txt
        echo "# No results generated" >> ${sample_id}_arcashla.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: \$(arcasHLA --version 2>&1 | head -1 || echo "0.5.0")
    END_VERSIONS
    """
}
