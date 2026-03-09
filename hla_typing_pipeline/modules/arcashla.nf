/*
 * arcasHLA Module
 * HLA typing from RNA-seq or WGS HLA-region data
 * Supports both BAM and FASTQ inputs
 */

process ARCASHLA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/arcashla", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)
    val reference

    output:
    tuple val(sample_id), path("${sample_id}_arcashla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}_arcashla.json"), emit: json_results, optional: true
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def min_count = (params.seq_type == 'dna') ? 5 : 75
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

    # Genotype — name FASTQs as sample_id so JSON output matches
    echo "[Step 2] Running HLA genotyping..."
    # Rename extracted FASTQs to use sample_id as stem
    if ls ${sample_id}/*.extracted.fq.gz 1>/dev/null 2>&1; then
        FQFILES=(\$(ls ${sample_id}/*.extracted.fq.gz))
        arcasHLA genotype "\${FQFILES[@]}" \
            -o ${sample_id} \
            -t ${task.cpus} \
            --min_count ${min_count} \
            -v
    else
        echo "No extracted FASTQ files found" >&2
    fi

    # Copy genotype JSON with consistent name for downstream parsing
    echo "[Step 3] Parsing results..."
    GENOTYPE_JSON=\$(find ${sample_id} -name "*.genotype.json" | head -1)
    if [ -n "\$GENOTYPE_JSON" ] && [ -f "\$GENOTYPE_JSON" ]; then
        cp "\$GENOTYPE_JSON" ${sample_id}_arcashla.json
        python3 - "\$GENOTYPE_JSON" "${sample_id}_arcashla.txt" "${sample_id}" << 'PYEOF'
import json, sys
CLASSICAL = {'A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5',
             'DQA1', 'DQB1', 'DPA1', 'DPB1'}
json_path, out_path, sample = sys.argv[1], sys.argv[2], sys.argv[3]
with open(json_path) as f:
    data = json.load(f)
with open(out_path, 'w') as out:
    out.write("# arcasHLA results for {}\\n".format(sample))
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene, alleles in sorted(data.items()):
        if gene not in CLASSICAL:
            continue
        if isinstance(alleles, list):
            a1 = alleles[0] if len(alleles) > 0 else '-'
            a2 = alleles[1] if len(alleles) > 1 else '-'
            out.write("HLA-{}\\t{}\\t{}\\n".format(gene, a1, a2))
PYEOF
    else
        echo "# arcasHLA results for ${sample_id}" > ${sample_id}_arcashla.txt
        echo "# No results generated" >> ${sample_id}_arcashla.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: \$(arcasHLA 2>&1 | grep -oP '(?<=arcasHLA )[0-9.]+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}

/*
 * arcasHLA from paired FASTQ files (pre-extracted HLA-region reads)
 */
process ARCASHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/arcashla", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_arcashla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}_arcashla.json"), emit: json_results, optional: true
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    // Lower min_count for WGS: genomic reads don't pseudoalign to intron-spanning cDNA
    // kmers, so effective per-gene read count is lower than RNA-seq. Use 5 for DNA.
    def min_count = (params.seq_type == 'dna') ? 5 : 75
    """
    mkdir -p ${sample_id}

    # Link FASTQs named as ${sample_id}.1/2.fq.gz so arcasHLA output JSON
    # is named ${sample_id}.genotype.json (parseable without glob)
    if [[ "${fastq1}" == *.gz ]]; then
        ln -sf "\$(realpath ${fastq1})" ${sample_id}/${sample_id}.1.fq.gz
        ln -sf "\$(realpath ${fastq2})" ${sample_id}/${sample_id}.2.fq.gz
    else
        gzip -c ${fastq1} > ${sample_id}/${sample_id}.1.fq.gz
        gzip -c ${fastq2} > ${sample_id}/${sample_id}.2.fq.gz
    fi

    echo "[Running arcasHLA genotype from FASTQ...]"
    arcasHLA genotype \
        ${sample_id}/${sample_id}.1.fq.gz \
        ${sample_id}/${sample_id}.2.fq.gz \
        -o ${sample_id} \
        -t ${task.cpus} \
        --min_count ${min_count} \
        -v

    # Copy genotype JSON with consistent name for downstream parsing
    echo "[Parsing results...]"
    GENOTYPE_JSON="${sample_id}/${sample_id}.genotype.json"
    if [ ! -f "\$GENOTYPE_JSON" ]; then
        # Fall back to glob in case arcasHLA used a different stem
        GENOTYPE_JSON=\$(find ${sample_id} -name "*.genotype.json" | head -1)
    fi

    if [ -n "\$GENOTYPE_JSON" ] && [ -f "\$GENOTYPE_JSON" ]; then
        cp "\$GENOTYPE_JSON" ${sample_id}_arcashla.json
        python3 - "\$GENOTYPE_JSON" "${sample_id}_arcashla.txt" "${sample_id}" << 'PYEOF'
import json, sys
# Classical HLA genes only — non-classical (DMA, DMB, DOA, DOB, E, F, G, H, J, K, L)
# are not in the calibration ground truth and not needed for consensus voting
CLASSICAL = {'A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5',
             'DQA1', 'DQB1', 'DPA1', 'DPB1'}
json_path, out_path, sample = sys.argv[1], sys.argv[2], sys.argv[3]
with open(json_path) as f:
    data = json.load(f)
with open(out_path, 'w') as out:
    out.write("# arcasHLA results for {}\\n".format(sample))
    out.write("Gene\\tAllele1\\tAllele2\\n")
    for gene, alleles in sorted(data.items()):
        if gene not in CLASSICAL:
            continue
        if isinstance(alleles, list):
            a1 = alleles[0] if len(alleles) > 0 else '-'
            a2 = alleles[1] if len(alleles) > 1 else '-'
            out.write("HLA-{}\\t{}\\t{}\\n".format(gene, a1, a2))
PYEOF
    else
        echo "# arcasHLA results for ${sample_id}" > ${sample_id}_arcashla.txt
        echo "# No results generated" >> ${sample_id}_arcashla.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: \$(arcasHLA 2>&1 | grep -oP '(?<=arcasHLA )[0-9.]+' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
