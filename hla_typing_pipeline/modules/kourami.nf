/*
 * Kourami Module
 * Assembly-graph HLA typing, Class I (A/B/C) + Class II (DQA1/DRB1), 4-field resolution
 * Adds algorithmic diversity: assembly-graph approach distinct from all current pipeline tools
 * BAM input only — requires hg38/GRCh38-aligned coordinate-sorted BAM
 *
 * Prerequisites (set in nextflow.config / run.config):
 *   params.kourami_dir  — path to Kourami install dir (contains build/Kourami.jar)
 *   params.kourami_db   — path to IMGT-formatted HLA DB (contains All_FINAL_with_Decoy.fa.gz + BWA index)
 *
 * HLA read extraction:
 *   By default: samtools extracts reads from the HLA region directly (no reference needed).
 *   Optional: set params.kourami_hs38_ref to use Kourami's alignAndExtract_hs38DH.sh script instead.
 *
 * Host tools used: samtools, bwa, java ≥11 (or provide via container)
 *   Container must include: bash, bwa, samtools, java ≥11
 *   Build from kourami_preprocess.dockerfile in linnil1/HLA_collections + add openjdk-11-jdk
 */

process KOURAMI {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/kourami", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_kourami.txt"), emit: results
    path("${sample_id}.kourami.result"), emit: raw, optional: true
    path "versions.yml", emit: versions

    script:
    def kourami_dir = params.kourami_dir ?: '/opt/kourami'
    def kourami_db  = params.kourami_db  ?: '/opt/kourami_db'
    def hs38_ref    = params.kourami_hs38_ref  // null = use samtools direct extraction
    """
    echo "[Kourami] Running on ${sample_id}..."

    # Ensure BAM index exists (required for region-based extraction)
    [ -f "${bam}.bai" ] || samtools index ${bam}

    # Step 1: Extract HLA-region reads.
    # If kourami_hs38_ref is set, use Kourami's alignAndExtract script for precise extraction.
    # Otherwise, use samtools to extract reads directly from the HLA region (chr6 28M-34M).
    # Auto-detect chromosome naming: check if BAM uses "chr6" or "6".
    if [ "${hs38_ref}" != "null" ] && [ -n "${hs38_ref}" ]; then
        bash ${kourami_dir}/scripts/alignAndExtract_hs38DH.sh \
            -d ${kourami_db} \
            -r ${hs38_ref} \
            ${sample_id}. ${bam} || true
        R1="${sample_id}._extract_1.fq.gz"
        R2="${sample_id}._extract_2.fq.gz"
    else
        # Direct extraction: detect chr naming (chr6 vs 6)
        CHR=\$(samtools view -H ${bam} | awk '/^@SQ.*SN:chr6\t/{print "chr6"; exit} /^@SQ.*SN:6\t/{print "6"; exit}')
        if [ -z "\$CHR" ]; then
            echo "[Kourami] WARNING: could not detect chr6 in BAM header for ${sample_id}"
            CHR="6"
        fi
        echo "[Kourami] Extracting HLA region from \${CHR}:28000000-34000000..."
        samtools view -b ${bam} "\${CHR}:28000000-34000000" | \
            samtools sort -n -@ ${task.cpus} | \
            samtools fastq -1 ${sample_id}._hla_1.fq.gz -2 ${sample_id}._hla_2.fq.gz -s /dev/null -
        R1="${sample_id}._hla_1.fq.gz"
        R2="${sample_id}._hla_2.fq.gz"
    fi

    if [ ! -f "\$R1" ] || [ ! -s "\$R1" ]; then
        echo "[Kourami] WARNING: HLA read extraction produced no output for ${sample_id}"
        printf "# Kourami results for ${sample_id}\n# WARNING: HLA read extraction failed\nGene\tAllele1\tAllele2\tReads1\tReads2\n" > ${sample_id}_kourami.txt
        exit 0
    fi

    # Step 2: Align extracted HLA reads to the augmented Kourami panel
    bwa mem -t ${task.cpus} \
        ${kourami_db}/All_FINAL_with_Decoy.fa.gz \
        \$R1 \$R2 \
        -o ${sample_id}.panel.sam

    samtools sort -@ ${task.cpus} -o ${sample_id}.panel.bam ${sample_id}.panel.sam
    samtools index ${sample_id}.panel.bam

    # Step 3: Kourami assembly-graph typing
    java -Xmx10g -jar ${kourami_dir}/build/Kourami.jar \
        -d ${kourami_db} \
        ${sample_id}.panel.bam \
        -o ${sample_id}.kourami || true

    # Parse .result file → standard pipeline TSV
    if [ -f "${sample_id}.kourami.result" ]; then
        python3 ${projectDir}/bin/parse_kourami_results.py \
            --input  ${sample_id}.kourami.result \
            --sample ${sample_id} \
            --output ${sample_id}_kourami.txt
    else
        printf "# Kourami results for ${sample_id}\n# WARNING: Kourami produced no output\nGene\tAllele1\tAllele2\tReads1\tReads2\n" > ${sample_id}_kourami.txt
    fi

    # Cleanup large intermediate files
    rm -f ${sample_id}.panel.sam ${sample_id}.panel.bam ${sample_id}.panel.bam.bai
    rm -f ${sample_id}._hla_1.fq.gz ${sample_id}._hla_2.fq.gz
    rm -f ${sample_id}._extract_1.fq.gz ${sample_id}._extract_2.fq.gz

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    kourami: "0.9.6"
	END_VERSIONS
    """
}
