/*
 * xHLA Module
 * HLA typing using Human Longevity Inc.'s xHLA algorithm
 * Supports BAM input (BAM recommended, FASTQ has limited support)
 * Uses k-mer based approach for HLA typing
 *
 * Note: xHLA was designed for hg19. For hg38 BAMs, we extract reads
 * from the hg38 HLA region and run xHLA's alignment/typing steps manually.
 */

process XHLA {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/xhla", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_xhla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index ${bam}
    fi

    # Determine chromosome naming convention (chr6 vs 6)
    CHR_PREFIX=\$(samtools view -H ${bam} | grep -m1 "^@SQ" | grep -o "SN:[^	]*" | cut -d: -f2 | grep -o "^chr" || echo "")

    # Detect reference genome version (hg19 vs hg38) based on chromosome 6 length
    CHR6_LEN=\$(samtools view -H ${bam} | grep "^@SQ.*SN:\${CHR_PREFIX}6" | grep -o "LN:[0-9]*" | cut -d: -f2)

    if [[ "\${CHR6_LEN}" =~ ^[0-9]+\$ ]] && [ "\${CHR6_LEN}" -gt "171000000" ]; then
        REF_VERSION="hg38"
        # hg38 HLA region
        HLA_REGION="\${CHR_PREFIX}6:28510120-33480577"
    else
        REF_VERSION="hg19"
        # hg19 HLA region (xHLA default)
        HLA_REGION="\${CHR_PREFIX}6:29886751-33090696"
    fi

    echo "Detected reference: \${REF_VERSION}, HLA region: \${HLA_REGION}"

    # Try running xHLA directly first (works best with hg19 + chr prefix)
    XHLA_SUCCESS=false
    if [ "\${REF_VERSION}" == "hg19" ] && [ "\${CHR_PREFIX}" == "chr" ]; then
        echo "[Running xHLA directly...]"
        python /opt/bin/run.py \
            --sample_id ${sample_id} \
            --input_bam_path ${bam} \
            --output_path ${sample_id}/ 2>&1 || true

        if [ -f "${sample_id}/report-${sample_id}-hla.json" ]; then
            XHLA_SUCCESS=true
        fi
    fi

    # If direct run failed, try manual extraction with correct coordinates
    if [ "\${XHLA_SUCCESS}" != "true" ]; then
        echo "[Running xHLA with manual extraction for \${REF_VERSION}...]"

        # Extract reads from HLA region
        samtools view ${bam} \${HLA_REGION} > ${sample_id}/temp.sam 2>/dev/null || true

        if [ -s "${sample_id}/temp.sam" ]; then
            # Convert to FASTQ using xHLA's preprocessor
            /opt/bin/preprocess.pl ${sample_id}/temp.sam | gzip > ${sample_id}/${sample_id}.fq.gz
            rm -f ${sample_id}/temp.sam

            # Run alignment
            /opt/bin/align.pl ${sample_id}/${sample_id}.fq.gz ${sample_id}/${sample_id}.tsv 2>&1 || true

            # Run typing with bash-based kill — timeout/SIGTERM ignored by R mclapply.
            # Run in background, poll every 10s, SIGKILL after 1800s.
            if [ -f "${sample_id}/${sample_id}.tsv" ]; then
                /opt/bin/typing.r ${sample_id}/${sample_id}.tsv ${sample_id}/${sample_id}.hla \
                    > ${sample_id}/typing_r.log 2>&1 &
                TYPING_PID=\$!
                echo "[typing.r started PID=\${TYPING_PID}]"
                WAITED=0
                while [ \${WAITED} -lt 1800 ]; do
                    kill -0 \${TYPING_PID} 2>/dev/null || { echo "[typing.r finished after \${WAITED}s]"; break; }
                    sleep 10
                    WAITED=\$((WAITED + 10))
                done
                if kill -0 \${TYPING_PID} 2>/dev/null; then
                    echo "[Killing typing.r after \${WAITED}s]"
                    kill -9 \${TYPING_PID} 2>/dev/null || true
                    # kill mclapply children
                    kill -9 \$(ps --ppid \${TYPING_PID} -o pid= 2>/dev/null) 2>/dev/null || true
                fi
                wait \${TYPING_PID} 2>/dev/null || true

                # Generate report if typing succeeded
                if [ -f "${sample_id}/${sample_id}.hla" ]; then
                    /opt/bin/report.py \
                        -in ${sample_id}/${sample_id}.hla \
                        -out ${sample_id}/report-${sample_id}-hla.json \
                        -subject ${sample_id} \
                        -sample ${sample_id} 2>&1 || true
                    XHLA_SUCCESS=true
                fi
            fi
        fi
    fi

    # Parse results to standard format
    echo "[Parsing xHLA results...]"
    if [ -f "${sample_id}/report-${sample_id}-hla.json" ]; then
        parse_xhla_results.py \
            --sample ${sample_id} \
            --input ${sample_id}/report-${sample_id}-hla.json \
            --output ${sample_id}_xhla.txt
    else
        # Try to recover allele calls from typing.r stdout (captured before OOM kill)
        RECOVERED=false
        if [ -f "${sample_id}/typing_r.log" ] && [ -s "${sample_id}/typing_r.log" ]; then
            echo "[Attempting allele recovery from typing.r output...]"
            grep -oE '[A-Z][A-Z0-9]*[*][0-9]+:[0-9]+' "${sample_id}/typing_r.log" | \
                grep -E '^(A|B|C|DRB1|DQA1|DQB1|DPA1|DPB1)[*]' | \
                sort -u > ${sample_id}/recovered_alleles.txt || true

            if [ -s "${sample_id}/recovered_alleles.txt" ]; then
                printf 'Gene\tAllele1\tAllele2\tReads1\tReads2\n' > ${sample_id}_xhla.txt
                for GENE in A B C DRB1 DQA1 DQB1 DPA1 DPB1; do
                    A1=\$(grep "^\${GENE}[*]" "${sample_id}/recovered_alleles.txt" | sed -n '1p' || true)
                    A2=\$(grep "^\${GENE}[*]" "${sample_id}/recovered_alleles.txt" | sed -n '2p' || true)
                    [ -z "\${A1}" ] && A1="NA" && A2="NA"
                    [ -z "\${A2}" ] && A2="\${A1}"
                    printf '%s\t%s\t%s\tNA\tNA\n' "\${GENE}" "\${A1}" "\${A2}"
                done >> ${sample_id}_xhla.txt
                RECOVERED=true
                echo "[Recovered alleles from typing.r output]"
            fi
        fi

        if [ "\${RECOVERED}" != "true" ]; then
            # Create empty results file if xHLA failed completely
            echo "# xHLA results for ${sample_id}" > ${sample_id}_xhla.txt
            echo "# xHLA typing failed - insufficient data or incompatible reference" >> ${sample_id}_xhla.txt
            printf 'Gene\tAllele1\tAllele2\tReads1\tReads2\n' >> ${sample_id}_xhla.txt
            printf 'A\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'B\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'C\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DRB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DQA1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DQB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DPA1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DPB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
        fi
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xhla: "1.0"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * xHLA from FASTQ files
 * Note: xHLA requires BAM/aligned reads. For FASTQ, we align to HLA reference first.
 * This is less reliable than BAM input.
 */
process XHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/xhla", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("${sample_id}_xhla.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}

    # xHLA works best with aligned BAM files
    # For FASTQ input, we'll process reads directly through xHLA's alignment

    echo "[Processing FASTQ for xHLA...]"

    # Combine and prepare FASTQ files
    if [[ "${fastq1}" == *.gz ]]; then
        zcat ${fastq1} ${fastq2} | gzip > ${sample_id}/${sample_id}.fq.gz
    else
        cat ${fastq1} ${fastq2} | gzip > ${sample_id}/${sample_id}.fq.gz
    fi

    # Run alignment
    echo "[Running xHLA alignment...]"
    /opt/bin/align.pl ${sample_id}/${sample_id}.fq.gz ${sample_id}/${sample_id}.tsv 2>&1 || true

    XHLA_SUCCESS=false

    # Run typing if alignment produced output
    if [ -f "${sample_id}/${sample_id}.tsv" ] && [ -s "${sample_id}/${sample_id}.tsv" ]; then
        echo "[Running xHLA typing...]"
        /opt/bin/typing.r ${sample_id}/${sample_id}.tsv ${sample_id}/${sample_id}.hla \
            > ${sample_id}/typing_r.log 2>&1 &
        TYPING_PID=\$!
        echo "[typing.r started PID=\${TYPING_PID}]"
        WAITED=0
        while [ \${WAITED} -lt 1800 ]; do
            kill -0 \${TYPING_PID} 2>/dev/null || { echo "[typing.r finished after \${WAITED}s]"; break; }
            sleep 10
            WAITED=\$((WAITED + 10))
        done
        if kill -0 \${TYPING_PID} 2>/dev/null; then
            echo "[Killing typing.r after \${WAITED}s]"
            kill -9 \${TYPING_PID} 2>/dev/null || true
            kill -9 \$(ps --ppid \${TYPING_PID} -o pid= 2>/dev/null) 2>/dev/null || true
        fi
        wait \${TYPING_PID} 2>/dev/null || true

        # Generate report if typing succeeded
        if [ -f "${sample_id}/${sample_id}.hla" ]; then
            /opt/bin/report.py \
                -in ${sample_id}/${sample_id}.hla \
                -out ${sample_id}/report-${sample_id}-hla.json \
                -subject ${sample_id} \
                -sample ${sample_id} 2>&1 || true
            XHLA_SUCCESS=true
        fi
    fi

    # Parse results to standard format
    echo "[Parsing xHLA results...]"
    if [ -f "${sample_id}/report-${sample_id}-hla.json" ]; then
        parse_xhla_results.py \
            --sample ${sample_id} \
            --input ${sample_id}/report-${sample_id}-hla.json \
            --output ${sample_id}_xhla.txt
    else
        # Try to recover allele calls from typing.r stdout (captured before OOM kill)
        RECOVERED=false
        if [ -f "${sample_id}/typing_r.log" ] && [ -s "${sample_id}/typing_r.log" ]; then
            echo "[Attempting allele recovery from typing.r output...]"
            grep -oE '[A-Z][A-Z0-9]*[*][0-9]+:[0-9]+' "${sample_id}/typing_r.log" | \
                grep -E '^(A|B|C|DRB1|DQA1|DQB1|DPA1|DPB1)[*]' | \
                sort -u > ${sample_id}/recovered_alleles.txt || true

            if [ -s "${sample_id}/recovered_alleles.txt" ]; then
                printf 'Gene\tAllele1\tAllele2\tReads1\tReads2\n' > ${sample_id}_xhla.txt
                for GENE in A B C DRB1 DQA1 DQB1 DPA1 DPB1; do
                    A1=\$(grep "^\${GENE}[*]" "${sample_id}/recovered_alleles.txt" | sed -n '1p' || true)
                    A2=\$(grep "^\${GENE}[*]" "${sample_id}/recovered_alleles.txt" | sed -n '2p' || true)
                    [ -z "\${A1}" ] && A1="NA" && A2="NA"
                    [ -z "\${A2}" ] && A2="\${A1}"
                    printf '%s\t%s\t%s\tNA\tNA\n' "\${GENE}" "\${A1}" "\${A2}"
                done >> ${sample_id}_xhla.txt
                RECOVERED=true
                echo "[Recovered alleles from typing.r output]"
            fi
        fi

        if [ "\${RECOVERED}" != "true" ]; then
            # Create empty results file if xHLA failed completely
            echo "# xHLA results for ${sample_id}" > ${sample_id}_xhla.txt
            echo "# xHLA typing failed - FASTQ input has limited support" >> ${sample_id}_xhla.txt
            printf 'Gene\tAllele1\tAllele2\tReads1\tReads2\n' >> ${sample_id}_xhla.txt
            printf 'A\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'B\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'C\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DRB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DQA1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DQB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DPA1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
            printf 'DPB1\tNA\tNA\tNA\tNA\n' >> ${sample_id}_xhla.txt
        fi
    fi

    # Cleanup large intermediate files
    rm -f ${sample_id}/${sample_id}.fq.gz

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xhla: "1.0"
    END_VERSIONS
    """
}
