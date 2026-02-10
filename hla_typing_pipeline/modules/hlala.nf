/*
 * HLA*LA Module
 * Graph-based HLA typing from WGS/WES data
 * Note: Requires ~16GB+ RAM for mapAgainstCompleteGenome option
 */

process HLALA {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/hlala", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val graph

    output:
    tuple val(sample_id), path("${sample_id}_hlala.txt"), emit: results
    tuple val(sample_id), path("${sample_id}/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    def graph_name = graph ?: 'PRG_MHC_GRCh38_withIMGT'
    """
    # Create output directory
    mkdir -p ${sample_id}

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Run HLA*LA
    echo "[Running HLA*LA...]"
    perl /usr/local/bin/HLA-LA/src/HLA-LA.pl \
        --BAM ${bam} \
        --graph ${graph_name} \
        --sampleID ${sample_id} \
        --maxThreads ${task.cpus} \
        --workingDir .

    # Parse results
    echo "[Parsing results...]"
    RESULT_FILE="${sample_id}/hla/R1_bestguess_G.txt"
    if [ -f "\$RESULT_FILE" ]; then
        # Convert HLA*LA format to standard format with quality scores
        python3 <<EOF
import sys
import re

with open('\$RESULT_FILE', 'r') as f:
    lines = f.readlines()

with open('${sample_id}_hlala.txt', 'w') as out:
    out.write("# HLA*LA results for ${sample_id}\\n")
    out.write("# Graph: ${graph_name}\\n")
    out.write("Gene\\tAllele1\\tAllele2\\tReads1\\tReads2\\tQuality\\n")

    for line in lines:
        line = line.strip()
        if not line or line.startswith('Locus'):
            continue
        parts = line.split('\\t')
        if len(parts) >= 3:
            gene = re.sub(r'^HLA-', '', parts[0])
            a1 = re.sub(r'^HLA-', '', parts[1]) if parts[1] != '?' else 'NA'
            a2 = re.sub(r'^HLA-', '', parts[2]) if parts[2] != '?' else 'NA'

            # Extract confidence score from remaining columns
            quality = 'NA'
            for col in parts[3:]:
                try:
                    q = float(col)
                    if 0.0 <= q <= 1.0:
                        quality = col
                        break
                except ValueError:
                    continue

            out.write(f"{gene}\\t{a1}\\t{a2}\\tNA\\tNA\\t{quality}\\n")
EOF
    else
        echo "# HLA*LA results for ${sample_id}" > ${sample_id}_hlala.txt
        echo "# HLA*LA typing failed - no results generated" >> ${sample_id}_hlala.txt
        echo "Gene\tAllele1\tAllele2\tReads1\tReads2\tQuality" >> ${sample_id}_hlala.txt
        echo "A\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "B\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "C\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "DRB1\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "DQA1\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "DQB1\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "DPA1\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
        echo "DPB1\tNA\tNA\tNA\tNA\tNA" >> ${sample_id}_hlala.txt
    fi

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlala: "1.0.3"
    END_VERSIONS
    """
}
