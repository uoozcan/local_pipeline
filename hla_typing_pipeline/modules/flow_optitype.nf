/*
 * flow-OptiType Module
 * HLA Class I typing from BAM files using NMDP's OptiType workflow
 * Uses nmdpbioinformatics/flow-optitype Docker/Singularity container
 *
 * This implementation wraps the flow-OptiType workflow for direct BAM processing
 * Reference: https://github.com/nmdp-bioinformatics/flow-OptiType
 */

process FLOW_OPTITYPE {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/flow_optitype", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_flow_optitype.txt"), emit: results
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

    # flow-optitype expects BAM files in a directory
    # Create a temporary directory and link the BAM
    mkdir -p input_bams
    ln -s \$(readlink -f ${bam}) input_bams/${bam}
    if [ -f "${bam}.bai" ]; then
        ln -s \$(readlink -f ${bam}.bai) input_bams/${bam}.bai
    elif [ -f "${bam.baseName}.bai" ]; then
        ln -s \$(readlink -f ${bam.baseName}.bai) input_bams/${bam.baseName}.bai
    fi

    # Run OptiType via the flow-optitype container
    echo "[Running flow-OptiType on ${sample_id}...]"

    # The container includes OptiType configured for BAM processing
    # Extract reads from HLA region and run OptiType
    python /usr/local/bin/OptiTypePipeline.py \\
        -i ${bam} \\
        --dna \\
        -v \\
        -o ${sample_id} \\
        -p ${sample_id} \\
        2>&1 || {
        echo "Warning: flow-OptiType encountered an error, attempting fallback..."

        # Fallback: extract HLA reads and run OptiType manually
        samtools view -b ${bam} chr6:28510120-33480577 > hla_region.bam
        samtools sort -n -@ ${task.cpus} hla_region.bam -o sorted.bam
        samtools fastq -@ ${task.cpus} -1 R1.fastq -2 R2.fastq -0 /dev/null -s /dev/null sorted.bam

        python /usr/local/bin/OptiTypePipeline.py \\
            -i R1.fastq R2.fastq \\
            --dna \\
            -v \\
            -o ${sample_id} \\
            -p ${sample_id}
    }

    # Parse results to standard format
    echo "[Parsing flow-OptiType results...]"
    result_tsv=\$(find ${sample_id} -name "*_result.tsv" | head -1)

    if [ -n "\$result_tsv" ] && [ -f "\$result_tsv" ]; then
        python3 <<EOF
import csv
import sys

with open('\${result_tsv}', 'r') as f:
    reader = csv.DictReader(f, delimiter='\\t')
    row = next(reader)

with open('${sample_id}_flow_optitype.txt', 'w') as out:
    out.write("# flow-OptiType results for ${sample_id}\\n")
    out.write("# Tool: OptiType via flow-optitype (NMDP Bioinformatics)\\n")
    out.write("Gene\\tAllele1\\tAllele2\\tReads1\\tReads2\\n")
    out.write(f"A\\t{row.get('A1', 'NA')}\\t{row.get('A2', 'NA')}\\tNA\\tNA\\n")
    out.write(f"B\\t{row.get('B1', 'NA')}\\t{row.get('B2', 'NA')}\\tNA\\tNA\\n")
    out.write(f"C\\t{row.get('C1', 'NA')}\\t{row.get('C2', 'NA')}\\tNA\\tNA\\n")
EOF
    else
        echo "# flow-OptiType results for ${sample_id}" > ${sample_id}_flow_optitype.txt
        echo "# No results generated" >> ${sample_id}_flow_optitype.txt
        echo "Gene\tAllele1\tAllele2\tReads1\tReads2" >> ${sample_id}_flow_optitype.txt
        echo "A\tNA\tNA\tNA\tNA" >> ${sample_id}_flow_optitype.txt
        echo "B\tNA\tNA\tNA\tNA" >> ${sample_id}_flow_optitype.txt
        echo "C\tNA\tNA\tNA\tNA" >> ${sample_id}_flow_optitype.txt
    fi

    # Cleanup
    rm -rf input_bams hla_region.bam sorted.bam R1.fastq R2.fastq

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flow-optitype: "latest"
        optitype: \$(python /usr/local/bin/OptiTypePipeline.py --version 2>&1 | head -1 || echo "1.3.5")
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
