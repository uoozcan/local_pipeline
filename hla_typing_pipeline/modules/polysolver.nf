/*
 * POLYSOLVER Module
 * HLA Class I typing optimised for tumor/normal BAMs
 * Broad Institute tool widely used in cancer genomics (AML and somatic studies)
 * BAM input only — hg19/GRCh37 OR hg38 aligned coordinate-sorted BAM
 *
 * Container: docker://sachet/polysolver:v4  (pull: singularity pull polysolver.sif docker://sachet/polysolver:v4)
 *
 * KNOWN ISSUE on WSL2 / modern kernels (tested Mar 2026):
 *   The novoalign binary inside the container (v2.07.18, built 2012) crashes with SIGSEGV
 *   on Linux kernels > ~4.x (prints "Interrupted..11 Stack Dump").
 *   POLYSOLVER will NOT produce results in this environment.
 *   Expected to work on native Linux (CentOS 7, Ubuntu 18.04) or HPC nodes.
 *
 * hg38 fix: SAMTOOLS_DIR env var is unset in the container for hg38 branch.
 *   Pass SINGULARITYENV_SAMTOOLS_DIR=/home/polysolver/binaries in run.config containerOptions.
 */

process POLYSOLVER {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/polysolver", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_polysolver.txt"), emit: results
    path("${sample_id}_polysolver_raw/winners.hla.nofreq.txt"), emit: raw, optional: true
    path "versions.yml", emit: versions

    script:
    def build = params.reference == 'hg38' ? 'hg38' : 'hg19'
    """
    mkdir -p ${sample_id}_polysolver_raw

    echo "[POLYSOLVER] Running on ${sample_id} (build: ${build})..."

    # Ensure BAM index exists
    [ -f "${bam}.bai" ] || samtools index ${bam}

    # Fix for hg38 branch: SAMTOOLS_DIR env var is unset in container
    export SAMTOOLS_DIR=/home/polysolver/binaries

    # POLYSOLVER args: BAM race includeFreq build format insertCalc outdir
    # race=Unknown (population-agnostic), includeFreq=0, insertCalc=0 (germline)
    bash /home/polysolver/scripts/shell_call_hla_type \
        ${bam} Unknown 0 ${build} STDFQ 0 ${sample_id}_polysolver_raw || true

    # Parse winners.hla.nofreq.txt → standard pipeline TSV
    if [ -f "${sample_id}_polysolver_raw/winners.hla.nofreq.txt" ]; then
        python3 ${projectDir}/bin/parse_polysolver_results.py \
            --input  ${sample_id}_polysolver_raw/winners.hla.nofreq.txt \
            --sample ${sample_id} \
            --output ${sample_id}_polysolver.txt
    else
        echo "# POLYSOLVER results for ${sample_id}" > ${sample_id}_polysolver.txt
        echo "# WARNING: POLYSOLVER produced no output" >> ${sample_id}_polysolver.txt
        echo "Gene\tAllele1\tAllele2\tReads1\tReads2" >> ${sample_id}_polysolver.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: "v4"
    END_VERSIONS
    """
}
