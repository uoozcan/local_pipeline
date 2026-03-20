/*
 * POLYSOLVER Module
 * HLA Class I typing optimised for tumor/normal BAMs
 * Broad Institute tool widely used in cancer genomics (AML and somatic studies)
 * BAM input only — hg19/GRCh37 OR hg38 aligned coordinate-sorted BAM
 *
 * Container: docker://sachet/polysolver:v4  (pull: singularity pull polysolver.sif docker://sachet/polysolver:v4)
 *
 * Fixes applied (effective on native Linux / Puhti; novoalign does NOT work on WSL2):
 *   1. novoalign SIGSEGV: `ulimit -s unlimited` added; helps on some systems with small default stacks
 *      NOTE: novoalign (v2 and v3) crashes with SIGSEGV on WSL2 kernel 6.6 regardless of ulimit —
 *      WSL2-specific kernel incompatibility; works on native Linux (CentOS 7, Ubuntu 18.04+) and Puhti
 *   2. hg38 SAMTOOLS_DIR: env var unset in container → `export SAMTOOLS_DIR=/home/polysolver/binaries`
 *   3. hg38 Picard "Illegal mate state": fixmate pre-processing fixes inconsistent mate flags
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

    # Fix 1: novoalign stack overflow on modern kernels (default 8MB stack too small)
    ulimit -s unlimited

    # Fix 2: SAMTOOLS_DIR env var is unset in container for hg38 branch
    export SAMTOOLS_DIR=/home/polysolver/binaries

    # Fix 3: Picard SamToFastq "Illegal mate state" — pre-sort by name + fixmate + re-sort
    # Required when BAM has inconsistent mate flags (common in some pipelines)
    echo "[POLYSOLVER] Applying fixmate pre-processing for ${sample_id}..."
    /home/polysolver/binaries/samtools sort -n -@ ${task.cpus} -o ${sample_id}_namesort.bam ${bam}
    /home/polysolver/binaries/samtools fixmate -m ${sample_id}_namesort.bam ${sample_id}_fixmate.bam
    /home/polysolver/binaries/samtools sort    -@ ${task.cpus} -o ${sample_id}_fixed.bam ${sample_id}_fixmate.bam
    /home/polysolver/binaries/samtools index ${sample_id}_fixed.bam
    POLYSOLVER_INPUT="${sample_id}_fixed.bam"

    # POLYSOLVER args: BAM race includeFreq build format insertCalc outdir
    # race=Unknown (population-agnostic), includeFreq=0, insertCalc=0 (germline)
    bash /home/polysolver/scripts/shell_call_hla_type \
        \$POLYSOLVER_INPUT Unknown 0 ${build} STDFQ 0 ${sample_id}_polysolver_raw || true

    # Cleanup intermediate BAMs
    rm -f ${sample_id}_namesort.bam ${sample_id}_fixmate.bam ${sample_id}_fixed.bam ${sample_id}_fixed.bam.bai

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
