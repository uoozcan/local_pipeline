/*
 * ArcasHLA module (DSL2)
 * Fixes:
 *  - Uses TMPDIR within task workdir to avoid /tmp "No space left on device"
 *  - Handles single-end extraction outputs (*.extracted.fq[.gz]) as well as paired (*.extracted.1/2.fq.gz)
 *  - Uses params.arcashla_genes (so --hla_genes mapping works)
 *  - Normalises outputs to <sample_id>.genotype.json for downstream steps
 */
nextflow.enable.dsl = 2

process ARCASHLA_BAM {
    tag "$sample_id"
    label 'process_medium'
    cpus { params.arcashla_cpus ?: 24 }
    memory { params.arcashla_memory ?: '32.GB' }
    time { params.arcashla_time ?: '8.h' }

    // Don't silently "succeed" with empty outputs
    maxRetries 2
    errorStrategy 'terminate'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id),
              path("${sample_id}.genotype.json"),
              path("${sample_id}.genotype.log"),
              emit: results
        path("${sample_id}.extract.log"), optional: true, emit: logs

    when:
        params.tools?.toString()?.contains('arcashla')

    shell:
        '''
        set -euo pipefail

        # Make sure temporary files do NOT go to /tmp (can be too small inside containers)
        mkdir -p tmp
        export TMPDIR="$PWD/tmp"
        export TMP="$TMPDIR" TEMP="$TMPDIR" TMPDIR="$TMPDIR"

        BAM_BASE=$(basename "!{bam}")
        BAM_BASE=${BAM_BASE%.bam}
        BAM_BASE=${BAM_BASE%.cram}

        echo "[ARCASHLA_BAM] sample_id=!{sample_id}" >&2
        echo "[ARCASHLA_BAM] bam=!{bam}" >&2
        echo "[ARCASHLA_BAM] bam_base=${BAM_BASE}" >&2
        echo "[ARCASHLA_BAM] genes=!{params.arcashla_genes}" >&2
        echo "[ARCASHLA_BAM] TMPDIR=${TMPDIR}" >&2

        # 1) Extract HLA-related reads from BAM into FASTQ(s)
        arcasHLA extract -o . -t !{task.cpus} -v "!{bam}" --log "${BAM_BASE}.extract.log"

        # 2) Genotype from extracted FASTQ(s)
        if [[ -f "${BAM_BASE}.extracted.1.fq.gz" && -f "${BAM_BASE}.extracted.2.fq.gz" ]]; then
            arcasHLA genotype -o . -t !{task.cpus} -v -g "!{params.arcashla_genes}" \
                "${BAM_BASE}.extracted.1.fq.gz" "${BAM_BASE}.extracted.2.fq.gz" \
                --log "${BAM_BASE}.genotype.log"

        elif [[ -f "${BAM_BASE}.extracted.fq.gz" ]]; then
            arcasHLA genotype -o . -t !{task.cpus} -v -g "!{params.arcashla_genes}" \
                "${BAM_BASE}.extracted.fq.gz" \
                --log "${BAM_BASE}.genotype.log"

        elif [[ -f "${BAM_BASE}.extracted.fq" ]]; then
            gzip -f "${BAM_BASE}.extracted.fq"
            arcasHLA genotype -o . -t !{task.cpus} -v -g "!{params.arcashla_genes}" \
                "${BAM_BASE}.extracted.fq.gz" \
                --log "${BAM_BASE}.genotype.log"

        else
            echo "ERROR: arcasHLA extract did not produce expected FASTQ outputs." >&2
            echo "Contents of workdir:" >&2
            ls -lah >&2 || true
            exit 2
        fi

        # 3) Sanity check + normalise filenames
        # ArcasHLA may use a shortened prefix (e.g., HG00096 instead of full BAM basename)
        # Find the actual genotype.json file created
        GENOTYPE_JSON=$(ls -1 *.genotype.json 2>/dev/null | head -1)
        if [[ -z "$GENOTYPE_JSON" || ! -s "$GENOTYPE_JSON" ]]; then
            echo "ERROR: No genotype.json file found" >&2
            ls -lah >&2 || true
            exit 1
        fi
        GENOTYPE_PREFIX=${GENOTYPE_JSON%.genotype.json}

        cp -f "${GENOTYPE_PREFIX}.genotype.json" "!{sample_id}.genotype.json"
        cp -f "${GENOTYPE_PREFIX}.genotype.log"  "!{sample_id}.genotype.log" 2>/dev/null || \
            cp -f "${BAM_BASE}.genotype.log"     "!{sample_id}.genotype.log" 2>/dev/null || true
        cp -f "${BAM_BASE}.extract.log"          "!{sample_id}.extract.log"  || true
        '''
}

process ARCASHLA_FASTQ {
    tag "$sample_id"
    label 'process_medium'
    cpus { params.arcashla_cpus ?: 24 }
    memory { params.arcashla_memory ?: '32.GB' }
    time { params.arcashla_time ?: '8.h' }

    maxRetries 2
    errorStrategy 'terminate'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id),
              path("${sample_id}.genotype.json"),
              path("${sample_id}.genotype.log"),
              emit: results

    when:
        params.tools?.toString()?.contains('arcashla')

    shell:
        '''
        set -euo pipefail
        mkdir -p tmp
        export TMPDIR="$PWD/tmp"
        export TMP="$TMPDIR" TEMP="$TMPDIR" TMPDIR="$TMPDIR"

        echo "[ARCASHLA_FASTQ] sample_id=!{sample_id}" >&2
        echo "[ARCASHLA_FASTQ] r1=!{read1}" >&2
        echo "[ARCASHLA_FASTQ] r2=!{read2}" >&2
        echo "[ARCASHLA_FASTQ] genes=!{params.arcashla_genes}" >&2

        arcasHLA genotype -o . -t !{task.cpus} -v -g "!{params.arcashla_genes}" \
            "!{read1}" "!{read2}" \
            --log "!{sample_id}.genotype.log"

        test -s "!{sample_id}.genotype.json"
        '''
}
