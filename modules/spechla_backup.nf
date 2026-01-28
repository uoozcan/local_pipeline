// modules/spechla.nf

process SPECHLA_BAM {

    tag "${sample_id}"
    label 'medium'

    input:
    tuple val(sample_id), path(bam), path(bai, optional: true)

    output:
    path("${sample_id}_spechla.tsv"), emit: results
    path("${sample_id}_spechla.log"), emit: logs

    script:
    // Pipeline-facing parameters
    def genes = (params.spechla_genes ?: params.hla_genes ?: 'A,B,C,DQA1,DQB1,DRB1').toString()
    def exon  = (params.spechla_exon_only != null ? params.spechla_exon_only : 0) as Integer

    // SpecHLA uses -r hg19|hg38 typically; allow reference_build override, else infer from reference_genome string
    def ref_build = (params.reference_build ?: '').toString()
    if( !ref_build ) {
        def rg = (params.reference_genome ?: '').toString()
        ref_build = (rg.toLowerCase().contains('hg19') || rg.toLowerCase().contains('grch37')) ? 'hg19' : 'hg38'
    }

    // Optional explicit paths for SpecHLA install
    def spechla_home   = (params.spechla_home   ?: '').toString()
    def spechla_script = (params.spechla_script ?: '').toString()

    """
    set -euo pipefail

    echo "==================================================" > ${sample_id}_spechla.log
    echo "SpecHLA typing (BAM) - ${sample_id}" >> ${sample_id}_spechla.log
    echo "==================================================" >> ${sample_id}_spechla.log
    echo "BAM: ${bam}" >> ${sample_id}_spechla.log
    if [[ -n "${bai:-}" ]]; then echo "BAI: ${bai}" >> ${sample_id}_spechla.log; fi
    echo "genes: ${genes}" >> ${sample_id}_spechla.log
    echo "ref_build: ${ref_build}" >> ${sample_id}_spechla.log
    echo "exon_only (-u): ${exon}" >> ${sample_id}_spechla.log
    echo "cpus: ${task.cpus}" >> ${sample_id}_spechla.log
    echo "" >> ${sample_id}_spechla.log

    # Ensure BAM index exists
    if [[ -n "${bai:-}" && -s "${bai}" ]]; then
        ln -sf "${bai}" "${bam}.bai" 2>> ${sample_id}_spechla.log || true
    fi
    if [[ ! -s "${bam}.bai" && ! -s "${bam%.*}.bai" ]]; then
        echo "Index missing; building BAM index with samtools index" >> ${sample_id}_spechla.log
        samtools index -@ ${task.cpus} "${bam}" >> ${sample_id}_spechla.log 2>&1
    fi

    # Resolve SpecHLA entrypoint (SpecHLA.sh)
    SPECHLA_SH="${spechla_script}"
    SPECHLA_EXTRACT=""

    if [[ -z "$SPECHLA_SH" && -n "${spechla_home}" ]]; then
        [[ -f "${spechla_home}/script/whole/SpecHLA.sh" ]] && SPECHLA_SH="${spechla_home}/script/whole/SpecHLA.sh"
        [[ -f "${spechla_home}/script/ExtractHLAread.sh" ]] && SPECHLA_EXTRACT="${spechla_home}/script/ExtractHLAread.sh"
    fi

    if [[ -z "$SPECHLA_SH" ]]; then
        for p in \
            "script/whole/SpecHLA.sh" \
            "./script/whole/SpecHLA.sh" \
            "./SpecHLA/script/whole/SpecHLA.sh" \
            "/opt/SpecHLA/script/whole/SpecHLA.sh" \
            "/usr/local/SpecHLA/script/whole/SpecHLA.sh" \
            "/SpecHLA/script/whole/SpecHLA.sh"
        do
            if [[ -f "$p" ]]; then SPECHLA_SH="$p"; break; fi
        done
    fi

    if [[ -z "$SPECHLA_EXTRACT" ]]; then
        for p in \
            "script/ExtractHLAread.sh" \
            "./script/ExtractHLAread.sh" \
            "./SpecHLA/script/ExtractHLAread.sh" \
            "/opt/SpecHLA/script/ExtractHLAread.sh" \
            "/usr/local/SpecHLA/script/ExtractHLAread.sh" \
            "/SpecHLA/script/ExtractHLAread.sh"
        do
            if [[ -f "$p" ]]; then SPECHLA_EXTRACT="$p"; break; fi
        done
    fi

    if [[ -z "$SPECHLA_SH" ]]; then
        # IMPORTANT: escape $(...) for Groovy/Nextflow compilation
        if command -v SpecHLA.sh >/dev/null 2>&1; then
            SPECHLA_SH="\$(command -v SpecHLA.sh)"
        fi
    fi

    if [[ -z "$SPECHLA_SH" || ! -f "$SPECHLA_SH" ]]; then
        echo "ERROR: SpecHLA.sh not found. Provide --spechla_home or --spechla_script, or ensure it exists in the container." >> ${sample_id}_spechla.log
        exit 127
    fi

    echo "Using SpecHLA.sh: $SPECHLA_SH" >> ${sample_id}_spechla.log
    if [[ -n "$SPECHLA_EXTRACT" ]]; then
        echo "Using ExtractHLAread.sh: $SPECHLA_EXTRACT" >> ${sample_id}_spechla.log
    else
        echo "ExtractHLAread.sh not found; will derive FASTQ from BAM via samtools fastq" >> ${sample_id}_spechla.log
    fi
    echo "" >> ${sample_id}_spechla.log

    OUTDIR="spechla_out"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"

    # Prepare FASTQ inputs for SpecHLA
    R1=""
    R2=""

    if [[ -n "$SPECHLA_EXTRACT" ]]; then
        echo "Extracting HLA reads from BAM..." >> ${sample_id}_spechla.log
        bash "$SPECHLA_EXTRACT" -s "${sample_id}" -b "${bam}" -r "${ref_build}" -o "$OUTDIR" >> ${sample_id}_spechla.log 2>&1

        # Try common output patterns produced by ExtractHLAread.sh
        R1=\$(ls -1 "$OUTDIR"/*_R1*.fq* "$OUTDIR"/*_1*.fq* 2>/dev/null | head -n 1 || true)
        R2=\$(ls -1 "$OUTDIR"/*_R2*.fq* "$OUTDIR"/*_2*.fq* 2>/dev/null | head -n 1 || true)
    fi

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "Generating FASTQs from BAM (samtools fastq)..." >> ${sample_id}_spechla.log
        R1="$OUTDIR/${sample_id}.R1.fq.gz"
        R2="$OUTDIR/${sample_id}.R2.fq.gz"
        samtools fastq -@ ${task.cpus} -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -n "${bam}" >> ${sample_id}_spechla.log 2>&1
    fi

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "ERROR: failed to produce paired FASTQs for SpecHLA." >> ${sample_id}_spechla.log
        exit 2
    fi

    echo "FASTQ R1: $R1" >> ${sample_id}_spechla.log
    echo "FASTQ R2: $R2" >> ${sample_id}_spechla.log
    echo "" >> ${sample_id}_spechla.log

    # Run SpecHLA
    echo "Running SpecHLA..." >> ${sample_id}_spechla.log
    bash "$SPECHLA_SH" \
        -N "${sample_id}" \
        -1 "$R1" \
        -2 "$R2" \
        -o "$OUTDIR" \
        -r "${ref_build}" \
        -u ${exon} \
        -j ${task.cpus} \
        >> ${sample_id}_spechla.log 2>&1

    # Validate output
    RES="$OUTDIR/${sample_id}/hla.result.txt"
    if [[ ! -s "$RES" ]]; then
        # Some versions write directly under OUTDIR
        RES2="$OUTDIR/hla.result.txt"
        if [[ -s "$RES2" ]]; then RES="$RES2"; fi
    fi

    if [[ ! -s "$RES" ]]; then
        echo "ERROR: SpecHLA produced no hla.result.txt (or it is empty)." >> ${sample_id}_spechla.log
        printf "sample\ttool\tstatus\tmessage\n%s\tspechla\tFAIL\tMissing hla.result.txt\n" "${sample_id}" > ${sample_id}_spechla.tsv
        exit 3
    fi

    # Convert to a simple TSV (sample, gene, allele1, allele2, tool)
    # Keep only requested genes if possible
    python3 - << 'PY' > ${sample_id}_spechla.tsv
import re
sample = "${sample_id}"
genes = set([g.strip() for g in "${genes}".split(",") if g.strip()])
res_path = "${RES}"

rows = []
with open(res_path, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+|\t", line)
        if len(parts) < 3:
            continue
        gene = parts[0].replace("HLA-", "")
        if genes and gene not in genes:
            continue
        rows.append((sample, gene, parts[1], parts[2], "spechla"))

print("sample\tgene\tallele1\tallele2\ttool")
for r in rows:
    print("\t".join(r))
PY

    if [[ \$(wc -l < ${sample_id}_spechla.tsv) -le 1 ]]; then
        echo "ERROR: SpecHLA output parsed but yielded no gene calls (check hla.result.txt format)." >> ${sample_id}_spechla.log
        printf "sample\ttool\tstatus\tmessage\n%s\tspechla\tFAIL\tParsed zero calls\n" "${sample_id}" > ${sample_id}_spechla.tsv
        exit 4
    fi

    echo "OK: SpecHLA results written to ${sample_id}_spechla.tsv" >> ${sample_id}_spechla.log
    """
}

process SPECHLA_FASTQ {

    tag "${sample_id}"
    label 'medium'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("${sample_id}_spechla.tsv"), emit: results
    path("${sample_id}_spechla.log"), emit: logs

    script:
    def genes = (params.spechla_genes ?: params.hla_genes ?: 'A,B,C,DQA1,DQB1,DRB1').toString()
    def exon  = (params.spechla_exon_only != null ? params.spechla_exon_only : 0) as Integer

    def ref_build = (params.reference_build ?: '').toString()
    if( !ref_build ) {
        def rg = (params.reference_genome ?: '').toString()
        ref_build = (rg.toLowerCase().contains('hg19') || rg.toLowerCase().contains('grch37')) ? 'hg19' : 'hg38'
    }

    def spechla_home   = (params.spechla_home   ?: '').toString()
    def spechla_script = (params.spechla_script ?: '').toString()

    """
    set -euo pipefail

    echo "==================================================" > ${sample_id}_spechla.log
    echo "SpecHLA typing (FASTQ) - ${sample_id}" >> ${sample_id}_spechla.log
    echo "==================================================" >> ${sample_id}_spechla.log
    echo "R1: ${read1}" >> ${sample_id}_spechla.log
    echo "R2: ${read2}" >> ${sample_id}_spechla.log
    echo "genes: ${genes}" >> ${sample_id}_spechla.log
    echo "ref_build: ${ref_build}" >> ${sample_id}_spechla.log
    echo "exon_only (-u): ${exon}" >> ${sample_id}_spechla.log
    echo "cpus: ${task.cpus}" >> ${sample_id}_spechla.log
    echo "" >> ${sample_id}_spechla.log

    SPECHLA_SH="${spechla_script}"

    if [[ -z "$SPECHLA_SH" && -n "${spechla_home}" ]]; then
        [[ -f "${spechla_home}/script/whole/SpecHLA.sh" ]] && SPECHLA_SH="${spechla_home}/script/whole/SpecHLA.sh"
    fi

    if [[ -z "$SPECHLA_SH" ]]; then
        for p in \
            "script/whole/SpecHLA.sh" \
            "./script/whole/SpecHLA.sh" \
            "./SpecHLA/script/whole/SpecHLA.sh" \
            "/opt/SpecHLA/script/whole/SpecHLA.sh" \
            "/usr/local/SpecHLA/script/whole/SpecHLA.sh" \
            "/SpecHLA/script/whole/SpecHLA.sh"
        do
            if [[ -f "$p" ]]; then SPECHLA_SH="$p"; break; fi
        done
    fi

    if [[ -z "$SPECHLA_SH" ]]; then
        if command -v SpecHLA.sh >/dev/null 2>&1; then
            SPECHLA_SH="\$(command -v SpecHLA.sh)"
        fi
    fi

    if [[ -z "$SPECHLA_SH" || ! -f "$SPECHLA_SH" ]]; then
        echo "ERROR: SpecHLA.sh not found. Provide --spechla_home or --spechla_script, or ensure it exists in the container." >> ${sample_id}_spechla.log
        exit 127
    fi

    echo "Using SpecHLA.sh: $SPECHLA_SH" >> ${sample_id}_spechla.log

    OUTDIR="spechla_out"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"

    echo "Running SpecHLA..." >> ${sample_id}_spechla.log
    bash "$SPECHLA_SH" \
        -N "${sample_id}" \
        -1 "${read1}" \
        -2 "${read2}" \
        -o "$OUTDIR" \
        -r "${ref_build}" \
        -u ${exon} \
        -j ${task.cpus} \
        >> ${sample_id}_spechla.log 2>&1

    RES="$OUTDIR/${sample_id}/hla.result.txt"
    if [[ ! -s "$RES" ]]; then
        RES2="$OUTDIR/hla.result.txt"
        if [[ -s "$RES2" ]]; then RES="$RES2"; fi
    fi

    if [[ ! -s "$RES" ]]; then
        echo "ERROR: SpecHLA produced no hla.result.txt (or it is empty)." >> ${sample_id}_spechla.log
        printf "sample\ttool\tstatus\tmessage\n%s\tspechla\tFAIL\tMissing hla.result.txt\n" "${sample_id}" > ${sample_id}_spechla.tsv
        exit 3
    fi

    python3 - << 'PY' > ${sample_id}_spechla.tsv
import re
sample = "${sample_id}"
genes = set([g.strip() for g in "${genes}".split(",") if g.strip()])
res_path = "${RES}"

rows = []
with open(res_path, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+|\t", line)
        if len(parts) < 3:
            continue
        gene = parts[0].replace("HLA-", "")
        if genes and gene not in genes:
            continue
        rows.append((sample, gene, parts[1], parts[2], "spechla"))

print("sample\tgene\tallele1\tallele2\ttool")
for r in rows:
    print("\t".join(r))
PY

    if [[ \$(wc -l < ${sample_id}_spechla.tsv) -le 1 ]]; then
        echo "ERROR: SpecHLA output parsed but yielded no gene calls." >> ${sample_id}_spechla.log
        printf "sample\ttool\tstatus\tmessage\n%s\tspechla\tFAIL\tParsed zero calls\n" "${sample_id}" > ${sample_id}_spechla.tsv
        exit 4
    fi

    echo "OK: SpecHLA results written to ${sample_id}_spechla.tsv" >> ${sample_id}_spechla.log
    """
}
