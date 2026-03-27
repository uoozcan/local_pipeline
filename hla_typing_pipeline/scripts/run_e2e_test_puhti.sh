#!/bin/bash
# =============================================================================
# run_e2e_test_puhti.sh
# Single-sample end-to-end pipeline test on CSC Puhti (WGS / WES / RNA-seq)
#
# Downloads one sample, runs all typing tools via Nextflow, collects results,
# and reports per-tool concordance against the Gourraud 2014 ground truth.
# Use this to verify the pipeline works end-to-end before full calibration.
#
# Usage:
#   bash scripts/run_e2e_test_puhti.sh --type wgs|wes|rna [OPTIONS]
#
# Options:
#   --type   TYPE    Analysis type: wgs, wes, or rna  (REQUIRED)
#   --project ID     CSC project account (default: project_2008084)
#   --sample SAMPLE  1KGP sample ID to test (defaults per type below)
#   --tools  TOOLS   Comma-separated tools (default: type-specific, no kourami)
#   --skip-download  Skip Phase 0 (input data already on scratch)
#   --skip-typing    Skip Phase 1 (Nextflow already ran)
#   --status         Show job/output status and exit
#   --dry-run        Print SLURM commands without submitting
#   -h, --help       Show this help
#
# Default test samples (all in Gourraud 2014 GT):
#   WGS: NA19238  (YRI, NYGC 30x CRAM from PRJEB31736)
#   WES: NA18501  (YRI, 1KGP Phase 3 WES BAM from EBI FTP)
#   RNA: NA18502  (YRI, Geuvadis ERP001942 paired FASTQs)
#
# Default tools per type (Kourami excluded — takes ~8h):
#   WGS: hlahd,spechla,arcashla,optitype,polysolver
#   WES: hlahd,spechla,arcashla,optitype,polysolver
#   RNA: hlahd,spechla,arcashla,optitype,t1k,seq2hla
#
# Output:
#   Phase 0 log:  ${LOGS_DIR}/test_phase0.out
#   Phase 1 log:  ${LOGS_DIR}/test_typing.out
#   Test report:  ${SCRATCH_BASE}/report_${SAMPLE}.txt
#
# Expected run time (small queue, no waiting):
#   WGS: ~3h extract + ~4h typing + ~5min report
#   WES: ~2h extract + ~3h typing + ~5min report
#   RNA: ~2h download + ~3h typing + ~5min report
# =============================================================================
set -euo pipefail

#-----------------------------------------------------------------------------
# Defaults
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALL_DIR="$(dirname "$SCRIPT_DIR")"   # hla_typing_pipeline/

TYPE=""
PROJECT_ID="${SLURM_JOB_ACCOUNT:-project_2008084}"
SAMPLE=""
TOOLS=""
REFERENCE=""
SKIP_DOWNLOAD=false
SKIP_TYPING=false
STATUS_ONLY=false
DRY_RUN=false

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --type)          TYPE="$2";            shift 2 ;;
        --project)       PROJECT_ID="$2";      shift 2 ;;
        --sample)        SAMPLE="$2";          shift 2 ;;
        --tools)         TOOLS="$2";           shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true;   shift ;;
        --skip-typing)   SKIP_TYPING=true;     shift ;;
        --status)        STATUS_ONLY=true;     shift ;;
        --dry-run)       DRY_RUN=true;         shift ;;
        -h|--help)
            sed -n '2,47p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$TYPE" ]]; then
    echo "[ERROR] --type wgs|wes|rna is required" >&2
    exit 1
fi

#-----------------------------------------------------------------------------
# Per-type configuration
#-----------------------------------------------------------------------------
case "$TYPE" in
    wgs)
        [[ -z "$SAMPLE" ]] && SAMPLE="NA19238"
        [[ -z "$TOOLS"  ]] && TOOLS="hlahd,spechla,arcashla,optitype,polysolver"
        SEQ_TYPE="dna"
        REFERENCE="hg38"
        HLA_REGION="chr6:28000000-34000000"
        INPUT_MODE="bam"   # extract CRAM → BAM, feed BAM to Nextflow
        ENA_PROJECT="PRJEB31736"
        ENA_META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=run_accession,submitted_ftp&format=tsv"
        PHASE0_TIME="03:00:00"
        PHASE0_JOB_NAME="test_wgs_extr"
        ;;
    wes)
        [[ -z "$SAMPLE" ]] && SAMPLE="NA18501"
        [[ -z "$TOOLS"  ]] && TOOLS="hlahd,spechla,arcashla,optitype,polysolver"
        SEQ_TYPE="wes"
        REFERENCE="hg19"
        HLA_REGION="6:28000000-34000000"   # ENSEMBL naming (no chr prefix)
        INPUT_MODE="bam"
        PHASE0_TIME="02:00:00"
        PHASE0_JOB_NAME="test_wes_extr"
        ;;
    rna)
        [[ -z "$SAMPLE" ]] && SAMPLE="NA18502"
        [[ -z "$TOOLS"  ]] && TOOLS="hlahd,spechla,arcashla,optitype,t1k,seq2hla"
        SEQ_TYPE="rna"
        REFERENCE="hg38"
        INPUT_MODE="fastq"
        ENA_PROJECT="ERP001942"
        ENA_META_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=sample_alias,sample_accession,run_accession,fastq_ftp&format=tsv"
        PHASE0_TIME="04:00:00"
        PHASE0_JOB_NAME="test_rna_dl"
        ;;
    *)
        echo "[ERROR] Unknown type '${TYPE}'. Use: wgs, wes, or rna" >&2
        exit 1
        ;;
esac

#-----------------------------------------------------------------------------
# Derived paths
#-----------------------------------------------------------------------------
SCRATCH_BASE="/scratch/${PROJECT_ID}/hla_calibration/test_${TYPE}"
BAM_DIR="${SCRATCH_BASE}/bams"
FASTQ_DIR="${SCRATCH_BASE}/fastqs"
RESULTS_DIR="${SCRATCH_BASE}/results"
BY_TOOL_DIR="${RESULTS_DIR}/by_tool"
INDEX_DIR="${SCRATCH_BASE}/index"
LOGS_DIR="${SCRATCH_BASE}/logs"
CRAM_REF_CACHE="${SCRATCH_BASE}/cram_ref_cache"
GT_FILE="/scratch/${PROJECT_ID}/hla_tools/hla_typing_pipeline/conf/1kgp_hla_gt.tsv"
ENA_META="${INDEX_DIR}/ena_run_table_${TYPE}.tsv"
URL_MAP="${INDEX_DIR}/sample_url_map_${TYPE}.tsv"
SAMPLESHEET="${SCRATCH_BASE}/test_samplesheet.csv"
REPORT_FILE="${SCRATCH_BASE}/report_${SAMPLE}.txt"
SAMPLE_LIST="${SCRATCH_BASE}/test_sample.txt"

if [[ "$DRY_RUN" == "false" ]]; then
    mkdir -p "${SCRATCH_BASE}" "${BAM_DIR}" "${FASTQ_DIR}" "${RESULTS_DIR}" \
             "${INDEX_DIR}" "${LOGS_DIR}" "${CRAM_REF_CACHE}"
    echo "${SAMPLE}" > "${SAMPLE_LIST}"
    SCRIPTS_DIR="${SCRATCH_BASE}"
else
    # In dry-run, write generated SLURM scripts to a temp dir
    SCRIPTS_DIR=$(mktemp -d)
    LOGS_DIR="${SCRIPTS_DIR}/logs"
    mkdir -p "${LOGS_DIR}"
fi

#-----------------------------------------------------------------------------
# --status: show job and output status then exit
#-----------------------------------------------------------------------------
if [[ "$STATUS_ONLY" == "true" ]]; then
    echo "=== E2E Test Status: ${TYPE^^} / ${SAMPLE} ==="
    echo "Scratch:  ${SCRATCH_BASE}"
    echo ""
    if [[ "$INPUT_MODE" == "bam" ]]; then
        BAM_OUT="${BAM_DIR}/${SAMPLE}_hla.bam"
        [[ -f "$BAM_OUT" ]] && echo "Phase 0: HLA BAM present ($(du -sh "$BAM_OUT" | cut -f1))" \
                            || echo "Phase 0: HLA BAM NOT present"
    else
        R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
        ([[ -f "$R1" ]] && [[ -f "$R2" ]]) \
            && echo "Phase 0: FASTQs present ($(du -sh "$R1" "$R2" | awk '{s+=$1}END{print s}') MB total)" \
            || echo "Phase 0: FASTQs NOT present"
    fi
    echo ""
    echo "Phase 1: Typing results:"
    for TOOL in $(echo "$TOOLS" | tr ',' ' '); do
        RES="${RESULTS_DIR}/${SAMPLE}/${TOOL}/${SAMPLE}_${TOOL}.txt"
        [[ -f "$RES" ]] \
            && printf "  %-12s PRESENT  %s\n" "$TOOL" "$RES" \
            || printf "  %-12s missing\n" "$TOOL"
    done
    echo ""
    [[ -f "$REPORT_FILE" ]] \
        && { echo "Test report: ${REPORT_FILE}"; echo ""; cat "$REPORT_FILE"; } \
        || echo "Test report: NOT YET GENERATED"
    exit 0
fi

echo "==================================================================="
echo " E2E Pipeline Test — CSC Puhti"
echo "==================================================================="
echo " Type:    ${TYPE^^}"
echo " Sample:  ${SAMPLE}"
echo " Tools:   ${TOOLS}"
echo " Scratch: ${SCRATCH_BASE}"
echo " Dry-run: ${DRY_RUN}"
echo "==================================================================="
echo ""

#-----------------------------------------------------------------------------
# Build URL map (type-specific, same logic as calibration scripts)
#-----------------------------------------------------------------------------
if [[ "$SKIP_DOWNLOAD" == "false" ]]; then

    # Download ENA metadata if needed (WGS/RNA only; WES uses the population panel)
    if [[ "$TYPE" != "wes" ]]; then
        if [[ ! -f "${ENA_META}" ]] && [[ "$DRY_RUN" == "false" ]]; then
            echo "[INFO] Downloading ENA metadata for ${ENA_PROJECT}..."
            curl -s "${ENA_META_URL}" > "${ENA_META}" \
                || wget -q -O "${ENA_META}" "${ENA_META_URL}" \
                || { echo "[ERROR] Could not download ENA run table"; exit 1; }
            echo "[OK] ENA metadata: ${ENA_META} ($(wc -l < "${ENA_META}") entries)"
        else
            [[ "$DRY_RUN" == "false" ]] && echo "[INFO] Using cached ENA metadata: ${ENA_META}"
        fi
    fi

    # Build URL map (type-specific)
    if [[ "$DRY_RUN" == "false" ]]; then

        if [[ "$TYPE" == "wgs" ]]; then
            echo "[INFO] Building CRAM URL map for ${SAMPLE}..."
            python3 - << PYEOF
import sys, re
meta_file   = "${ENA_META}"
url_map_out = "${URL_MAP}"
sample      = "${SAMPLE}"
sample_map = {}
with open(meta_file) as fh:
    header = None
    for line in fh:
        line = line.rstrip("\n")
        if not line: continue
        parts = line.split("\t")
        if header is None:
            header = parts
            try:
                i_run = header.index("run_accession")
                i_ftp = header.index("submitted_ftp")
            except ValueError:
                print(f"[ERROR] Unexpected ENA header: {header}", file=sys.stderr); sys.exit(1)
            continue
        if len(parts) <= max(i_run, i_ftp): continue
        err = parts[i_run].strip()
        ftp = parts[i_ftp].strip()
        m = re.search(r'/([A-Z0-9]+)\.final\.cram(?:;|\$)', ftp)
        if not m: continue
        s = m.group(1)
        if s != sample: continue
        cram_url = f"https://ftp.sra.ebi.ac.uk/vol1/run/{err[:6]}/{err}/{s}.final.cram"
        if s not in sample_map:
            sample_map[s] = (err, cram_url)
with open(url_map_out, "w") as fout:
    if sample in sample_map:
        err, url = sample_map[sample]
        fout.write(f"{sample}\t{err}\t{url}\n")
        print(f"[OK] CRAM URL: {url}")
    else:
        print(f"[ERROR] {sample} not found in PRJEB31736 ENA metadata", file=sys.stderr)
        sys.exit(1)
PYEOF

        elif [[ "$TYPE" == "wes" ]]; then
            echo "[INFO] Building WES BAM URL map for ${SAMPLE}..."
            # Download population panel for URL construction
            PANEL_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
            PANEL_FILE="${INDEX_DIR}/1kgp_panel.tsv"
            if [[ ! -f "${PANEL_FILE}" ]]; then
                wget -q -O "${PANEL_FILE}" "${PANEL_URL}" \
                    || { echo "[ERROR] Could not download population panel"; exit 1; }
            fi
            python3 - << PYEOF
import sys
panel_file  = "${PANEL_FILE}"
url_map_out = "${URL_MAP}"
sample      = "${SAMPLE}"
ftp_base    = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
DATES       = ["20121211", "20130415", "20120522"]
pop_map = {}
with open(panel_file) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("sample"): continue
        parts = line.split()
        if len(parts) >= 2:
            pop_map[parts[0]] = parts[1]
pop = pop_map.get(sample)
if not pop:
    print(f"[ERROR] {sample} not in population panel", file=sys.stderr); sys.exit(1)
url = f"{ftp_base}/{sample}/exome_alignment/{sample}.mapped.ILLUMINA.bwa.{pop}.exome.{DATES[0]}.bam"
alts = [f"{ftp_base}/{sample}/exome_alignment/{sample}.mapped.ILLUMINA.bwa.{pop}.exome.{d}.bam"
        for d in DATES[1:]]
with open(url_map_out, "w") as fout:
    fout.write(f"{sample}\t{url}\t{'|'.join(alts)}\n")
print(f"[OK] WES BAM URL: {url}")
PYEOF

        elif [[ "$TYPE" == "rna" ]]; then
            echo "[INFO] Building RNA FASTQ URL map for ${SAMPLE}..."
            python3 - << PYEOF
import sys
meta_file   = "${ENA_META}"
url_map_out = "${URL_MAP}"
sample      = "${SAMPLE}"
sample_runs = {}
with open(meta_file) as fh:
    header = None
    for line in fh:
        line = line.rstrip("\n")
        if not line: continue
        parts = line.split("\t")
        if header is None:
            header = parts
            try:
                i_sample = header.index("sample_alias")
                i_ftp    = header.index("fastq_ftp")
            except ValueError:
                print(f"[ERROR] Unexpected ENA header: {header}", file=sys.stderr); sys.exit(1)
            continue
        if len(parts) <= max(i_sample, i_ftp): continue
        sid = parts[i_sample].strip()
        if ':' in sid: sid = sid.split(':', 1)[1]
        if sid != sample: continue
        ftp_paths = [p.strip() for p in parts[i_ftp].split(";") if p.strip()]
        r1 = next((p for p in ftp_paths if "_1.fastq.gz" in p), None)
        r2 = next((p for p in ftp_paths if "_2.fastq.gz" in p), None)
        if r1 and r2:
            if sid not in sample_runs: sample_runs[sid] = []
            sample_runs[sid].append((r1, r2))
runs = sample_runs.get(sample)
if not runs:
    print(f"[ERROR] {sample} not found in Geuvadis ENA metadata", file=sys.stderr); sys.exit(1)
r1_primary, r2_primary = runs[0]
r1_alts = "|".join(r for r, _ in runs[1:]) if len(runs) > 1 else ""
r2_alts = "|".join(r for _, r in runs[1:]) if len(runs) > 1 else ""
with open(url_map_out, "w") as fout:
    fout.write(f"{sample}\t{r1_primary}\t{r2_primary}\t{r1_alts}\t{r2_alts}\n")
print(f"[OK] FASTQ URLs: R1={r1_primary}")
PYEOF
        fi
    else
        echo "[DRY-RUN] Would build URL map from ENA metadata → ${URL_MAP}"
    fi
fi   # end SKIP_DOWNLOAD==false URL-map block

#-----------------------------------------------------------------------------
# Phase 0: SLURM job — download/extract input data for the test sample
#-----------------------------------------------------------------------------
PHASE0_DEP=""

if [[ "$SKIP_DOWNLOAD" == "false" ]]; then
    echo ""
    echo "=== Phase 0: Submitting input acquisition job ==="

    PHASE0_SCRIPT="${SCRIPTS_DIR}/phase0_test.sh"

    if [[ "$INPUT_MODE" == "bam" ]]; then
        # WGS: CRAM → HLA BAM extraction
        # WES: BAM → HLA BAM extraction (stream remote URL)
        if [[ "$TYPE" == "wgs" ]]; then
            PHASE0_BODY=$(cat << 'WGS_P0'
set -euo pipefail
module purge
module load samtools 2>/dev/null || true

BAM_OUT="BAM_DIR_PH/${SAMPLE}_hla.bam"
TMP_BAM="BAM_DIR_PH/${SAMPLE}_hla_tmp.bam"

if [[ -f "${BAM_OUT}" ]] && [[ -f "${BAM_OUT}.bai" ]]; then
    echo "[SKIP] HLA BAM already present"; exit 0
fi

LINE=$(grep -P "^${SAMPLE}\t" "URL_MAP_PH" || true)
[[ -z "$LINE" ]] && { echo "[WARN] No CRAM URL for ${SAMPLE} — skipping" >&2; exit 0; }
CRAM_URL=$(echo "$LINE" | cut -f3)
echo "[INFO] CRAM URL: ${CRAM_URL}"

export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="CRAM_CACHE_PH/%2s/%2s/%s"
mkdir -p "CRAM_CACHE_PH" "BAM_DIR_PH"

echo "[INFO] Streaming HLA region: HLA_REGION_PH"
if samtools view -b -h -o "${TMP_BAM}" "${CRAM_URL}" "HLA_REGION_PH" 2>/dev/null; then
    NREADS=$(samtools view -c "${TMP_BAM}" 2>/dev/null || echo 0)
    [[ "${NREADS}" -gt 0 ]] || { echo "[ERROR] 0 reads extracted" >&2; rm -f "${TMP_BAM}"; exit 1; }
    echo "[OK] HLA reads: ${NREADS}"
else
    echo "[ERROR] samtools failed for ${CRAM_URL}" >&2; rm -f "${TMP_BAM}"; exit 1
fi

samtools sort -@ 3 -m 2G -o "${BAM_OUT}" "${TMP_BAM}"
samtools index "${BAM_OUT}"
rm -f "${TMP_BAM}"
echo "[OK] HLA BAM: ${BAM_OUT} ($(samtools view -c "${BAM_OUT}") reads, indexed)"
WGS_P0
)
        else
            # WES
            PHASE0_BODY=$(cat << 'WES_P0'
set -euo pipefail
module purge
module load samtools 2>/dev/null || true

BAM_OUT="BAM_DIR_PH/${SAMPLE}_hla.bam"
TMP_BAM="BAM_DIR_PH/${SAMPLE}_hla_tmp.bam"

if [[ -f "${BAM_OUT}" ]] && [[ -f "${BAM_OUT}.bai" ]]; then
    echo "[SKIP] HLA BAM already present"; exit 0
fi

LINE=$(grep -P "^${SAMPLE}\t" "URL_MAP_PH" || true)
[[ -z "$LINE" ]] && { echo "[ERROR] No BAM URL for ${SAMPLE}" >&2; exit 1; }
PRIMARY_URL=$(echo "$LINE" | cut -f2)
ALT_URLS=$(echo "$LINE" | cut -f3 | tr '|' ' ')
mkdir -p "BAM_DIR_PH"

SUCCESS=false
for BAM_URL in ${PRIMARY_URL} ${ALT_URLS}; do
    echo "[INFO] Trying: ${BAM_URL}"
    if samtools view -b -h -o "${TMP_BAM}" "${BAM_URL}" "HLA_REGION_PH" 2>/dev/null; then
        NREADS=$(samtools view -c "${TMP_BAM}" 2>/dev/null || echo 0)
        if [[ "${NREADS}" -gt 0 ]]; then
            echo "[OK] HLA reads: ${NREADS}"
            SUCCESS=true
            break
        fi
    fi
    echo "[WARN] Failed or 0 reads for ${BAM_URL} — trying next"
    rm -f "${TMP_BAM}"
done

[[ "$SUCCESS" == "true" ]] || { echo "[ERROR] All URLs failed for ${SAMPLE}" >&2; exit 1; }

samtools sort -@ 3 -m 2G -o "${BAM_OUT}" "${TMP_BAM}"
samtools index "${BAM_OUT}"
rm -f "${TMP_BAM}"
echo "[OK] HLA BAM: ${BAM_OUT} ($(samtools view -c "${BAM_OUT}") reads, indexed)"
WES_P0
)
        fi

        cat > "$PHASE0_SCRIPT" << PHASE0EOF
#!/bin/bash
#SBATCH --job-name=${PHASE0_JOB_NAME}
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=${PHASE0_TIME}
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=${LOGS_DIR}/test_phase0.out
#SBATCH --error=${LOGS_DIR}/test_phase0.err

SAMPLE="${SAMPLE}"
echo "[INFO] Phase 0 — input acquisition for \${SAMPLE}  (${TYPE^^})"
echo "[INFO] Date: \$(date)"

${PHASE0_BODY}
PHASE0EOF

        sed -i \
            -e "s|BAM_DIR_PH|${BAM_DIR}|g" \
            -e "s|URL_MAP_PH|${URL_MAP}|g" \
            -e "s|HLA_REGION_PH|${HLA_REGION:-}|g" \
            -e "s|CRAM_CACHE_PH|${CRAM_REF_CACHE}|g" \
            "$PHASE0_SCRIPT"

    else
        # RNA: download FASTQ pair
        cat > "$PHASE0_SCRIPT" << PHASE0EOF
#!/bin/bash
#SBATCH --job-name=${PHASE0_JOB_NAME}
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=${PHASE0_TIME}
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/test_phase0.out
#SBATCH --error=${LOGS_DIR}/test_phase0.err

set -euo pipefail
SAMPLE="${SAMPLE}"
FASTQ_DIR="${FASTQ_DIR}"
URL_MAP="${URL_MAP}"
R1_OUT="\${FASTQ_DIR}/\${SAMPLE}_R1.fastq.gz"
R2_OUT="\${FASTQ_DIR}/\${SAMPLE}_R2.fastq.gz"

echo "[INFO] Phase 0 — FASTQ download for \${SAMPLE}  (RNA)"
echo "[INFO] Date: \$(date)"

if [[ -f "\${R1_OUT}" ]] && [[ -s "\${R1_OUT}" ]] && [[ -f "\${R2_OUT}" ]] && [[ -s "\${R2_OUT}" ]]; then
    echo "[SKIP] FASTQs already present for \${SAMPLE}"; exit 0
fi

LINE=\$(grep -P "^\${SAMPLE}\t" "\$URL_MAP" || true)
[[ -z "\$LINE" ]] && { echo "[ERROR] No FASTQ URL for \${SAMPLE}" >&2; exit 1; }
R1_URL=\$(echo "\$LINE" | cut -f2)
R2_URL=\$(echo "\$LINE" | cut -f3)
R1_ALTS=\$(echo "\$LINE" | cut -f4 | tr '|' ' ')
R2_ALTS=\$(echo "\$LINE" | cut -f5 | tr '|' ' ')

mkdir -p "\${FASTQ_DIR}"

download_file() {
    local URL="\$1" OUT="\$2"
    wget -q --tries=3 --timeout=120 -O "\${OUT}.tmp" "https://\${URL}" && mv "\${OUT}.tmp" "\${OUT}"
}

R1_OK=false
for URL in \${R1_URL} \${R1_ALTS}; do
    download_file "\${URL}" "\${R1_OUT}" && R1_OK=true && echo "[OK] R1: \${R1_OUT}" && break \
        || echo "[WARN] Failed: \${URL}"
    rm -f "\${R1_OUT}.tmp"
done

R2_OK=false
for URL in \${R2_URL} \${R2_ALTS}; do
    download_file "\${URL}" "\${R2_OUT}" && R2_OK=true && echo "[OK] R2: \${R2_OUT}" && break \
        || echo "[WARN] Failed: \${URL}"
    rm -f "\${R2_OUT}.tmp"
done

[[ "\$R1_OK" == "true" ]] && [[ "\$R2_OK" == "true" ]] \
    || { echo "[ERROR] FASTQ download failed for \${SAMPLE}" >&2; exit 1; }
echo "[OK] \${SAMPLE}: paired FASTQs downloaded"
PHASE0EOF
    fi   # end INPUT_MODE==fastq

    SBATCH_P0="sbatch --parsable --account=${PROJECT_ID} ${PHASE0_SCRIPT}"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_P0}"
        PHASE0_JOB="DRY_RUN_JOB"
    else
        PHASE0_JOB=$(eval "$SBATCH_P0")
        echo "[INFO] Phase 0 submitted: ${PHASE0_JOB}"
        PHASE0_DEP="--dependency=afterany:${PHASE0_JOB}"  # afterany: proceed even if extraction partially fails
    fi

else
    echo "[SKIP] Phase 0 — using existing input data"
fi

#-----------------------------------------------------------------------------
# Phase 1: Nextflow typing
#-----------------------------------------------------------------------------
PHASE1_DEP=""
TYPING_JOB=""

if [[ "$SKIP_TYPING" == "false" ]]; then
    echo ""
    echo "=== Phase 1: Submitting Nextflow typing job ==="

    # Samplesheet generation snippet (embedded in Phase 1 script)
    if [[ "$INPUT_MODE" == "bam" ]]; then
        GEN_SHEET="echo 'sample_id,bam_path' > \"${SAMPLESHEET}\"
BAM=\"${BAM_DIR}/${SAMPLE}_hla.bam\"
if [[ -f \"\$BAM\" ]] && [[ -f \"\${BAM}.bai\" ]]; then
    echo \"${SAMPLE},\${BAM}\" >> \"${SAMPLESHEET}\"
else
    echo \"[ERROR] HLA BAM not found: \${BAM}\" >&2; exit 1
fi"
    else
        GEN_SHEET="echo 'sample_id,fastq_1,fastq_2' > \"${SAMPLESHEET}\"
R1=\"${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz\"
R2=\"${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz\"
if [[ -f \"\$R1\" ]] && [[ -s \"\$R1\" ]] && [[ -f \"\$R2\" ]] && [[ -s \"\$R2\" ]]; then
    echo \"${SAMPLE},\${R1},\${R2}\" >> \"${SAMPLESHEET}\"
else
    echo \"[ERROR] FASTQs not found: \${R1} / \${R2}\" >&2; exit 1
fi"
    fi

    TYPING_SCRIPT="${SCRIPTS_DIR}/phase1_test_typing.sh"
    cat > "$TYPING_SCRIPT" << TYPINGEOF
#!/bin/bash
#SBATCH --job-name=test_${TYPE}_typi
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=24
#SBATCH --output=${LOGS_DIR}/test_typing.out
#SBATCH --error=${LOGS_DIR}/test_typing.err

set -euo pipefail
module purge
module load nextflow 2>/dev/null || module load nextflow/23.10.0
module load apptainer 2>/dev/null || true

echo "=== Phase 1: HLA Typing (${TYPE^^} / ${SAMPLE}) ==="
echo "Date: \$(date)"

cd "${INSTALL_DIR}"

# Build samplesheet
${GEN_SHEET}
echo "[OK] Samplesheet: ${SAMPLESHEET}"

nextflow run main.nf \\
    --input_samplesheet "${SAMPLESHEET}" \\
    --seq_type ${SEQ_TYPE} \\
    --reference ${REFERENCE} \\
    --tools "${TOOLS}" \\
    --weighting equal \\
    --skip_qc true \\
    --project "${PROJECT_ID}" \\
    --outdir "${RESULTS_DIR}" \\
    -c conf/puhti.config \\
    -profile apptainer \\
    -with-trace "${RESULTS_DIR}/pipeline_info/trace_test_${TYPE}.txt" \\
    -with-report "${RESULTS_DIR}/pipeline_info/report_test_${TYPE}.html" \\
    -work-dir "${SCRATCH_BASE}/work"

echo ""
echo "[OK] Phase 1 complete — results in ${RESULTS_DIR}"
TYPINGEOF

    SBATCH_TYPING="sbatch --parsable ${PHASE0_DEP} ${TYPING_SCRIPT}"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] ${SBATCH_TYPING}"
        TYPING_JOB="DRY_RUN_JOB"
    else
        TYPING_JOB=$(eval "$SBATCH_TYPING")
        echo "[INFO] Phase 1 typing job submitted: ${TYPING_JOB}"
        PHASE1_DEP="--dependency=afterany:${TYPING_JOB}"  # afterany: report even if some tools failed
    fi
else
    echo "[SKIP] Phase 1 — using existing typing results"
fi

#-----------------------------------------------------------------------------
# Phase 2: Collect results + generate test report
#-----------------------------------------------------------------------------
echo ""
echo "=== Phase 2: Submitting result collection + test report job ==="

REPORT_SCRIPT="${SCRIPTS_DIR}/phase2_report.sh"
cat > "$REPORT_SCRIPT" << REPORTEOF
#!/bin/bash
#SBATCH --job-name=test_${TYPE}_rpt
#SBATCH --account=${PROJECT_ID}
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=${LOGS_DIR}/test_report.out
#SBATCH --error=${LOGS_DIR}/test_report.err

set -euo pipefail
module purge
module load python-data 2>/dev/null || module load python/3.9 2>/dev/null || true

echo "=== Phase 2: Collecting results + generating test report ==="
echo "Sample:  ${SAMPLE}  (${TYPE^^})"
echo "Date:    \$(date)"
echo ""

# Collect results into by_tool layout
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    mkdir -p "${BY_TOOL_DIR}/\${TOOL}"
    SRC="${RESULTS_DIR}/${SAMPLE}/\${TOOL}/${SAMPLE}_\${TOOL}.txt"
    DST="${BY_TOOL_DIR}/\${TOOL}/${SAMPLE}_\${TOOL}.txt"
    if [[ -f "\$SRC" ]] && [[ ! -e "\$DST" ]]; then
        ln -sf "\$SRC" "\$DST"
        echo "[OK] Collected: \${TOOL}"
    elif [[ -f "\$SRC" ]]; then
        echo "[OK] Already collected: \${TOOL}"
    else
        echo "[WARN] No result for \${TOOL}: \${SRC}"
    fi
done

echo ""
echo "=== Test Report ==="

# Generate report using per-tool inspect output filtered to this sample
TMP_REPORT_DIR=\$(mktemp -d)
{
    if [[ -f "${GT_FILE}" ]]; then
        echo "[INFO] GT-backed concordance for ${SAMPLE}"
        for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
            SRC="${RESULTS_DIR}/${SAMPLE}/\${TOOL}/${SAMPLE}_\${TOOL}.txt"
            echo ""
            echo "--- \${TOOL} ---"
            if [[ ! -f "\$SRC" ]]; then
                echo "NO OUTPUT"
                continue
            fi

            TOOL_TSV="\${TMP_REPORT_DIR}/\${TOOL}.tsv"
            python3 "${INSTALL_DIR}/bin/calibrate_tool_weights.py" inspect \\
                --ground-truth "${GT_FILE}" \\
                --results-dir  "${BY_TOOL_DIR}" \\
                --tool         "\${TOOL}" \\
                --genes        A,B,C,DRB1,DQB1 \\
                --output       "\${TOOL_TSV}" >/dev/null

            awk -F '\t' -v sample="${SAMPLE}" 'NR==1 || \$1==sample' "\${TOOL_TSV}" || true
        done
    else
        echo "[WARN] GT file not found: ${GT_FILE}"
        echo "[INFO] Showing raw tool outputs:"
        for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
            SRC="${RESULTS_DIR}/${SAMPLE}/\${TOOL}/${SAMPLE}_\${TOOL}.txt"
            if [[ -f "\$SRC" ]]; then
                echo ""
                echo "--- \${TOOL} ---"
                cat "\$SRC"
            else
                echo "--- \${TOOL}: NO OUTPUT ---"
            fi
        done
    fi
} | tee "${REPORT_FILE}"
rm -rf "\${TMP_REPORT_DIR}"

echo ""
echo "=== Summary ==="
N_PASS=0; N_FAIL=0
TOTAL_TOOLS=\$(echo "${TOOLS}" | tr ',' '\n' | wc -l)
for TOOL in \$(echo "${TOOLS}" | tr ',' ' '); do
    SRC="${RESULTS_DIR}/${SAMPLE}/\${TOOL}/${SAMPLE}_\${TOOL}.txt"
    if [[ -f "\$SRC" ]]; then
        echo "  [PASS] \${TOOL}"
        N_PASS=\$((N_PASS+1))
    else
        echo "  [FAIL] \${TOOL} — no output"
        N_FAIL=\$((N_FAIL+1))
    fi
done
echo ""
echo "Result: \${N_PASS}/\${TOTAL_TOOLS} tools produced output"
if [[ "\${N_PASS}" -eq 0 ]]; then
    echo "STATUS: FAIL (0 tools produced output)"
    echo "Reason: Phase 1 completed, but every requested tool failed or was ignored."
    echo "Action: inspect ${LOGS_DIR}/test_typing.out and ${LOGS_DIR}/test_typing.err for tool-level failures."
    echo ""
    echo "Full report saved: ${REPORT_FILE}"
    echo "SLURM logs:        ${LOGS_DIR}/"
    exit 2
elif [[ "\${N_FAIL}" -eq 0 ]]; then
    echo "STATUS: PASS"
else
    echo "STATUS: WARN (\${N_FAIL} tools missing)"
fi
echo ""
echo "Full report saved: ${REPORT_FILE}"
echo "SLURM logs:        ${LOGS_DIR}/"
REPORTEOF

SBATCH_REPORT="sbatch --parsable ${PHASE1_DEP} ${REPORT_SCRIPT}"
if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] ${SBATCH_REPORT}"
else
    REPORT_JOB=$(eval "$SBATCH_REPORT")
    echo "[INFO] Phase 2 report job submitted: ${REPORT_JOB}"
fi

echo ""
echo "==================================================================="
echo " Jobs submitted. Monitor with:"
echo "   squeue -u \$USER"
echo ""
echo " Check status any time with:"
echo "   bash scripts/run_e2e_test_puhti.sh --type ${TYPE} --project ${PROJECT_ID} --sample ${SAMPLE} --status"
echo ""
echo " After completion, view report:"
echo "   cat ${REPORT_FILE}"
echo "==================================================================="
