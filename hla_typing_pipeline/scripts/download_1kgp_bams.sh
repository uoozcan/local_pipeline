#!/bin/bash
#SBATCH --job-name=1kgp_download
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/project_2008084/ozcanumu/hla_calibration/logs/download_%j.out
#SBATCH --error=/scratch/project_2008084/ozcanumu/hla_calibration/logs/download_%j.err
#
# Download 1KGP Phase 3 low-coverage BAMs from EBI for HLA calibration
#
# Submit:
#   sbatch scripts/download_1kgp_bams.sh
#
# Uses lftp mget with lcd — correctly handles FTP glob patterns.
# Downloads to a per-sample temp dir then renames to {SAMPLE}.bam.
# Approximate size: ~3-5 GB per sample × 50 samples = ~150-250 GB
#=============================================================================

set -euo pipefail

BAM_DIR="/scratch/project_2008084/ozcanumu/hla_calibration/1kgp_bams"
mkdir -p "$BAM_DIR"

# EBI 1KGP Phase 3 alignment base path (on the FTP server)
EBI_HOST="ftp://ftp.1000genomes.ebi.ac.uk"

echo "=== 1KGP BAM Download ==="
echo "Target directory: $BAM_DIR"
echo "Start: $(date)"
echo ""

# Verify lftp is available (available by default on Puhti login/compute nodes)
if ! command -v lftp &>/dev/null; then
    echo "ERROR: lftp not found. Try: module load lftp"
    exit 1
fi

# Population codes for each sample — required for exact EBI filename
declare -A SAMPLE_POP=(
    [NA12878]=CEU [NA12891]=CEU [NA12892]=CEU
    [NA12877]=CEU [NA12889]=CEU [NA12890]=CEU
    [NA12873]=CEU [NA12874]=CEU [NA12875]=CEU
    [NA19240]=YRI [NA19238]=YRI [NA19239]=YRI
    [NA19209]=YRI [NA19210]=YRI [NA19129]=YRI
    [NA19130]=YRI [NA19131]=YRI [NA19152]=YRI
    [NA18526]=CHB [NA18524]=CHB [NA18529]=CHB
    [NA18532]=CHB [NA18537]=CHB [NA18542]=CHB
    [NA18561]=CHB [NA18562]=CHB [NA18563]=CHB
    [HG00096]=GBR [HG00097]=GBR [HG00099]=GBR
    [HG00100]=GBR [HG00101]=GBR [HG00102]=GBR
    [HG00103]=GBR [HG00105]=GBR [HG00106]=GBR
    [HG00107]=GBR [HG00108]=GBR [HG00109]=GBR
    [HG00110]=GBR [NA20502]=TSI [NA20503]=TSI
    [NA20504]=TSI [NA20505]=TSI [NA20506]=TSI
    [NA20507]=TSI [NA20508]=TSI [NA20509]=TSI
    [NA20510]=TSI [NA20511]=TSI
)

OK=0
SKIP=0
FAIL=0

for SAMPLE in "${!SAMPLE_POP[@]}"; do
    BAM_OUT="${BAM_DIR}/${SAMPLE}.bam"
    BAI_OUT="${BAM_DIR}/${SAMPLE}.bam.bai"

    if [[ -f "$BAM_OUT" ]] && [[ -f "$BAI_OUT" ]]; then
        echo "[SKIP] $SAMPLE — already exists"
        (( SKIP++ )) || true
        continue
    fi

    POP="${SAMPLE_POP[$SAMPLE]}"
    REMOTE_DIR="/vol1/ftp/phase3/data/${SAMPLE}/alignment"
    GLOB_BAM="${SAMPLE}.mapped.ILLUMINA.bwa.${POP}.low_coverage.*.bam"
    GLOB_BAI="${GLOB_BAM}.bai"

    echo "[DL]   $SAMPLE ($POP)"

    # Use lftp mget with lcd:
    #   - lcd sets the local download directory
    #   - cd sets the remote directory
    #   - mget expands the glob server-side and downloads matching files
    # Files land in BAM_DIR with their original names (e.g. NA12878.mapped...bam)
    TMPDIR="${BAM_DIR}/.tmp_${SAMPLE}"
    mkdir -p "$TMPDIR"

    lftp "${EBI_HOST}" << LFTPEOF 2>&1 || { echo "[WARN] lftp failed for $SAMPLE"; (( FAIL++ )) || true; rm -rf "$TMPDIR"; continue; }
set net:max-retries 5
set net:reconnect-interval-base 10
set ftp:passive-mode true
lcd ${TMPDIR}
cd ${REMOTE_DIR}
mget -c ${GLOB_BAM}
mget -c ${GLOB_BAI}
quit
LFTPEOF

    # Find the downloaded BAM (original name includes pop code and date)
    DOWNLOADED_BAM=$(ls "${TMPDIR}/${SAMPLE}.mapped."*.bam 2>/dev/null | grep -v '\.bai$' | head -1 || true)
    DOWNLOADED_BAI=$(ls "${TMPDIR}/${SAMPLE}.mapped."*.bam.bai 2>/dev/null | head -1 || true)

    if [[ -z "$DOWNLOADED_BAM" ]]; then
        echo "[WARN] $SAMPLE: BAM not found after download (check logs)"
        (( FAIL++ )) || true
        rm -rf "$TMPDIR"
        continue
    fi

    # Rename to standard {SAMPLE}.bam for pipeline
    mv "$DOWNLOADED_BAM" "$BAM_OUT"
    echo "       $(basename "$DOWNLOADED_BAM") -> $(basename "$BAM_OUT")"

    if [[ -n "$DOWNLOADED_BAI" ]]; then
        mv "$DOWNLOADED_BAI" "$BAI_OUT"
    else
        # Index BAM if .bai was not downloaded or not available
        echo "       Indexing (no .bai downloaded)..."
        module load samtools 2>/dev/null || module load biokit 2>/dev/null || true
        samtools index "$BAM_OUT"
    fi

    rm -rf "$TMPDIR"
    echo "       [OK] $SAMPLE  ($(du -sh "$BAM_OUT" | cut -f1))"
    (( OK++ )) || true
done

echo ""
echo "=== Download complete ==="
echo "  OK:      $OK"
echo "  Skipped: $SKIP"
echo "  Failed:  $FAIL"
echo "  End: $(date)"
echo ""
BAM_COUNT=$(ls "$BAM_DIR"/*.bam 2>/dev/null | wc -l || echo 0)
echo "Total BAMs in $BAM_DIR: $BAM_COUNT"
echo ""
if [[ $FAIL -gt 0 ]]; then
    echo "NOTE: $FAIL samples failed. Re-run this script to retry (already-downloaded files are skipped)."
    echo ""
fi
echo "Next step — submit typing + calibration:"
echo "  cd /projappl/project_2008084/hla_typing/pipeline"
echo "  bash scripts/submit_calibration_puhti.sh --project project_2008084"
