#!/bin/bash
#SBATCH --job-name=transfer_bam
#SBATCH --account=project_2008084
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=small
#SBATCH --output=logs/transfer_bam_%j.log
#SBATCH --error=logs/transfer_bam_%j.err

# STEP 1: Transfer Encrypted BAM Files from Allas to Puhti
# Handles Crypt4GH decryption automatically

set -uo pipefail  # Continue on errors

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/config_batch_bam.sh

# Set passphrase
export C4GH_PASSPHRASE

# Configuration
BAM_FILE_LIST="${WORK_DIR}/batch1_file_list.txt"
DEST_PATH="${RAW_BAM_DIR}"
TRANSFER_LOG="${LOGS_DIR}/transfer_bam_$(date +%Y%m%d_%H%M%S).log"
FAILED_LIST="${LOGS_DIR}/failed_bam_transfers_$(date +%Y%m%d_%H%M%S).txt"

# Create directories
mkdir -p "${DEST_PATH}"
mkdir -p "${LOGS_DIR}"

# Initialize counters
total_files=0
successful_transfers=0
failed_transfers=0
skipped_files=0

# Logging function
log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${TRANSFER_LOG}"
}

log_msg "=========================================="
log_msg "BAM Files Transfer Started"
log_msg "=========================================="
log_msg "Batch: ${BATCH_NAME}"
log_msg "File list: ${BAM_FILE_LIST}"
log_msg "Destination: ${DEST_PATH}"
log_msg "Secret key: ${SECRET_KEY}"
log_msg ""

# Verify file list exists
if [ ! -f "${BAM_FILE_LIST}" ]; then
    log_msg "ERROR: File list not found: ${BAM_FILE_LIST}"
    exit 1
fi

# Count total files
total_files=$(wc -l < "${BAM_FILE_LIST}")
log_msg "Total BAM files to transfer: ${total_files}"
log_msg ""

# Change to destination directory
cd "${DEST_PATH}"

# Process each file
while IFS= read -r file_path; do
    # Skip empty lines
    [ -z "$file_path" ] && continue
    
    # Extract file information
    file_name=$(basename "$file_path")
    decrypted_name="${file_name%.c4gh}"  # Remove .c4gh extension
    
    # Build the relative path for organization
    relative_path="${file_path#2014061-DNA_seq/DNA_seq/DRAGEN/}"
    subdir=$(dirname "$relative_path")
    
    # Expected final location
    output_file="${relative_path%.c4gh}"
    
    # Create subdirectory if needed
    if [ "$subdir" != "." ]; then
        mkdir -p "$subdir"
    fi
    
    # Check if file already exists at final location
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        log_msg "SKIP: ${output_file} already exists"
        ((successful_transfers++))
        ((skipped_files++))
        continue
    fi
    
    # Transfer with decryption
    log_msg "Transferring: ${file_path}"
    
    # Retry logic: try up to 3 times
    transfer_success=false
    for attempt in 1 2 3; do
        if [ $attempt -gt 1 ]; then
            log_msg "  Retry attempt ${attempt}/3..."
            sleep 30  # Wait before retry (BAM files are large)
        fi
        
        # a-get downloads to current directory with just the filename
        if a-get "$file_path" --sk "${SECRET_KEY}" 2>> "${TRANSFER_LOG}"; then
            # File should now be in current directory as decrypted_name
            if [ -f "$decrypted_name" ] && [ -s "$decrypted_name" ]; then
                # Move to proper subdirectory
                if mv "$decrypted_name" "$output_file" 2>> "${TRANSFER_LOG}"; then
                    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
                        file_size=$(du -h "$output_file" | cut -f1)
                        log_msg "SUCCESS: ${output_file} (${file_size})"
                        
                        # Create BAM index
                        log_msg "  Creating BAM index..."
                        if module load samtools 2>/dev/null; then
                            if samtools index "$output_file" 2>> "${TRANSFER_LOG}"; then
                                log_msg "  Index created: ${output_file}.bai"
                            else
                                log_msg "  WARNING: Index creation failed (will retry later)"
                            fi
                        fi
                        
                        ((successful_transfers++))
                        transfer_success=true
                        break
                    fi
                fi
            fi
        fi
        
        # Clean up partial downloads
        if [ -f "$decrypted_name" ]; then
            rm -f "$decrypted_name"
        fi
    done
    
    # If all retries failed, log it
    if [ "$transfer_success" = false ]; then
        log_msg "ERROR: Failed to transfer ${file_path} after 3 attempts"
        echo "$file_path" >> "${FAILED_LIST}"
        ((failed_transfers++))
    fi
    
    log_msg ""
    
done < "${BAM_FILE_LIST}"

# Clear passphrase
unset C4GH_PASSPHRASE

# Summary
log_msg "=========================================="
log_msg "Transfer Summary"
log_msg "=========================================="
log_msg "Completed: $(date)"
log_msg "Total files: ${total_files}"
log_msg "Successful: ${successful_transfers}"
log_msg "Skipped (already present): ${skipped_files}"
log_msg "Failed: ${failed_transfers}"
log_msg "=========================================="

if [ $failed_transfers -gt 0 ]; then
    log_msg "Failed files list saved to: ${FAILED_LIST}"
    log_msg ""
    log_msg "To retry failed transfers:"
    log_msg "  Edit batch1_bam_file_list.txt with failed files"
    log_msg "  sbatch step1_transfer_bam.sh"
else
    # Create transfer completion marker
    touch "${DEST_PATH}/.transfer_complete"
    log_msg "All transfers completed successfully!"
    log_msg "Marker file created: ${DEST_PATH}/.transfer_complete"
fi

log_msg ""
log_msg "NEXT STEP:"
log_msg "  sbatch step2_prepare_batch_bam.sh"

exit 0