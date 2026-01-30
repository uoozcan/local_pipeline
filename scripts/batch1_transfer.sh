#!/bin/bash
#SBATCH --job-name=batch1_transfer_v2
#SBATCH --account=project_2008084
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=small
#SBATCH --output=logs/batch1_transfer_v2_%j.log
#SBATCH --error=logs/batch1_transfer_v2_%j.err

# Batch 1 Transfer Script - DNA-seq VenEx samples (ROBUST)
# Continue on errors, retry failed transfers

set -uo pipefail  # Removed -e to continue on errors

export C4GH_PASSPHRASE="kontrogroup2025"

# Configuration
BATCH_NAME="batch1_VenEx_DNA"
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
FILE_LIST="${WORK_DIR}/batch1_file_list.txt"
DEST_PATH="${WORK_DIR}/raw_fastq/${BATCH_NAME}"
SECRET_KEY="${WORK_DIR}/all_data_csc_key.sec"
LOG_DIR="${WORK_DIR}/logs"
TRANSFER_LOG="${LOG_DIR}/transfer_v2_${BATCH_NAME}_$(date +%Y%m%d_%H%M%S).log"
FAILED_LIST="${LOG_DIR}/failed_transfers_$(date +%Y%m%d_%H%M%S).txt"

# Create directories
mkdir -p "${DEST_PATH}"
mkdir -p "${LOG_DIR}"

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
log_msg "Batch 1 Transfer Started (ROBUST v2)"
log_msg "=========================================="
log_msg "Batch: ${BATCH_NAME}"
log_msg "File list: ${FILE_LIST}"
log_msg "Destination: ${DEST_PATH}"
log_msg "Secret key: ${SECRET_KEY}"
log_msg ""

# Verify file list exists
if [ ! -f "${FILE_LIST}" ]; then
    log_msg "ERROR: File list not found: ${FILE_LIST}"
    exit 1
fi

# Verify secret key exists
if [ ! -f "${SECRET_KEY}" ]; then
    log_msg "ERROR: Secret key not found: ${SECRET_KEY}"
    exit 1
fi

# Count total files
total_files=$(wc -l < "${FILE_LIST}")
log_msg "Total files to transfer: ${total_files}"
log_msg ""

# Change to destination directory
cd "${DEST_PATH}"

# Process each file
while IFS= read -r file_path; do
    # Skip empty lines
    [ -z "$file_path" ] && continue
    
    # Extract file information
    file_name=$(basename "$file_path")
    decrypted_name="${file_name%.c4gh}"
    
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
            sleep 10  # Wait before retry
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
    
done < "${FILE_LIST}"

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
    log_msg "To retry failed transfers, run:"
    log_msg "  sbatch --export=FILE_LIST=${FAILED_LIST} batch1_transfer_v2.sh"
else
    # Create transfer completion marker
    touch "${DEST_PATH}/.transfer_complete"
    log_msg "All transfers completed successfully!"
    log_msg "Marker file created: ${DEST_PATH}/.transfer_complete"
fi

exit 0  # Always exit with 0 to allow reruns