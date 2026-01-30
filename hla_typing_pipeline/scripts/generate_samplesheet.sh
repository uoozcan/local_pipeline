#!/bin/bash
#
# Generate samplesheet from input files
#
# Usage:
#   ./generate_samplesheet.sh /path/to/bam_dir > samples_bam.csv
#   ./generate_samplesheet.sh /path/to/fastq_dir --fastq > samples_fastq.csv
#

set -euo pipefail

show_help() {
    cat << EOF
Generate Samplesheet for HLA Typing Pipeline
=============================================

Usage: $(basename "$0") <input_dir> [options]

Options:
  --fastq           Generate FASTQ samplesheet (default: BAM)
  --pattern <pat>   File pattern for sample ID extraction
  --output <file>   Write to file instead of stdout
  --help            Show this help

Examples:
  # BAM files
  $(basename "$0") /path/to/bam_files > samples_bam.csv

  # FASTQ files
  $(basename "$0") /path/to/fastq_files --fastq > samples_fastq.csv

  # Custom pattern
  $(basename "$0") /path/to/files --pattern '_sorted.bam' > samples.csv

EOF
    exit 0
}

# Defaults
INPUT_DIR=""
INPUT_TYPE="bam"
PATTERN=""
OUTPUT=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fastq)
            INPUT_TYPE="fastq"
            shift
            ;;
        --pattern)
            PATTERN="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            ;;
        -*)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
        *)
            INPUT_DIR="$1"
            shift
            ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory required" >&2
    echo "Use --help for usage" >&2
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Directory not found: $INPUT_DIR" >&2
    exit 1
fi

# Make path absolute
INPUT_DIR=$(cd "$INPUT_DIR" && pwd)

# Generate samplesheet
generate_bam_samplesheet() {
    echo "sample_id,bam_path"

    for bam in "$INPUT_DIR"/*.bam; do
        [[ -f "$bam" ]] || continue

        # Extract sample ID
        filename=$(basename "$bam")
        if [[ -n "$PATTERN" ]]; then
            sample_id="${filename%$PATTERN}"
        else
            sample_id="${filename%.bam}"
            # Remove common suffixes
            sample_id="${sample_id%.sorted}"
            sample_id="${sample_id%.dedup}"
            sample_id="${sample_id%.aligned}"
        fi

        echo "${sample_id},${bam}"
    done
}

generate_fastq_samplesheet() {
    echo "sample_id,fastq_1,fastq_2"

    # Find R1 files and match with R2
    for r1 in "$INPUT_DIR"/*_R1*.fastq.gz "$INPUT_DIR"/*_1.fastq.gz "$INPUT_DIR"/*_R1*.fq.gz "$INPUT_DIR"/*_1.fq.gz; do
        [[ -f "$r1" ]] || continue

        filename=$(basename "$r1")

        # Determine R2 file
        if [[ "$filename" == *"_R1"* ]]; then
            r2="${r1/_R1/_R2}"
        elif [[ "$filename" == *"_1.fastq"* ]]; then
            r2="${r1/_1.fastq/_2.fastq}"
        elif [[ "$filename" == *"_1.fq"* ]]; then
            r2="${r1/_1.fq/_2.fq}"
        else
            continue
        fi

        if [[ ! -f "$r2" ]]; then
            echo "WARNING: No R2 found for $r1" >&2
            continue
        fi

        # Extract sample ID
        sample_id=$(basename "$r1")
        sample_id="${sample_id%_R1*}"
        sample_id="${sample_id%_1.fastq*}"
        sample_id="${sample_id%_1.fq*}"

        echo "${sample_id},${r1},${r2}"
    done
}

# Output
if [[ "$INPUT_TYPE" == "bam" ]]; then
    if [[ -n "$OUTPUT" ]]; then
        generate_bam_samplesheet > "$OUTPUT"
        echo "Generated BAM samplesheet: $OUTPUT" >&2
        echo "Samples: $(tail -n +2 "$OUTPUT" | wc -l)" >&2
    else
        generate_bam_samplesheet
    fi
else
    if [[ -n "$OUTPUT" ]]; then
        generate_fastq_samplesheet > "$OUTPUT"
        echo "Generated FASTQ samplesheet: $OUTPUT" >&2
        echo "Samples: $(tail -n +2 "$OUTPUT" | wc -l)" >&2
    else
        generate_fastq_samplesheet
    fi
fi
