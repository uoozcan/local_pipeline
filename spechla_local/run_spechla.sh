#!/bin/bash
# Wrapper script to run SpecHLA with proper PATH settings
# Usage: ./run_spechla.sh [SpecHLA options]
# Example: ./run_spechla.sh -n SAMPLE -b /path/to/sample.bam -o ./output -j 8 -u 1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set PATH to include local tools and miniconda
export PATH="${SCRIPT_DIR}/local_bin:/home/umut/miniconda3/bin:$PATH"

# Run SpecHLA
bash "${SCRIPT_DIR}/script/whole/SpecHLA.sh" "$@"
