#!/bin/bash
# Build basetools container for HLA typing pipeline
# Contains: samtools, bc, fastqc, python3, multiqc

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_DIR="${1:-${SCRIPT_DIR}/../../hla_references/containers}"

# Create container directory if it doesn't exist
mkdir -p "${CONTAINER_DIR}"

# Create Singularity definition file
cat > "${CONTAINER_DIR}/basetools.def" << 'EOF'
Bootstrap: docker
From: ubuntu:22.04

%labels
    Author HLA Pipeline Team
    Description Base tools container for HLA typing pipeline
    Version 1.0

%post
    apt-get update && apt-get install -y --no-install-recommends \
        bc \
        gawk \
        wget \
        curl \
        unzip \
        bzip2 \
        ca-certificates \
        default-jre \
        python3 \
        python3-pip \
        python3-dev \
        samtools \
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

    # Install FastQC
    cd /opt
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
    unzip fastqc_v0.12.1.zip
    chmod +x FastQC/fastqc
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc
    rm fastqc_v0.12.1.zip

    # Install Python packages
    pip3 install --no-cache-dir \
        numpy \
        pandas \
        matplotlib \
        seaborn \
        multiqc \
        scipy

%environment
    export LC_ALL=C
    export PATH=/opt/FastQC:$PATH

%runscript
    exec "$@"

%test
    samtools --version
    fastqc --version
    python3 --version
    multiqc --version
EOF

echo "Building basetools container..."
cd "${CONTAINER_DIR}"

# Build the container
if command -v singularity &> /dev/null; then
    singularity build basetools.sif basetools.def
    echo "Container built successfully: ${CONTAINER_DIR}/basetools.sif"
else
    echo "Singularity is not installed. Definition file created at: ${CONTAINER_DIR}/basetools.def"
    echo "Build manually with: singularity build basetools.sif basetools.def"
fi
