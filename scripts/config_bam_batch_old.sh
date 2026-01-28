#!/bin/bash
# Configuration for Batch HLA Typing Pipeline - BAM Analysis
# VenEx DNA-seq Batch Processing

# ==============================================================================
# CSC Puhti Settings
# ==============================================================================
export PROJECT_ID="project_2008084"
export USER_NAME="ozcanumu"
export BASE_DIR="/scratch/${PROJECT_ID}/${USER_NAME}/hla_rnaseq_analysis"

# ==============================================================================
# Batch Information
# ==============================================================================
export BATCH_NAME="batch1_VenEx_DNA_BAM"
export DATA_TYPE="DNA"  # WGS DNA data

# ==============================================================================
# Paths
# ==============================================================================
export RAW_BAM_DIR="${BASE_DIR}/raw_bam/batch1_VenEx_DNA"
export PIPELINE_DIR="${BASE_DIR}/hla_typing_pipeline"
export INPUT_DIR="${BASE_DIR}/pipeline_input_bam/${BATCH_NAME}"
export RESULTS_DIR="${BASE_DIR}/results_bam/${BATCH_NAME}"
export LOGS_DIR="${BASE_DIR}/logs"
export SINGULARITY_CACHE="${BASE_DIR}/singularity_cache"

# ==============================================================================
# File Lists
# ==============================================================================
export BAM_FILE_LIST="${BASE_DIR}/batch1_bam_file_list.txt"
export SAMPLE_SHEET="${INPUT_DIR}/samplesheet.csv"

# ==============================================================================
# Pipeline Settings
# ==============================================================================
# Tools optimized for WGS BAM analysis
export TOOLS="optitype,arcashla,spechla"  # ArcasHLA, xHLA, HLA-HD work well with BAM
export SEQ_TYPE="dna"  # DNA sequencing data

# HLA genes to analyze
export HLA_GENES="A,B,C,DRB1,DQB1,DPB1"

# Majority voting settings
export ENABLE_CONSENSUS="true"
export MIN_TOOLS_CONSENSUS=1
export CONSENSUS_RESOLUTION="2field"

# ==============================================================================
# Resource Settings
# ==============================================================================
export MAX_CPUS=40
export MAX_MEMORY="180.GB"
export MAX_TIME="24.h"

# SLURM settings
export SLURM_ACCOUNT="${PROJECT_ID}"
export SLURM_PARTITION="small"

# ==============================================================================
# Reference Data
# ==============================================================================
export HLA_REFERENCES="/scratch/${PROJECT_ID}/hla_references"
export REFERENCE_GENOME="${HLA_REFERENCES}/hs38DH.fa"  # HLA-aware reference

# ==============================================================================
# Helper Functions
# ==============================================================================

# Extract sample ID from BAM filename
# Input: VX_36_4_D1_tumor.bam -> Output: VX_36_4_D1
get_sample_id() {
    local bam_file="$1"
    echo "$(basename ${bam_file} | sed 's/_tumor\.bam\.c4gh$//')"
}

# Get subdirectory path for sample
# Input: VX_36_4_D1 -> Output: VenEx_DNA_NONHUS_Batch1/VX_36_4_D1
get_sample_subdir() {
    local sample_id="$1"
    echo "VenEx_DNA_NONHUS_Batch1/${sample_id}"
}