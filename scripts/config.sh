#!/bin/bash
# Configuration for HLA typing pipeline - BAM analysis

# CSC Puhti settings
export PROJECT_ID="project_2008084"
export BASE_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"

# Sample information
export SAMPLE_ID="NA12273"
export BAM_FILE="NA12273.bam"
export BAI_FILE="NA12273.bam.bai"

# Pipeline settings
export TOOLS="optitype,arcashla,spechla"
export ENABLE_CONSENSUS="true"
export SEQ_TYPE="dna"  # Important: WGS DNA data!

# Resource settings
export MAX_CPUS=40
export MAX_MEMORY="180.GB"
export MAX_TIME="24.h"

# Paths
export RAW_BAM_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/raw_bam"
export PIPELINE_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
export INPUT_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/pipeline_input"
export RESULTS_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/results"
export LOGS_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/logs"
export SINGULARITY_CACHE="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/singularity_cache"

# Reference data
export HLA_REFERENCES="/scratch/project_2008084/hla_references"

# SLURM settings
export SLURM_ACCOUNT="project_2008084"
export SLURM_PARTITION="small"
