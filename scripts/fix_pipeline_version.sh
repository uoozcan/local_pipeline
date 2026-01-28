#!/bin/bash
# Check and Fix Pipeline Version

echo "=========================================="
echo "Pipeline Version Check & Fix"
echo "=========================================="
echo ""

# Load configuration
WORK_DIR="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis"
source ${WORK_DIR}/scripts/config_bam_batch.sh

cd "${PIPELINE_DIR}"

echo "Pipeline directory: ${PIPELINE_DIR}"
echo ""

# Check if it's a Git repository
if [ ! -d ".git" ]; then
    echo "✗ Not a Git repository"
    echo ""
    echo "This pipeline needs to be a Git repository for version management."
    echo "Options:"
    echo "  1. Clone a fresh version of the pipeline"
    echo "  2. Initialize Git: git init && git add . && git commit -m 'Initial'"
    exit 1
fi

echo "✓ Git repository detected"
echo ""

# Check current status
echo "Current Status:"
echo "---------------"
git status --short
echo ""

# Get current branch/tag
CURRENT_BRANCH=$(git branch --show-current 2>/dev/null)
CURRENT_COMMIT=$(git rev-parse --short HEAD 2>/dev/null)

echo "Current branch: ${CURRENT_BRANCH}"
echo "Current commit: ${CURRENT_COMMIT}"
echo ""

# Check for available tags
echo "=========================================="
echo "Available Versions (Tags)"
echo "=========================================="
echo ""

TAGS=$(git tag -l | sort -V)

if [ -z "${TAGS}" ]; then
    echo "No release tags found"
    echo ""
    echo "Checking for remote tags..."
    git fetch --tags 2>/dev/null || true
    TAGS=$(git tag -l | sort -V)
    
    if [ -z "${TAGS}" ]; then
        echo "Still no tags available"
        echo "You may be using a development version"
    fi
fi

if [ ! -z "${TAGS}" ]; then
    echo "Available release versions:"
    echo "${TAGS}" | tail -10
    echo ""
    
    LATEST_TAG=$(echo "${TAGS}" | tail -1)
    echo "Latest release: ${LATEST_TAG}"
    echo ""
fi

# Check for uncommitted changes
echo "=========================================="
echo "Repository State"
echo "=========================================="
echo ""

if [ -n "$(git status --porcelain)" ]; then
    echo "⚠ WARNING: You have uncommitted changes"
    echo ""
    git status --short
    echo ""
    echo "These changes might be causing issues."
    echo "Consider stashing them before switching versions."
else
    echo "✓ No uncommitted changes"
fi

# Offer fix options
echo ""
echo "=========================================="
echo "Fix Options"
echo "=========================================="
echo ""

echo "Choose an action:"
echo "  1) Switch to latest stable release (${LATEST_TAG})"
echo "  2) Switch to main/master branch"
echo "  3) Pull latest changes from current branch"
echo "  4) Show detailed Git log"
echo "  5) Check for specific version"
echo "  6) Cancel (no changes)"
echo ""

read -p "Select option (1-6): " choice

case $choice in
    1)
        if [ -z "${LATEST_TAG}" ]; then
            echo "No tags available to switch to"
            exit 1
        fi
        
        echo ""
        echo "Switching to ${LATEST_TAG}..."
        
        # Stash changes if any
        if [ -n "$(git status --porcelain)" ]; then
            echo "Stashing uncommitted changes..."
            git stash
        fi
        
        # Checkout tag
        if git checkout "${LATEST_TAG}"; then
            echo ""
            echo "✓ Switched to ${LATEST_TAG}"
            echo ""
            echo "Pipeline is now at stable release ${LATEST_TAG}"
            echo ""
            echo "Try running the pipeline again:"
            echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
        else
            echo "✗ Failed to checkout ${LATEST_TAG}"
            exit 1
        fi
        ;;
        
    2)
        echo ""
        echo "Attempting to switch to main/master branch..."
        
        # Try main first, then master
        if git rev-parse --verify main >/dev/null 2>&1; then
            TARGET_BRANCH="main"
        elif git rev-parse --verify master >/dev/null 2>&1; then
            TARGET_BRANCH="master"
        else
            echo "✗ Neither 'main' nor 'master' branch found"
            exit 1
        fi
        
        echo "Switching to ${TARGET_BRANCH}..."
        
        # Stash changes if any
        if [ -n "$(git status --porcelain)" ]; then
            git stash
        fi
        
        if git checkout "${TARGET_BRANCH}"; then
            git pull origin "${TARGET_BRANCH}" 2>/dev/null || true
            echo ""
            echo "✓ Switched to ${TARGET_BRANCH} and pulled latest"
            echo ""
            echo "Try running the pipeline again:"
            echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
        else
            echo "✗ Failed to checkout ${TARGET_BRANCH}"
            exit 1
        fi
        ;;
        
    3)
        echo ""
        echo "Pulling latest changes for ${CURRENT_BRANCH}..."
        
        if git pull origin "${CURRENT_BRANCH}" 2>/dev/null; then
            echo ""
            echo "✓ Updated to latest ${CURRENT_BRANCH}"
            echo ""
            echo "Try running the pipeline again:"
            echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
        else
            echo "✗ Failed to pull updates"
            echo "  You may need to check remote configuration"
        fi
        ;;
        
    4)
        echo ""
        echo "Recent commits:"
        git log --oneline --graph --decorate -20
        ;;
        
    5)
        echo ""
        read -p "Enter tag/branch/commit to checkout: " TARGET
        
        if [ -z "${TARGET}" ]; then
            echo "No target specified"
            exit 1
        fi
        
        echo "Switching to ${TARGET}..."
        
        if [ -n "$(git status --porcelain)" ]; then
            git stash
        fi
        
        if git checkout "${TARGET}"; then
            echo ""
            echo "✓ Switched to ${TARGET}"
            echo ""
            echo "Try running the pipeline again:"
            echo "  sbatch ${WORK_DIR}/scripts/step3_run_pipeline_bam.sh"
        else
            echo "✗ Failed to checkout ${TARGET}"
            exit 1
        fi
        ;;
        
    6|*)
        echo "No changes made"
        exit 0
        ;;
esac

echo ""
echo "=========================================="
echo "Verification"
echo "=========================================="
echo ""

echo "Current version:"
git describe --tags --always 2>/dev/null || git rev-parse --short HEAD

echo ""
echo "Key files:"
ls -lh main.nf nextflow.config 2>/dev/null

echo ""
echo "If the error persists, run:"
echo "  bash diagnose_pipeline_code.sh"
