#!/bin/bash

# Initialize Git repository for HLA typing pipeline
# Usage: bash initialize_git_repo.sh

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Configuration
REPO_DIR="/scratch/project_2008084/hla-typing-pipeline"

if [ ! -d "$REPO_DIR" ]; then
    log_error "Repository directory not found: $REPO_DIR"
    log_info "Please run reorganize_for_github.sh first"
    exit 1
fi

cd "$REPO_DIR"

# ============================================================================
# STEP 1: Initialize Git
# ============================================================================
log_info "Initializing Git repository..."

if [ -d ".git" ]; then
    log_warn "Git repository already initialized"
else
    git init
    log_info "Git repository initialized"
fi

# ============================================================================
# STEP 2: Configure Git (CSC Puhti specific)
# ============================================================================
log_info "Configuring Git for CSC Puhti..."

# Set core.sharedRepository for group collaboration
git config core.sharedRepository group

# Increase buffer for large files
git config http.postBuffer 524288000

# Set up Git LFS if available (for large files)
if command -v git-lfs &> /dev/null; then
    log_info "Git LFS is available, initializing..."
    git lfs install
else
    log_warn "Git LFS not found. Large files should be stored separately."
fi

# ============================================================================
# STEP 3: Check .gitignore
# ============================================================================
log_info "Verifying .gitignore..."

if [ ! -f ".gitignore" ]; then
    log_error ".gitignore not found!"
    exit 1
fi

# Count lines in .gitignore
GITIGNORE_LINES=$(wc -l < .gitignore)
log_info ".gitignore has $GITIGNORE_LINES lines"

# ============================================================================
# STEP 4: Check for sensitive files
# ============================================================================
log_info "Checking for sensitive files..."

SENSITIVE_PATTERNS=(
    "*.sec"
    "*.key"
    "*_key.*"
    ".passphrase"
    ".c4gh_passphrase"
    "*.c4gh"
)

SENSITIVE_FOUND=0
for pattern in "${SENSITIVE_PATTERNS[@]}"; do
    if find . -name "$pattern" -type f | grep -q .; then
        log_warn "Found sensitive files matching: $pattern"
        SENSITIVE_FOUND=1
    fi
done

if [ $SENSITIVE_FOUND -eq 1 ]; then
    log_error "Sensitive files found! Remove them before committing."
    log_info "Run: find . -name '*.sec' -o -name '*.key' | xargs rm"
    exit 1
fi

# ============================================================================
# STEP 5: Check for large files
# ============================================================================
log_info "Checking for large files (>10MB)..."

LARGE_FILES=$(find . -type f -size +10M ! -path "./.git/*" ! -path "./work/*" 2>/dev/null | head -20)

if [ -n "$LARGE_FILES" ]; then
    log_warn "Large files found that may not be suitable for Git:"
    echo "$LARGE_FILES"
    log_info "Consider adding these to .gitignore or using Git LFS"
fi

# ============================================================================
# STEP 6: Initial commit
# ============================================================================
log_info "Creating initial commit..."

# Stage all files
git add .

# Check what will be committed
log_info "Files to be committed:"
git status --short | head -20
echo "..."
TOTAL_FILES=$(git status --short | wc -l)
log_info "Total files staged: $TOTAL_FILES"

# Ask for confirmation
echo ""
echo -e "${BLUE}Ready to create initial commit.${NC}"
read -p "Proceed? (yes/no): " confirm

if [ "$confirm" != "yes" ]; then
    log_warn "Commit cancelled"
    exit 0
fi

# Create initial commit
git commit -m "Initial commit: HLA typing pipeline

- Nextflow DSL2 pipeline for HLA typing
- Multiple HLA typing tools integrated
- Confidence-weighted majority voting
- Optimized for CSC Puhti HPC environment
- Containerized with Singularity
- Comprehensive documentation

Pipeline version: 2.0.0"

log_info "Initial commit created successfully!"

# ============================================================================
# STEP 7: Create tags
# ============================================================================
log_info "Creating version tag..."

git tag -a v2.0.0 -m "Version 2.0.0 - Initial GitHub release

Features:
- Multi-tool HLA typing (9 tools)
- Majority voting consensus
- BAM and FASTQ support
- SLURM integration
- Comprehensive documentation"

log_info "Tag v2.0.0 created"

# ============================================================================
# STEP 8: GitHub preparation
# ============================================================================
log_info "Preparing GitHub instructions..."

cat > GITHUB_SETUP.md << 'EOF'
# GitHub Repository Setup

## Prerequisites

1. **Create GitHub account** (if you don't have one)
   - Go to https://github.com
   - Sign up for a free account

2. **Create new repository on GitHub**
   - Go to https://github.com/new
   - Repository name: `hla-typing-pipeline`
   - Description: "Nextflow pipeline for HLA typing from NGS data with multi-tool consensus calling"
   - Choose: Public or Private
   - **DO NOT** initialize with README, .gitignore, or license (we already have these)
   - Click "Create repository"

## Push to GitHub

After creating the repository on GitHub, run these commands from Puhti:

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Add GitHub as remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/hla-typing-pipeline.git

# Verify remote
git remote -v

# Push to GitHub (first time)
git branch -M main
git push -u origin main

# Push tags
git push origin --tags
```

## Authentication

### Option 1: Personal Access Token (Recommended)

1. **Generate token on GitHub**:
   - Go to GitHub Settings → Developer settings → Personal access tokens → Tokens (classic)
   - Generate new token
   - Select scopes: `repo` (full control of private repositories)
   - Copy the token (you won't see it again!)

2. **Use token for authentication**:
   ```bash
   # When prompted for password, paste your token
   git push -u origin main
   ```

3. **Cache credentials** (optional, for convenience):
   ```bash
   git config --global credential.helper 'cache --timeout=3600'
   ```

### Option 2: SSH Key

1. **Generate SSH key on Puhti**:
   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   cat ~/.ssh/id_ed25519.pub
   ```

2. **Add key to GitHub**:
   - Go to GitHub Settings → SSH and GPG keys → New SSH key
   - Paste the public key

3. **Change remote to SSH**:
   ```bash
   git remote set-url origin git@github.com:YOUR_USERNAME/hla-typing-pipeline.git
   ```

## Post-Push Steps

1. **Add repository description** on GitHub
2. **Add topics/tags**: nextflow, hla-typing, bioinformatics, genomics, ngs
3. **Set up GitHub Pages** (optional, for documentation)
4. **Add collaborators** (if applicable)
5. **Create issues** for known TODOs
6. **Set up releases** for versioning

## Repository Settings (Recommended)

On GitHub repository settings:

1. **Branches**:
   - Set `main` as default branch
   - Add branch protection rules (optional)

2. **Actions** (if using CI/CD):
   - Enable GitHub Actions

3. **Security**:
   - Enable vulnerability alerts
   - Enable Dependabot (for dependency updates)

## Making Changes

After initial push, normal Git workflow:

```bash
# Make changes
vim main.nf

# Stage and commit
git add main.nf
git commit -m "Update: description of changes"

# Push to GitHub
git push
```

## Troubleshooting

### Large files error
If you get "file too large" error:
```bash
# Add to .gitignore and remove from git
echo "large_file.fa" >> .gitignore
git rm --cached large_file.fa
git commit -m "Remove large file"
```

### Authentication failed
- Check your token/SSH key is correct
- For token: use the token as password, not your GitHub password
- For SSH: ensure key is added to ssh-agent

### Push rejected
If you have divergent branches:
```bash
git pull --rebase origin main
git push
```
EOF

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=================================================="
echo -e "${GREEN}Git Repository Setup Complete!${NC}"
echo "=================================================="
echo ""
echo "Repository: $REPO_DIR"
echo "Branch: $(git branch --show-current)"
echo "Commits: $(git log --oneline | wc -l)"
echo "Tags: $(git tag | tr '\n' ' ')"
echo ""
echo -e "${BLUE}Next Steps:${NC}"
echo "1. Review GITHUB_SETUP.md for GitHub instructions"
echo "2. Create repository on GitHub"
echo "3. Push: git remote add origin https://github.com/YOUR_USERNAME/hla-typing-pipeline.git"
echo "4. Push: git push -u origin main"
echo "5. Push tags: git push origin --tags"
echo ""
echo "Documentation:"
echo "  cat GITHUB_SETUP.md"
echo ""
echo "=================================================="
