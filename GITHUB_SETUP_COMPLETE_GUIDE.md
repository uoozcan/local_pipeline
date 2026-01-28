# HLA Typing Pipeline - GitHub Repository Setup Guide

Complete guide for reorganizing the HLA typing pipeline on CSC Puhti and publishing to GitHub.

## Overview

This guide will help you:
1. Reorganize your current pipeline structure
2. Prepare it for GitHub
3. Initialize Git repository
4. Push to GitHub
5. Maintain the repository

**Time required:** ~1-2 hours  
**Skills needed:** Basic Linux commands, Git basics

---

## Prerequisites

### On CSC Puhti

```bash
# Check you're in the project directory
cd /scratch/project_2008084/hla_rnaseq_analysis

# Load Git module
module load git

# Verify Git is available
git --version
```

### On GitHub

1. Create a GitHub account (if you don't have one): https://github.com/join
2. Have your authentication ready (Personal Access Token or SSH key)

---

## Phase 1: Reorganize Pipeline Structure

### Step 1.1: Upload Scripts to Puhti

Transfer the reorganization scripts to Puhti:

```bash
cd /scratch/project_2008084/hla_rnaseq_analysis

# Upload the scripts:
# - reorganize_for_github.sh
# - initialize_git_repo.sh
# - verify_repo_structure.sh

# Make them executable
chmod +x reorganize_for_github.sh initialize_git_repo.sh verify_repo_structure.sh
```

### Step 1.2: Review Current Structure

```bash
# See what you currently have
tree -L 2 -d

# Check for sensitive files (should be removed before GitHub)
find . -name "*.sec" -o -name "*.key" -o -name ".passphrase"
```

**⚠️ IMPORTANT:** Before proceeding, backup your current work:

```bash
# Create a backup
tar -czf ~/hla_pipeline_backup_$(date +%Y%m%d).tar.gz \
  main.nf nextflow.config modules/ bin/ conf/ docs/ README.md

# Verify backup
tar -tzf ~/hla_pipeline_backup_$(date +%Y%m%d).tar.gz | head -20
```

### Step 1.3: Run Reorganization

```bash
# Submit reorganization job
sbatch reorganize_for_github.sh

# Monitor job
squeue -u $USER

# Check results when complete
tail logs/slurm-*.out
```

This will:
- Create new repository at `/scratch/project_2008084/hla-typing-pipeline`
- Archive old files to `/scratch/project_2008084/hla_archive_YYYYMMDD`
- Organize files in GitHub-friendly structure
- Create necessary documentation

### Step 1.4: Verify New Structure

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Check structure
tree -L 2

# Should see:
# .
# ├── main.nf
# ├── nextflow.config
# ├── modules/
# ├── bin/
# ├── conf/
# ├── docs/
# ├── examples/
# ├── scripts/
# ├── README.md
# ├── LICENSE
# ├── .gitignore
# └── REORGANIZATION_NOTES.md
```

---

## Phase 2: Prepare for Git

### Step 2.1: Review Files

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Check .gitignore is comprehensive
cat .gitignore

# Verify no sensitive files
find . -name "*.sec" -o -name "*.key" -o -name ".passphrase"

# Check for large files (>10MB)
find . -type f -size +10M ! -path "./.git/*"
```

### Step 2.2: Update Configuration Paths

Edit configuration files to use absolute paths or environment variables:

```bash
# Edit conf/puhti.config
nano conf/puhti.config

# Update paths:
# - Singularity cache: /scratch/project_2008084/hla_references/singularity_cache/containers
# - Reference genomes: /scratch/project_2008084/hla_references

# Edit conf/params.config
nano conf/params.config
```

### Step 2.3: Update Documentation

```bash
# Update installation guide
nano docs/INSTALLATION.md

# Make sure paths are correct:
# - Reference data location
# - Container cache location
# - Project directory structure
```

### Step 2.4: Run Verification

```bash
# Run verification script
bash verify_repo_structure.sh

# Review output
# Fix any errors before proceeding
```

---

## Phase 3: Initialize Git Repository

### Step 3.1: Load Git Module

```bash
module load git

# Configure Git (first time only)
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### Step 3.2: Run Git Initialization

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Run initialization script
bash initialize_git_repo.sh

# This will:
# - Initialize Git repository
# - Create initial commit
# - Create version tag (v2.0.0)
# - Generate GITHUB_SETUP.md
```

### Step 3.3: Review Git Status

```bash
# Check Git status
git status

# View commit history
git log --oneline

# Check tags
git tag -l

# See what files are tracked
git ls-files | head -20
```

---

## Phase 4: Create GitHub Repository

### Step 4.1: Create Repository on GitHub

1. Go to https://github.com/new
2. Fill in details:
   - **Repository name:** `hla-typing-pipeline`
   - **Description:** "Nextflow pipeline for HLA typing from NGS data with multi-tool consensus calling"
   - **Visibility:** Choose Public or Private
   - ⚠️ **DO NOT** initialize with README, .gitignore, or license
3. Click "Create repository"

### Step 4.2: Set Up Authentication

#### Option A: Personal Access Token (Recommended)

1. On GitHub: Settings → Developer settings → Personal access tokens → Tokens (classic)
2. Click "Generate new token (classic)"
3. Select scopes: `repo` (full control)
4. Copy the token (save it securely!)

#### Option B: SSH Key

```bash
# Generate SSH key on Puhti
ssh-keygen -t ed25519 -C "your.email@example.com"

# Display public key
cat ~/.ssh/id_ed25519.pub

# Copy and add to GitHub:
# Settings → SSH and GPG keys → New SSH key
```

---

## Phase 5: Push to GitHub

### Step 5.1: Add Remote

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# For HTTPS (with token):
git remote add origin https://github.com/YOUR_USERNAME/hla-typing-pipeline.git

# OR for SSH:
git remote add origin git@github.com:YOUR_USERNAME/hla-typing-pipeline.git

# Verify remote
git remote -v
```

### Step 5.2: Push Code

```bash
# Push main branch
git push -u origin main

# If using HTTPS token, enter:
# Username: your-github-username
# Password: your-personal-access-token

# Push tags
git push origin --tags
```

### Step 5.3: Verify on GitHub

Visit your repository: `https://github.com/YOUR_USERNAME/hla-typing-pipeline`

Check:
- ✓ Files are visible
- ✓ README displays properly
- ✓ No sensitive files present
- ✓ Tags are visible

---

## Phase 6: Configure GitHub Repository

### Step 6.1: Add Repository Description

On GitHub repository page:
1. Click ⚙️ (Settings icon) next to "About"
2. Add description
3. Add topics: `nextflow`, `bioinformatics`, `hla-typing`, `genomics`, `ngs`, `immunogenetics`
4. Add website (if applicable)

### Step 6.2: Configure Settings

Go to Settings tab:

1. **General:**
   - Set default branch: `main`
   - Enable issues
   - Enable wikis (optional)

2. **Branches:**
   - Add branch protection for `main` (optional):
     - Require pull request reviews
     - Require status checks to pass

3. **Security:**
   - Enable Dependabot alerts
   - Enable secret scanning

### Step 6.3: Create Release

1. Go to "Releases" → "Create a new release"
2. Choose tag: `v2.0.0`
3. Release title: "Version 2.0.0 - Initial Release"
4. Description:
   ```markdown
   # HLA Typing Pipeline v2.0.0
   
   Initial GitHub release of the comprehensive HLA typing pipeline.
   
   ## Features
   - Multi-tool HLA typing (OptiType, ArcasHLA, SpecHLA, xHLA, HLA-HD, HLA*LA, Seq2HLA, HLAProfiler, Kourami)
   - Confidence-weighted majority voting consensus
   - Support for BAM and FASTQ inputs
   - Optimized for CSC Puhti HPC
   - Containerized with Singularity
   
   ## Installation
   See [INSTALLATION.md](docs/INSTALLATION.md) for setup instructions.
   
   ## Usage
   See [USAGE.md](docs/USAGE.md) for usage examples.
   ```
5. Click "Publish release"

---

## Phase 7: Post-Setup

### Step 7.1: Test Clone

Test that others can clone your repository:

```bash
# From a different directory
cd /tmp
git clone https://github.com/YOUR_USERNAME/hla-typing-pipeline.git
cd hla-typing-pipeline
tree -L 2
```

### Step 7.2: Add Collaborators (Optional)

If working with others:
1. Settings → Collaborators
2. Add by username or email
3. Set appropriate permissions

### Step 7.3: Set Up Issues

Create initial issues for known TODOs:

1. Go to "Issues" tab
2. Click "New issue"
3. Examples:
   - "Add GitHub Actions CI/CD"
   - "Create container building workflow"
   - "Add comprehensive test dataset"

---

## Maintaining the Repository

### Daily Workflow

```bash
cd /scratch/project_2008084/hla-typing-pipeline

# Make changes
nano main.nf

# Check status
git status

# Stage changes
git add main.nf

# Commit with descriptive message
git commit -m "Fix: Correct HLA extraction parameters"

# Push to GitHub
git push
```

### Creating New Versions

```bash
# Update CHANGELOG.md
nano CHANGELOG.md

# Commit changes
git add CHANGELOG.md
git commit -m "Release version 2.1.0"

# Create tag
git tag -a v2.1.0 -m "Version 2.1.0 - Bug fixes and improvements"

# Push changes and tags
git push
git push origin v2.1.0

# Create release on GitHub
```

### Branch Workflow (Advanced)

```bash
# Create feature branch
git checkout -b feature/new-tool

# Make changes
# ... edit files ...

# Commit changes
git add .
git commit -m "Add new HLA typing tool"

# Push branch
git push -u origin feature/new-tool

# Create pull request on GitHub
```

---

## Troubleshooting

### Large File Error

```bash
# Error: file too large
# Solution: Add to .gitignore
echo "large_file.fa" >> .gitignore
git rm --cached large_file.fa
git commit -m "Remove large file from tracking"
git push
```

### Authentication Failed

```bash
# Check remote URL
git remote -v

# Update remote URL
git remote set-url origin https://github.com/YOUR_USERNAME/hla-typing-pipeline.git

# For HTTPS, use token as password, not GitHub password
# For SSH, check key is added: ssh -T git@github.com
```

### Divergent Branches

```bash
# Pull with rebase
git pull --rebase origin main

# Resolve conflicts if any
# ... edit files ...
git add .
git rebase --continue

# Push
git push
```

### Sensitive File Committed

```bash
# Remove from history (use with caution!)
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch path/to/sensitive/file" \
  --prune-empty --tag-name-filter cat -- --all

# Force push (dangerous!)
git push origin --force --all
```

---

## Best Practices

### Commit Messages

Use clear, descriptive commit messages:

```bash
# Good
git commit -m "Fix: Correct chromosome naming in SpecHLA module"
git commit -m "Add: Support for WGS data in OptiType"
git commit -m "Update: Documentation for CSC Puhti setup"

# Bad
git commit -m "fix"
git commit -m "updates"
git commit -m "changes"
```

### .gitignore Management

Regularly review and update `.gitignore`:

```bash
# Add patterns for new file types
echo "*.tmp" >> .gitignore
echo "new_results_dir/" >> .gitignore

# Commit .gitignore changes
git add .gitignore
git commit -m "Update .gitignore: Add temporary file patterns"
git push
```

### Documentation

Keep documentation up to date:
- Update README when adding features
- Update CHANGELOG for each version
- Add comments to code
- Create issues for known limitations

---

## Support

### Getting Help

1. **Documentation:** Check `docs/` directory
2. **Issues:** Create issue on GitHub
3. **CSC Support:** servicedesk@csc.fi
4. **Git Help:** `git help <command>`

### Useful Resources

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Git Documentation](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)
- [CSC Puhti User Guide](https://docs.csc.fi/computing/systems-puhti/)

---

## Checklist

Use this checklist to track your progress:

- [ ] Phase 1: Reorganize pipeline structure
  - [ ] Upload scripts to Puhti
  - [ ] Create backup
  - [ ] Run reorganization script
  - [ ] Verify new structure
  
- [ ] Phase 2: Prepare for Git
  - [ ] Review .gitignore
  - [ ] Check for sensitive files
  - [ ] Update configuration paths
  - [ ] Update documentation
  - [ ] Run verification script
  
- [ ] Phase 3: Initialize Git
  - [ ] Configure Git
  - [ ] Run initialization script
  - [ ] Review Git status
  
- [ ] Phase 4: Create GitHub repository
  - [ ] Create repository on GitHub
  - [ ] Set up authentication (token or SSH)
  
- [ ] Phase 5: Push to GitHub
  - [ ] Add remote
  - [ ] Push code
  - [ ] Push tags
  - [ ] Verify on GitHub
  
- [ ] Phase 6: Configure GitHub
  - [ ] Add repository description
  - [ ] Configure settings
  - [ ] Create release
  
- [ ] Phase 7: Post-setup
  - [ ] Test clone
  - [ ] Add collaborators (if applicable)
  - [ ] Create initial issues

---

**Last Updated:** January 2026  
**Pipeline Version:** 2.0.0  
**Author:** Umut Özcan
