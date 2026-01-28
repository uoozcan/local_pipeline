# Deployment Checklist for GitHub Release

This checklist will guide you through deploying your HLA typing pipeline to GitHub.

## ✅ Pre-Release Checklist

### 1. Repository Setup
- [ ] Create new GitHub repository (public or private)
- [ ] Initialize local git repository
  ```bash
  cd hla-typing-pipeline
  git init
  ```
- [ ] Add remote origin
  ```bash
  git remote add origin https://github.com/yourusername/hla-typing-pipeline.git
  ```

### 2. Update Documentation
- [ ] Replace `yourusername` with your GitHub username in:
  - [ ] README.md
  - [ ] nextflow.config (manifest section)
  - [ ] params.yaml.example
- [ ] Update contact information in README.md
- [ ] Add your email address
- [ ] Verify all links work

### 3. Update Container References
If NOT using CSC Puhti containers:
- [ ] Update container paths in nextflow.config
- [ ] Build or pull Docker images
- [ ] Convert to Singularity if needed
- [ ] Update documentation with new paths

### 4. Test the Pipeline

**Local Test (Docker)**:
```bash
# With small test dataset
nextflow run main.nf \
    --input test_data/ \
    --input_type fastq \
    --tools optitype \
    -profile docker
```

**CSC Puhti Test**:
```bash
# Submit test job
sbatch submit_slurm.sh
```

- [ ] Test with FASTQ input
- [ ] Test with BAM input
- [ ] Test with each tool individually
- [ ] Test with majority voting
- [ ] Verify all output files are created
- [ ] Check log files for errors

### 5. Version Control
- [ ] Update version number in:
  - [ ] nextflow.config (manifest.version)
  - [ ] README.md
  - [ ] CHANGELOG.md
- [ ] Verify CHANGELOG.md is up to date
- [ ] Tag release appropriately

### 6. Code Quality
- [ ] Remove any hardcoded paths
- [ ] Remove sensitive information
- [ ] Check for TODO comments
- [ ] Verify all modules are properly included
- [ ] Test all execution profiles

## 🚀 Release Steps

### Step 1: Initial Commit

```bash
# Add all files
git add .

# Check what's being added
git status

# Make initial commit
git commit -m "Initial release v2.0.0

- Multi-tool HLA typing pipeline
- Support for OptiType, ArcasHLA, SpecHLA
- Majority voting consensus calling
- Optimized for CSC Puhti HPC
- Comprehensive documentation"

# Push to GitHub
git push -u origin main
```

### Step 2: Create Release on GitHub

1. Go to your repository on GitHub
2. Click "Releases" → "Create a new release"
3. Tag version: `v2.0.0`
4. Release title: `HLA Typing Pipeline v2.0.0`
5. Description:
   ```
   # HLA Typing Multi-Tool Pipeline v2.0.0
   
   First stable release of the comprehensive HLA typing pipeline.
   
   ## Features
   - Multi-tool integration (OptiType, ArcasHLA, SpecHLA)
   - Majority voting consensus calling
   - Flexible input support (BAM, CRAM, FASTQ)
   - Optimized for HPC environments (CSC Puhti/SLURM)
   - Docker and Singularity support
   - Comprehensive documentation
   
   ## Quick Start
   ```bash
   nextflow run yourusername/hla-typing-pipeline \
       --input samples/ \
       --input_type fastq \
       --tools optitype,arcashla \
       -profile docker
   ```
   
   See [README.md](README.md) for complete documentation.
   
   ## Known Issues
   - SpecHLA requires `--spechla_exon_only 1` for exome data
   
   ## Citation
   If you use this pipeline, please cite the individual tools (see README).
   ```
6. Attach any additional files if needed
7. Click "Publish release"

### Step 3: Post-Release

- [ ] Test installation from GitHub:
  ```bash
  nextflow run yourusername/hla-typing-pipeline \
      --input test_data/ \
      --input_type fastq \
      --tools optitype \
      -profile docker
  ```
- [ ] Update any external documentation
- [ ] Announce release (if applicable)
- [ ] Monitor issue tracker for bug reports

## 📦 Optional: Container Registry Setup

If you want to host your own containers:

### Docker Hub

1. Create account on https://hub.docker.com
2. Create repository for each tool
3. Build and push containers:
   ```bash
   docker build -t yourusername/optitype:latest optitype/
   docker push yourusername/optitype:latest
   ```
4. Update nextflow.config with new paths

### GitHub Container Registry

1. Create GitHub Personal Access Token
2. Login to GHCR:
   ```bash
   echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin
   ```
3. Build and push:
   ```bash
   docker build -t ghcr.io/yourusername/optitype:latest optitype/
   docker push ghcr.io/yourusername/optitype:latest
   ```
4. Update nextflow.config

## 🔧 Maintenance

### Version Updates

For future updates:
1. Create new branch: `git checkout -b feature/new-feature`
2. Make changes
3. Update CHANGELOG.md
4. Update version numbers
5. Test thoroughly
6. Create pull request
7. Merge and create new release

### Bug Fixes

For bug fixes:
1. Create issue on GitHub
2. Create branch: `git checkout -b bugfix/issue-123`
3. Fix bug
4. Test fix
5. Update CHANGELOG.md
6. Create pull request
7. Create patch release (e.g., v2.0.1)

## 📋 Useful Git Commands

```bash
# Check status
git status

# View commit history
git log --oneline

# Create new branch
git checkout -b branch-name

# Switch branches
git checkout main

# Pull latest changes
git pull origin main

# View remote
git remote -v

# Tag a release
git tag -a v2.0.0 -m "Version 2.0.0"
git push origin v2.0.0

# Undo last commit (local only)
git reset --soft HEAD~1

# Discard all local changes
git reset --hard HEAD
```

## ✨ Success Checklist

Before announcing your release, verify:
- [ ] Pipeline installs cleanly from GitHub
- [ ] All documentation links work
- [ ] README examples are tested and work
- [ ] Container images are accessible
- [ ] License file is present
- [ ] .gitignore prevents committing sensitive data
- [ ] CHANGELOG is complete
- [ ] Version numbers are consistent
- [ ] Contact information is correct

## 🎉 You're Ready to Release!

Once all items are checked, your pipeline is ready for release. Good luck! 🚀

---

**Questions?**
- GitHub Docs: https://docs.github.com
- Nextflow Best Practices: https://www.nextflow.io/docs/latest/
- nf-core Standards: https://nf-co.re/developers/guidelines
