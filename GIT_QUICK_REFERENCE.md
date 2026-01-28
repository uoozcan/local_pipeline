# Git Quick Reference for CSC Puhti

Quick command reference for managing the HLA typing pipeline repository on CSC Puhti.

## Initial Setup (One-time)

```bash
# Load Git module
module load git

# Configure Git
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# Set up credential caching (optional, 1 hour)
git config --global credential.helper 'cache --timeout=3600'

# Enable color output
git config --global color.ui auto
```

## Daily Commands

### Check Status
```bash
# See what's changed
git status

# See differences
git diff

# See specific file diff
git diff main.nf
```

### Stage and Commit
```bash
# Stage specific files
git add main.nf nextflow.config

# Stage all changes
git add .

# Stage only modified files (not new files)
git add -u

# Commit with message
git commit -m "Fix: Description of changes"

# Commit all modified files (skip staging)
git commit -am "Update: Description"
```

### Push Changes
```bash
# Push to GitHub
git push

# Push new branch
git push -u origin branch-name

# Force push (dangerous! use only if necessary)
git push --force
```

### Pull Changes
```bash
# Pull latest changes
git pull

# Pull with rebase
git pull --rebase

# Fetch without merging
git fetch origin
```

## Branch Operations

### Create and Switch
```bash
# Create new branch
git branch feature-name

# Switch to branch
git checkout feature-name

# Create and switch in one command
git checkout -b feature-name

# Switch back to main
git checkout main
```

### Merge and Delete
```bash
# Merge branch into current branch
git merge feature-name

# Delete local branch
git branch -d feature-name

# Delete remote branch
git push origin --delete feature-name

# List all branches
git branch -a
```

## Viewing History

```bash
# View commit history
git log

# Compact history
git log --oneline

# Last 10 commits
git log --oneline -10

# Graphical view
git log --graph --oneline --all

# Changes in last commit
git show

# Changes in specific commit
git show abc123
```

## Undoing Changes

### Before Commit
```bash
# Discard changes in working directory
git checkout -- main.nf

# Unstage file (keep changes)
git reset HEAD main.nf

# Discard all changes
git reset --hard
```

### After Commit
```bash
# Undo last commit (keep changes)
git reset --soft HEAD~1

# Undo last commit (discard changes)
git reset --hard HEAD~1

# Amend last commit message
git commit --amend -m "New message"

# Revert a specific commit
git revert abc123
```

## Remote Operations

```bash
# List remotes
git remote -v

# Add remote
git remote add origin https://github.com/user/repo.git

# Change remote URL
git remote set-url origin https://github.com/user/repo.git

# Remove remote
git remote remove origin
```

## Tags

```bash
# List tags
git tag

# Create tag
git tag -a v2.1.0 -m "Version 2.1.0"

# Push tag
git push origin v2.1.0

# Push all tags
git push origin --tags

# Delete local tag
git tag -d v2.1.0

# Delete remote tag
git push origin --delete v2.1.0

# Checkout tag
git checkout v2.1.0
```

## File Management

```bash
# Remove file from Git (keep local)
git rm --cached file.txt

# Remove file from Git and filesystem
git rm file.txt

# Move/rename file
git mv old-name.txt new-name.txt

# Track empty directory
touch dir/.gitkeep
git add dir/.gitkeep
```

## Search and Find

```bash
# Search in files
git grep "search term"

# Search in specific file type
git grep "search term" -- "*.nf"

# Show who modified each line
git blame main.nf

# Find commits by message
git log --grep="bug fix"
```

## Stash (Temporary Save)

```bash
# Save current changes
git stash

# Save with message
git stash save "WIP: feature description"

# List stashes
git stash list

# Apply most recent stash
git stash apply

# Apply and remove stash
git stash pop

# Apply specific stash
git stash apply stash@{2}

# Delete stash
git stash drop stash@{0}

# Clear all stashes
git stash clear
```

## Pipeline-Specific Commands

### Quick Update Workflow
```bash
# Edit files
nano main.nf

# Stage, commit, push
git add main.nf
git commit -m "Fix: issue description"
git push
```

### Before Making Changes
```bash
# Always start with latest code
git pull

# Create feature branch
git checkout -b fix/issue-description

# Make changes...
```

### After Testing
```bash
# Add all changes
git add .

# Commit with good message
git commit -m "Add: New HLA typing tool integration"

# Push to GitHub
git push -u origin fix/issue-description

# Create pull request on GitHub
```

### Updating Documentation
```bash
# Edit docs
nano docs/USAGE.md

# Stage and commit
git add docs/USAGE.md
git commit -m "Update: Usage documentation"
git push
```

## Common Scenarios

### Scenario 1: Fixed a Bug
```bash
git add modules/optitype.nf
git commit -m "Fix: Correct OptiType parameter parsing"
git push
```

### Scenario 2: Added New Feature
```bash
git checkout -b feature/new-tool
# ... make changes ...
git add modules/newtool.nf
git commit -m "Add: Integration of new HLA typing tool"
git push -u origin feature/new-tool
# Create pull request on GitHub
```

### Scenario 3: Update Configuration
```bash
git add conf/puhti.config
git commit -m "Update: Increase memory for SpecHLA"
git push
```

### Scenario 4: Accidentally Committed Sensitive File
```bash
# Remove from staging (before push)
git rm --cached secret.key
echo "secret.key" >> .gitignore
git add .gitignore
git commit -m "Remove sensitive file"

# If already pushed, contact admin immediately
```

### Scenario 5: Merge Conflict
```bash
# After git pull shows conflict
git status  # See conflicted files

# Edit conflicted files, resolve <<<<< ===== >>>>> markers
nano main.nf

# Mark as resolved
git add main.nf

# Complete merge
git commit -m "Merge: Resolve conflict in main.nf"
git push
```

## CSC Puhti Specific

### Load Git
```bash
# Add to .bashrc for automatic loading
echo 'module load git' >> ~/.bashrc

# Or load manually each session
module load git
```

### Working with Large Files
```bash
# Check size before adding
du -h filename

# Don't commit files > 50MB
# Instead, document their location in README
```

### Batch Job Commits
```bash
# Create commit script for automated runs
cat > commit_results.sh << 'EOF'
#!/bin/bash
module load git
cd /scratch/project_2008084/hla-typing-pipeline
git add results/
git commit -m "Results: Batch $(date +%Y%m%d)"
git push
EOF
```

## Troubleshooting

### Problem: "Permission denied"
```bash
# Check file permissions
ls -la ~/.ssh/

# Regenerate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"
```

### Problem: "Authentication failed"
```bash
# For HTTPS: Use Personal Access Token, not password
# For SSH: Check key is added
ssh -T git@github.com
```

### Problem: "Merge conflict"
```bash
# See conflicted files
git status

# Option 1: Resolve manually
nano conflicted-file.txt
git add conflicted-file.txt
git commit

# Option 2: Use theirs
git checkout --theirs conflicted-file.txt
git add conflicted-file.txt
git commit

# Option 3: Use ours
git checkout --ours conflicted-file.txt
git add conflicted-file.txt
git commit
```

### Problem: "Detached HEAD"
```bash
# Create branch from current state
git checkout -b recovery-branch

# Or return to main
git checkout main
```

## Help Commands

```bash
# General help
git help

# Command-specific help
git help commit
git help push

# Quick help
git commit --help
```

## Pro Tips

1. **Always commit before pulling**: Avoid merge conflicts
2. **Write descriptive messages**: Future you will thank you
3. **Commit often**: Small, focused commits are better
4. **Use branches**: Keep main stable
5. **Pull before push**: Stay synchronized
6. **Test before commit**: Broken code shouldn't be committed
7. **Review before push**: Check what you're about to share

## Emergency Commands

### Undo Everything (Nuclear Option)
```bash
# Discard ALL local changes
git reset --hard HEAD
git clean -fd

# Sync with remote completely
git fetch origin
git reset --hard origin/main
```

### Recover Deleted File
```bash
# Find commit where file existed
git log -- path/to/file

# Restore from specific commit
git checkout abc123 -- path/to/file
```

### Recover Deleted Branch
```bash
# Find the commit
git reflog

# Recreate branch
git branch recovered-branch abc123
```

---

## Quick Reference Card Summary

```bash
# Daily workflow
git pull                          # Get latest
git add .                         # Stage all
git commit -m "message"           # Commit
git push                          # Share

# Check status
git status                        # What changed
git log --oneline                 # History

# Undo
git checkout -- file              # Discard changes
git reset HEAD file               # Unstage
git reset --soft HEAD~1           # Undo commit

# Branches
git checkout -b name              # Create & switch
git checkout main                 # Switch to main
git merge name                    # Merge branch

# Remote
git remote -v                     # List remotes
git push origin branch            # Push branch
git pull                          # Get updates
```

---

**Keep this file handy for quick reference!**

For detailed help: `git help <command>`  
For emergency: Read "Troubleshooting" section above
