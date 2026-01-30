#!/bin/bash

# Verify repository structure and readiness for GitHub
# Usage: bash verify_repo_structure.sh

set -euo pipefail

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

REPO_DIR="/scratch/project_2008084/hla-typing-pipeline"
ERRORS=0
WARNINGS=0

check_pass() {
    echo -e "  ${GREEN}✓${NC} $1"
}

check_fail() {
    echo -e "  ${RED}✗${NC} $1"
    ((ERRORS++))
}

check_warn() {
    echo -e "  ${YELLOW}⚠${NC} $1"
    ((WARNINGS++))
}

check_info() {
    echo -e "  ${BLUE}ℹ${NC} $1"
}

echo "=================================================="
echo "HLA Typing Pipeline - Repository Verification"
echo "=================================================="
echo ""

if [ ! -d "$REPO_DIR" ]; then
    echo -e "${RED}ERROR: Repository directory not found: $REPO_DIR${NC}"
    exit 1
fi

cd "$REPO_DIR"

# ============================================================================
# 1. Essential Files
# ============================================================================
echo -e "${BLUE}[1] Checking essential files...${NC}"

ESSENTIAL_FILES=(
    "main.nf"
    "nextflow.config"
    "README.md"
    "LICENSE"
    ".gitignore"
)

for file in "${ESSENTIAL_FILES[@]}"; do
    if [ -f "$file" ]; then
        check_pass "$file exists"
    else
        check_fail "$file missing"
    fi
done

echo ""

# ============================================================================
# 2. Directory Structure
# ============================================================================
echo -e "${BLUE}[2] Checking directory structure...${NC}"

REQUIRED_DIRS=(
    "modules"
    "bin"
    "conf"
    "docs"
)

OPTIONAL_DIRS=(
    "examples"
    "scripts"
)

for dir in "${REQUIRED_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        COUNT=$(find "$dir" -type f | wc -l)
        check_pass "$dir/ exists ($COUNT files)"
    else
        check_fail "$dir/ missing"
    fi
done

for dir in "${OPTIONAL_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        COUNT=$(find "$dir" -type f | wc -l)
        check_pass "$dir/ exists ($COUNT files)"
    else
        check_warn "$dir/ not found (optional)"
    fi
done

echo ""

# ============================================================================
# 3. Modules
# ============================================================================
echo -e "${BLUE}[3] Checking Nextflow modules...${NC}"

if [ -d "modules" ]; then
    MODULE_COUNT=$(find modules/ -name "*.nf" | wc -l)
    check_info "Found $MODULE_COUNT Nextflow modules"
    
    EXPECTED_MODULES=(
        "optitype.nf"
        "arcashla.nf"
        "spechla.nf"
        "majority_voting.nf"
    )
    
    for module in "${EXPECTED_MODULES[@]}"; do
        if [ -f "modules/$module" ]; then
            check_pass "modules/$module"
        else
            check_warn "modules/$module not found"
        fi
    done
else
    check_fail "modules/ directory missing"
fi

echo ""

# ============================================================================
# 4. Documentation
# ============================================================================
echo -e "${BLUE}[4] Checking documentation...${NC}"

# Check README
if [ -f "README.md" ]; then
    README_LINES=$(wc -l < README.md)
    if [ $README_LINES -gt 20 ]; then
        check_pass "README.md is comprehensive ($README_LINES lines)"
    else
        check_warn "README.md is short ($README_LINES lines)"
    fi
else
    check_fail "README.md missing"
fi

# Check documentation files
if [ -d "docs" ]; then
    DOC_COUNT=$(find docs/ -name "*.md" | wc -l)
    check_info "Found $DOC_COUNT documentation files"
    
    IMPORTANT_DOCS=(
        "INSTALLATION.md"
        "USAGE.md"
        "QUICKSTART.md"
    )
    
    for doc in "${IMPORTANT_DOCS[@]}"; do
        if [ -f "docs/$doc" ]; then
            check_pass "docs/$doc"
        else
            check_warn "docs/$doc not found"
        fi
    done
fi

echo ""

# ============================================================================
# 5. Configuration Files
# ============================================================================
echo -e "${BLUE}[5] Checking configuration files...${NC}"

if [ -f "nextflow.config" ]; then
    if grep -q "profiles" nextflow.config; then
        check_pass "nextflow.config has profiles"
    else
        check_warn "nextflow.config missing profiles section"
    fi
fi

if [ -d "conf" ]; then
    if [ -f "conf/puhti.config" ]; then
        check_pass "conf/puhti.config exists"
    else
        check_warn "conf/puhti.config not found"
    fi
    
    if [ -f "conf/params.config" ]; then
        check_pass "conf/params.config exists"
    else
        check_warn "conf/params.config not found"
    fi
fi

echo ""

# ============================================================================
# 6. Sensitive Files Check
# ============================================================================
echo -e "${BLUE}[6] Checking for sensitive files...${NC}"

SENSITIVE_FOUND=0

# Check for key files
if find . -type f \( -name "*.sec" -o -name "*.key" -o -name ".passphrase" \) ! -path "./.git/*" | grep -q .; then
    check_fail "Found sensitive key files (*.sec, *.key, .passphrase)"
    SENSITIVE_FOUND=1
else
    check_pass "No sensitive key files found"
fi

# Check for large data files
if find . -type f \( -name "*.bam" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) ! -path "./.git/*" | grep -q .; then
    check_warn "Found large data files (BAM/FASTQ)"
    check_info "These should be added to .gitignore"
else
    check_pass "No large data files in repository"
fi

echo ""

# ============================================================================
# 7. .gitignore Verification
# ============================================================================
echo -e "${BLUE}[7] Verifying .gitignore...${NC}"

if [ -f ".gitignore" ]; then
    GITIGNORE_LINES=$(wc -l < .gitignore)
    check_pass ".gitignore exists ($GITIGNORE_LINES lines)"
    
    # Check for essential patterns
    ESSENTIAL_PATTERNS=(
        "work/"
        ".nextflow/"
        "*.log"
        "results/"
        "*.sec"
        "*.key"
    )
    
    for pattern in "${ESSENTIAL_PATTERNS[@]}"; do
        if grep -q "$pattern" .gitignore; then
            check_pass "Ignores: $pattern"
        else
            check_warn "Missing in .gitignore: $pattern"
        fi
    done
else
    check_fail ".gitignore missing"
fi

echo ""

# ============================================================================
# 8. File Size Check
# ============================================================================
echo -e "${BLUE}[8] Checking file sizes...${NC}"

# Find files larger than 10MB
LARGE_FILES=$(find . -type f -size +10M ! -path "./.git/*" 2>/dev/null | wc -l)

if [ $LARGE_FILES -eq 0 ]; then
    check_pass "No files larger than 10MB"
else
    check_warn "Found $LARGE_FILES files larger than 10MB"
    find . -type f -size +10M ! -path "./.git/*" 2>/dev/null | head -5 | while read file; do
        SIZE=$(du -h "$file" | cut -f1)
        check_info "  $file ($SIZE)"
    done
fi

# Total repository size
TOTAL_SIZE=$(du -sh . 2>/dev/null | cut -f1)
check_info "Total repository size: $TOTAL_SIZE"

echo ""

# ============================================================================
# 9. Git Status
# ============================================================================
echo -e "${BLUE}[9] Git status...${NC}"

if [ -d ".git" ]; then
    check_pass "Git repository initialized"
    
    BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
    check_info "Current branch: $BRANCH"
    
    COMMITS=$(git log --oneline 2>/dev/null | wc -l)
    check_info "Commits: $COMMITS"
    
    TAGS=$(git tag 2>/dev/null | wc -l)
    check_info "Tags: $TAGS"
    
    # Check for uncommitted changes
    if git diff-index --quiet HEAD -- 2>/dev/null; then
        check_pass "No uncommitted changes"
    else
        MODIFIED=$(git status --short | wc -l)
        check_warn "$MODIFIED files with uncommitted changes"
    fi
    
    # Check for remote
    if git remote -v | grep -q "origin"; then
        REMOTE=$(git remote get-url origin 2>/dev/null)
        check_pass "Remote configured: $REMOTE"
    else
        check_info "No remote configured yet (run git remote add origin <url>)"
    fi
else
    check_warn "Not a Git repository (run: git init)"
fi

echo ""

# ============================================================================
# 10. Example Files
# ============================================================================
echo -e "${BLUE}[10] Checking example files...${NC}"

if [ -d "examples" ]; then
    if [ -f "examples/samplesheet_example.csv" ]; then
        check_pass "Example samplesheet exists"
    else
        check_warn "No example samplesheet"
    fi
    
    if [ -f "examples/params.yaml.example" ]; then
        check_pass "Example params file exists"
    else
        check_warn "No example params file"
    fi
else
    check_warn "examples/ directory not found"
fi

echo ""

# ============================================================================
# Summary
# ============================================================================
echo "=================================================="
echo "Verification Summary"
echo "=================================================="

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}✓ Repository is ready for GitHub!${NC}"
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}⚠ Repository is mostly ready ($WARNINGS warnings)${NC}"
    echo "  Review warnings above before pushing to GitHub"
else
    echo -e "${RED}✗ Repository has issues ($ERRORS errors, $WARNINGS warnings)${NC}"
    echo "  Fix errors before pushing to GitHub"
fi

echo ""
echo "Statistics:"
echo "  Errors:   $ERRORS"
echo "  Warnings: $WARNINGS"
echo ""

if [ $ERRORS -eq 0 ]; then
    echo "Next steps:"
    echo "  1. Review warnings (if any)"
    echo "  2. Run: bash initialize_git_repo.sh"
    echo "  3. Follow instructions in GITHUB_SETUP.md"
fi

echo "=================================================="

# Exit with error if there are errors
if [ $ERRORS -gt 0 ]; then
    exit 1
fi

exit 0
