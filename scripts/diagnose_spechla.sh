#!/bin/bash
# Diagnostic script for SpecHLA module syntax error

echo "=========================================="
echo "Diagnosing SpecHLA Module Syntax Error"
echo "=========================================="
echo ""

SPECHLA_FILE="/scratch/project_2008084/ozcanumu/hla_rnaseq_analysis/modules/spechla.nf"

if [ ! -f "$SPECHLA_FILE" ]; then
    echo "ERROR: File not found: $SPECHLA_FILE"
    exit 1
fi

echo "File found: $SPECHLA_FILE"
echo ""
echo "First 30 lines of the file:"
echo "----------------------------------------"
head -n 30 "$SPECHLA_FILE" | cat -n
echo "----------------------------------------"
echo ""

echo "Checking for common syntax issues:"
echo "1. Looking for unclosed braces/parentheses..."
echo "2. Checking process declarations..."
echo "3. Verifying directive syntax..."
echo ""

# Count braces
OPEN_BRACES=$(grep -o '{' "$SPECHLA_FILE" | wc -l)
CLOSE_BRACES=$(grep -o '}' "$SPECHLA_FILE" | wc -l)
echo "Open braces: $OPEN_BRACES"
echo "Close braces: $CLOSE_BRACES"
echo ""

# Show lines around line 7
echo "Context around line 7 (lines 1-15):"
echo "----------------------------------------"
sed -n '1,15p' "$SPECHLA_FILE" | cat -n
echo "----------------------------------------"
