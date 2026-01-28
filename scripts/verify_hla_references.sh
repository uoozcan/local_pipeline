#!/bin/bash
# verify_hla_references.sh
# Comprehensive verification of HLA typing tool references

set -e

echo "======================================================================"
echo "   HLA TYPING PIPELINE - REFERENCE VERIFICATION"
echo "======================================================================"
echo ""
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

# Configuration
CONTAINER_DIR="/scratch/project_2008084/hla_references/singularity_cache/containers"
REF_DIR="/scratch/project_2008084/hla_references"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $2"
    else
        echo -e "${RED}✗${NC} $2"
    fi
}

# Function to print warning
print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

echo "======================================================================"
echo "1. CHECKING SINGULARITY CONTAINERS"
echo "======================================================================"
echo ""

# Check each container
for tool in optitype arcashla spechla; do
    SIF="${CONTAINER_DIR}/${tool}.sif"
    if [ -f "$SIF" ]; then
        SIZE=$(du -h "$SIF" | cut -f1)
        print_status 0 "${tool}.sif found (${SIZE})"
        
        # Test container execution
        if singularity exec "$SIF" echo "Container test OK" &>/dev/null; then
            echo "  → Container executable: YES"
        else
            print_warning "Container may have execution issues"
        fi
    else
        print_status 1 "${tool}.sif NOT FOUND at $SIF"
    fi
    echo ""
done

echo "======================================================================"
echo "2. CHECKING ARCASHLA REFERENCE"
echo "======================================================================"
echo ""

echo "Testing ArcasHLA container and reference..."
singularity exec ${CONTAINER_DIR}/arcashla.sif bash -c '
    echo "Container environment:"
    echo "  Python: $(python --version 2>&1)"
    echo "  ArcasHLA location: $(which arcasHLA 2>/dev/null || echo "NOT IN PATH")"
    echo ""
    
    # Check for reference directory
    for ref_path in "/usr/local/bin/dat/ref" "/home/arcasHLA-master/dat/ref" "/app/arcasHLA/dat/ref"; do
        if [ -d "$ref_path" ]; then
            echo "  ✓ Reference directory: $ref_path"
            ls -lh "$ref_path" 2>/dev/null | head -5
            
            # Check for Kallisto index
            if [ -f "$ref_path/hla.idx" ]; then
                echo "  ✓ Kallisto index found: $ref_path/hla.idx"
                ls -lh "$ref_path/hla.idx"
            else
                echo "  ✗ Kallisto index NOT found in $ref_path"
            fi
            echo ""
        fi
    done
    
    # Try to run arcasHLA help
    echo "Testing arcasHLA command:"
    arcasHLA --help 2>&1 | head -10 || echo "  ✗ arcasHLA command failed"
' || echo "  ✗ Cannot access ArcasHLA container"

echo ""

echo "======================================================================"
echo "3. CHECKING SPECHLA REFERENCE"
echo "======================================================================"
echo ""

echo "Testing SpecHLA container and database..."
singularity exec ${CONTAINER_DIR}/spechla.sif bash -c '
    echo "Container environment:"
    echo "  SpecHLA location: $(which SpecHLA.sh 2>/dev/null || echo "NOT IN PATH")"
    echo ""
    
    # Check for database
    for db_path in "/app/SpecHLA/database" "/opt/SpecHLA/database" "/usr/local/SpecHLA/database"; do
        if [ -d "$db_path" ]; then
            echo "  ✓ Database directory: $db_path"
            
            # Check for hg38 database
            if [ -d "$db_path/hg38" ]; then
                echo "  ✓ hg38 database found"
                echo "  Contents:"
                ls -lh "$db_path/hg38/" 2>/dev/null | head -10
            else
                echo "  ✗ hg38 database NOT found"
            fi
            
            # Check for hg19 database
            if [ -d "$db_path/hg19" ]; then
                echo "  ✓ hg19 database found"
            fi
            echo ""
        fi
    done
    
    # Check for required tools
    echo "Checking required tools:"
    for tool in bwa samtools bcftools freebayes; do
        if command -v $tool &>/dev/null; then
            echo "  ✓ $tool: $(which $tool)"
        else
            echo "  ✗ $tool: NOT FOUND"
        fi
    done
' || echo "  ✗ Cannot access SpecHLA container"

echo ""

echo "======================================================================"
echo "4. CHECKING OPTITYPE REFERENCE"
echo "======================================================================"
echo ""

echo "Testing OptiType container and reference..."
singularity exec ${CONTAINER_DIR}/optitype.sif bash -c '
    echo "Container environment:"
    echo "  OptiTypePipeline location: $(which OptiTypePipeline.py 2>/dev/null || echo "NOT IN PATH")"
    echo ""
    
    # Check for data directory
    for data_path in "/opt/OptiType/data" "/usr/local/bin/OptiType/data" "/app/OptiType/data"; do
        if [ -d "$data_path" ]; then
            echo "  ✓ Data directory: $data_path"
            ls -lh "$data_path/" 2>/dev/null | head -10
            echo ""
        fi
    done
    
    # Check for config file
    for config in "/opt/OptiType/config.ini" "/usr/local/bin/OptiType/config.ini"; do
        if [ -f "$config" ]; then
            echo "  ✓ Config file: $config"
            echo "  Contents:"
            cat "$config" 2>/dev/null | head -20
            echo ""
        fi
    done
' || echo "  ✗ Cannot access OptiType container"

echo ""

echo "======================================================================"
echo "5. TEST READ SIMULATION (if tools are working)"
echo "======================================================================"
echo ""

# Create a small test FASTQ
echo "Creating test FASTQ files..."
cat > test_R1.fastq << 'EOF'
@test_read_1
GCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
EOF

cat > test_R2.fastq << 'EOF'
@test_read_1
CGGCCTCGCTCTGGTTGTAGTAGCCGCGCAGGGTCCCAGGTTCGCTCGGTCAGTCTGTGAGTGGGCCTTCATATTCCGTGTCTCCTGGTCCCAATACTCCGGCCCCTCCTGCTCTATCCACGGCGCCCGCGGCTCCATCTTCTGGCTCGCGGCGTCGCTGTCGAACCGCACGAACTGCGTGTCGTCCACGTAGCCCACGGCGATGAAGCGGGGCTCCCCGCGGCCGGGCCGGGACACGGAGGTGTAGAAATACCTCATGGAGTGGGAGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
EOF

echo "Test files created."
echo ""

# Test ArcasHLA extraction
echo "Testing ArcasHLA extract with test reads..."
singularity exec ${CONTAINER_DIR}/arcashla.sif bash -c '
    arcasHLA extract test_R1.fastq test_R2.fastq -o test_arcas -v 2>&1 | head -30
' 2>&1 | tee arcashla_test.log || print_warning "ArcasHLA test failed"

echo ""

# Cleanup
rm -f test_R1.fastq test_R2.fastq

echo "======================================================================"
echo "6. SUMMARY"
echo "======================================================================"
echo ""

echo "Container Status:"
ls -lh ${CONTAINER_DIR}/*.sif 2>/dev/null || echo "No containers found"
echo ""

echo "Reference Check Complete!"
echo ""
echo "Next steps:"
echo "1. If any references are missing, rebuild them using tool-specific commands"
echo "2. If containers are missing, download or build them"
echo "3. Review log files for detailed error messages"
echo ""
echo "For detailed logs, check:"
echo "  - arcashla_test.log (if created)"
echo "  - Individual container outputs above"
echo ""
echo "======================================================================"
