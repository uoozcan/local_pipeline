#!/bin/bash
#
# Script to download and set up databases for HLA typing tools
# Run this once before using HLAscan or HLA*LA
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DB_DIR="${PROJECT_DIR}/hla_references/databases"

mkdir -p "$DB_DIR"

echo "=============================================="
echo "HLA Typing Tools - Database Setup"
echo "=============================================="
echo "Database directory: $DB_DIR"
echo ""

# ============================================
# HLAscan - IMGT/HLA Database
# ============================================
setup_hlascan_db() {
    echo "[1/2] Setting up HLAscan IMGT/HLA database..."

    HLASCAN_DB="${DB_DIR}/hlascan_db"
    mkdir -p "$HLASCAN_DB"

    if [ -d "$HLASCAN_DB/IMGT_HLA" ] && [ "$(ls -A $HLASCAN_DB/IMGT_HLA 2>/dev/null)" ]; then
        echo "  HLAscan database already exists at $HLASCAN_DB"
        return 0
    fi

    echo "  Downloading IMGT/HLA database from GitHub..."
    cd "$HLASCAN_DB"

    # Clone IMGT/HLA database
    if [ ! -d "IMGTHLA" ]; then
        git clone --depth 1 https://github.com/ANHIG/IMGTHLA.git
    fi

    # Create the expected structure
    mkdir -p IMGT_HLA

    # Copy required files
    if [ -d "IMGTHLA/fasta" ]; then
        cp -r IMGTHLA/fasta/* IMGT_HLA/ 2>/dev/null || true
    fi
    if [ -d "IMGTHLA/alignments" ]; then
        cp -r IMGTHLA/alignments/* IMGT_HLA/ 2>/dev/null || true
    fi

    echo "  HLAscan database setup complete: $HLASCAN_DB/IMGT_HLA"
}

# ============================================
# HLA*LA - Graph Database
# ============================================
setup_hlala_db() {
    echo ""
    echo "[2/2] Setting up HLA*LA graph database..."

    HLALA_DB="${DB_DIR}/hlala_graphs"
    mkdir -p "$HLALA_DB"

    if [ -d "$HLALA_DB/PRG_MHC_GRCh38_withIMGT" ]; then
        echo "  HLA*LA graph already exists at $HLALA_DB"
        return 0
    fi

    echo "  Downloading HLA*LA PRG_MHC_GRCh38_withIMGT graph..."
    echo "  (This is a large download ~10GB, may take some time)"

    cd "$HLALA_DB"

    # Download from Zenodo or HLA*LA repository
    # The graph files are hosted at: http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz

    GRAPH_URL="http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz"
    GRAPH_FILE="PRG_MHC_GRCh38_withIMGT.tar.gz"

    if [ ! -f "$GRAPH_FILE" ]; then
        echo "  Downloading graph from $GRAPH_URL"
        wget -c "$GRAPH_URL" -O "$GRAPH_FILE" || {
            echo "  ERROR: Failed to download HLA*LA graph."
            echo "  Please download manually from:"
            echo "    $GRAPH_URL"
            echo "  And extract to: $HLALA_DB/"
            return 1
        }
    fi

    echo "  Extracting graph..."
    tar -xzf "$GRAPH_FILE"

    echo "  HLA*LA graph setup complete: $HLALA_DB/PRG_MHC_GRCh38_withIMGT"
}

# Run setup
echo ""
setup_hlascan_db

echo ""
setup_hlala_db

echo ""
echo "=============================================="
echo "Database Setup Complete!"
echo "=============================================="
echo ""
echo "Database locations:"
echo "  HLAscan: ${DB_DIR}/hlascan_db/IMGT_HLA"
echo "  HLA*LA:  ${DB_DIR}/hlala_graphs/PRG_MHC_GRCh38_withIMGT"
echo ""
echo "Now update your wrapper scripts to use these databases."
echo "=============================================="
