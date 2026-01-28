/*
 * SpecHLA module (DSL2) - Supports both container and local installation modes.
 *
 * For local (HPC) mode without containers:
 *   - Set params.spechla_home to the local SpecHLA installation directory
 *   - Or set params.spechla_script to full path to SpecHLA.sh
 *   - Use profile 'local' or 'local_hpc' which disables containers for SpecHLA
 *
 * For container mode:
 *   - Use profile 'singularity' or 'docker'
 */

process SPECHLA_BAM {
    tag "${sample_id}"
    memory { params.max_memory ?: '14.GB' }
    cpus { params.max_cpus ?: 4 }
    publishDir "${params.outdir}/spechla", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}_spechla.tsv"), emit: results
    path("${sample_id}_spechla.log"), emit: logs

    script:
    def genes = (params.spechla_genes ?: params.hla_genes ?: 'A,B,C,DRB1,DQB1,DPB1').toString()
    def exon_only = (params.spechla_exon_only != null ? params.spechla_exon_only : 0) as Integer
    def ref_build = (params.reference_build ?: 'hg38').toString()
    def spechla_home   = (params.spechla_home   ?: '').toString()
    def spechla_script = (params.spechla_script ?: '').toString()
    def use_local = (params.spechla_use_local != null ? params.spechla_use_local : false)

    """
    set -euo pipefail

    SAMPLE="${sample_id}"
    BAM="${bam}"
    BAI="${bai}"
    GENES="${genes}"
    REF="${ref_build}"
    EXON_ONLY="${exon_only}"
    SPECHLA_HOME="${spechla_home}"
    SPECHLA_SCRIPT="${spechla_script}"
    USE_LOCAL="${use_local}"
    CPUS="${task.cpus}"

    LOG="\${SAMPLE}_spechla.log"
    echo "SpecHLA (BAM) sample=\${SAMPLE}" > "\$LOG"
    echo "BAM=\${BAM}" >> "\$LOG"
    [[ -n "\$BAI" ]] && echo "BAI=\${BAI}" >> "\$LOG" || true
    echo "GENES=\${GENES}" >> "\$LOG"
    echo "REF=\${REF}" >> "\$LOG"
    echo "EXON_ONLY=\${EXON_ONLY}" >> "\$LOG"
    echo "CPUS=\${CPUS}" >> "\$LOG"
    echo "USE_LOCAL=\${USE_LOCAL}" >> "\$LOG"
    echo "SPECHLA_HOME=\${SPECHLA_HOME}" >> "\$LOG"

    # For local mode, set up PATH early (before any tool usage)
    if [[ "\$USE_LOCAL" == "true" && -n "\$SPECHLA_HOME" && -d "\$SPECHLA_HOME/local_bin" ]]; then
      export PATH="\$SPECHLA_HOME/local_bin:\$PATH"
      echo "Added local_bin to PATH early" >> "\$LOG"
      if [[ -d "\$SPECHLA_HOME/local_tools/samtools-1.17/htslib-1.17" ]]; then
        export LD_LIBRARY_PATH="\$SPECHLA_HOME/local_tools/samtools-1.17/htslib-1.17:\${LD_LIBRARY_PATH:-}"
      fi
    fi

    # Ensure BAM index exists
    if [[ -n "\$BAI" && -s "\$BAI" ]]; then
      ln -sf "\$BAI" "\${BAM}.bai" 2>>"\$LOG" || true
    fi
    if [[ ! -s "\${BAM}.bai" && ! -s "\${BAM%.*}.bai" ]]; then
      echo "Index missing; creating with samtools index" >> "\$LOG"
      samtools index "\$BAM" >> "\$LOG" 2>&1
    fi

    # Resolve SpecHLA.sh - prefer local installation if available
    SPECHLA_SH="\$SPECHLA_SCRIPT"
    if [[ -z "\$SPECHLA_SH" && -n "\$SPECHLA_HOME" && -f "\$SPECHLA_HOME/script/whole/SpecHLA.sh" ]]; then
      SPECHLA_SH="\$SPECHLA_HOME/script/whole/SpecHLA.sh"
    fi
    if [[ -z "\$SPECHLA_SH" ]]; then
      for p in \\
        "/opt/SpecHLA/script/whole/SpecHLA.sh" \\
        "script/whole/SpecHLA.sh" \\
        "./script/whole/SpecHLA.sh" \\
        "./SpecHLA/script/whole/SpecHLA.sh" \\
        "/opt/SpecHLA/script/whole/SpecHLA.sh" \\
        "/usr/local/SpecHLA/script/whole/SpecHLA.sh" \\
        "/SpecHLA/script/whole/SpecHLA.sh"
      do
        [[ -f "\$p" ]] && SPECHLA_SH="\$p" && break
      done
    fi
    if [[ -z "\$SPECHLA_SH" ]] && command -v SpecHLA.sh &>/dev/null; then
      SPECHLA_SH="\$(command -v SpecHLA.sh)"
    fi

    if [[ -z "\$SPECHLA_SH" || ! -f "\$SPECHLA_SH" ]]; then
      echo "ERROR: SpecHLA.sh not found. Set params.spechla_home or params.spechla_script." >> "\$LOG"
      exit 127
    fi

    # Get SpecHLA base directory
    ORIG_SPECHLA_DIR="\$(dirname "\$(dirname "\$(dirname "\$SPECHLA_SH")")")"
    echo "Original SpecHLA directory: \$ORIG_SPECHLA_DIR" >> "\$LOG"

    # For local mode, use SpecHLA in place; for container mode, copy if read-only
    if [[ "\$USE_LOCAL" == "true" ]]; then
      # Local mode: use SpecHLA directory as-is
      SPECHLA_DIR="\$ORIG_SPECHLA_DIR"

      # Set PATH to include local tools
      if [[ -d "\$SPECHLA_DIR/local_bin" ]]; then
        export PATH="\$SPECHLA_DIR/local_bin:\$PATH"
        echo "Added local_bin to PATH" >> "\$LOG"
      fi

      # Set library paths for local tools
      if [[ -d "\$SPECHLA_DIR/local_tools/samtools-1.17/htslib-1.17" ]]; then
        export LD_LIBRARY_PATH="\$SPECHLA_DIR/local_tools/samtools-1.17/htslib-1.17:\${LD_LIBRARY_PATH:-}"
        echo "Set LD_LIBRARY_PATH for htslib" >> "\$LOG"
      fi
    else
      # Container mode: copy if read-only
      LOCAL_SPECHLA="\$PWD/SpecHLA_local"
      if [[ ! -w "\$ORIG_SPECHLA_DIR" ]]; then
        echo "Container is read-only, copying SpecHLA to local directory..." >> "\$LOG"
        rm -rf "\$LOCAL_SPECHLA"
        cp -r "\$ORIG_SPECHLA_DIR" "\$LOCAL_SPECHLA" 2>> "\$LOG"
        SPECHLA_DIR="\$LOCAL_SPECHLA"
        SPECHLA_SH="\$LOCAL_SPECHLA/script/whole/SpecHLA.sh"
      else
        SPECHLA_DIR="\$ORIG_SPECHLA_DIR"
      fi

      # Create spechla_env symlink to conda environment (container mode)
      if [[ -d "/opt/conda/envs/specHLA" ]]; then
        echo "Creating symlink for spechla_env to conda environment" >> "\$LOG"
        rm -f "\$SPECHLA_DIR/spechla_env" 2>/dev/null || true
        ln -sf /opt/conda/envs/specHLA "\$SPECHLA_DIR/spechla_env" 2>> "\$LOG"
        export PATH="/opt/conda/envs/specHLA/bin:\$PATH"
      fi
    fi

    echo "Working SpecHLA directory: \$SPECHLA_DIR" >> "\$LOG"
    echo "Using SpecHLA.sh=\${SPECHLA_SH}" >> "\$LOG"

    # Build bowtie2 indexes if they don't exist
    DB_REF="\$SPECHLA_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
    DB_REF_V2="\$SPECHLA_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"

    if [[ -f "\$DB_REF" && ! -f "\${DB_REF}.1.bt2" ]]; then
      echo "Building bowtie2 index for HLA database..." >> "\$LOG"
      bowtie2-build "\$DB_REF" "\$DB_REF" >> "\$LOG" 2>&1 || echo "Warning: bowtie2-build failed for DB_REF" >> "\$LOG"
    fi

    if [[ -f "\$DB_REF_V2" && ! -f "\${DB_REF_V2}.1.bt2" ]]; then
      echo "Building bowtie2 index for HLA database v2..." >> "\$LOG"
      bowtie2-build "\$DB_REF_V2" "\$DB_REF_V2" >> "\$LOG" 2>&1 || echo "Warning: bowtie2-build failed for DB_REF_V2" >> "\$LOG"
    fi

    # Build SpecHap binary if missing (only in container mode)
    SPECHAP_BIN="\$SPECHLA_DIR/bin/SpecHap/build/SpecHap"
    if [[ ! -f "\$SPECHAP_BIN" && "\$USE_LOCAL" != "true" ]]; then
      echo "SpecHap binary missing, attempting to compile..." >> "\$LOG"
      SPECHAP_SRC="\$SPECHLA_DIR/bin/SpecHap"
      if [[ -f "\$SPECHAP_SRC/CMakeLists.txt" ]]; then
        # Set environment for htslib from conda
        export CPLUS_INCLUDE_PATH="/opt/conda/envs/specHLA/include:\${CPLUS_INCLUDE_PATH:-}"
        export C_INCLUDE_PATH="/opt/conda/envs/specHLA/include:\${C_INCLUDE_PATH:-}"
        export LIBRARY_PATH="/opt/conda/envs/specHLA/lib:\${LIBRARY_PATH:-}"
        export LD_LIBRARY_PATH="/opt/conda/envs/specHLA/lib:\${LD_LIBRARY_PATH:-}"

        mkdir -p "\$SPECHAP_SRC/build" 2>> "\$LOG"
        cd "\$SPECHAP_SRC/build"
        cmake .. \\
          -DCMAKE_PREFIX_PATH="/opt/conda/envs/specHLA" \\
          -DCMAKE_CXX_FLAGS="-I/opt/conda/envs/specHLA/include" \\
          -DCMAKE_EXE_LINKER_FLAGS="-L/opt/conda/envs/specHLA/lib" \\
          -DHTSlib_INCLUDE_DIR="/opt/conda/envs/specHLA/include" \\
          -DHTSlib_LIBRARY="/opt/conda/envs/specHLA/lib/libhts.so" \\
          >> "\$LOG" 2>&1 && make -j\${CPUS:-4} VERBOSE=1 >> "\$LOG" 2>&1 || echo "Warning: SpecHap compilation failed" >> "\$LOG"
        cd - > /dev/null
        if [[ -f "\$SPECHAP_BIN" ]]; then
          echo "SpecHap compiled successfully" >> "\$LOG"
        else
          echo "Warning: SpecHap binary not found after compilation" >> "\$LOG"
        fi
      else
        echo "Warning: SpecHap CMakeLists.txt not found" >> "\$LOG"
      fi
    fi

    OUTDIR="spechla_out"
    rm -rf "\$OUTDIR" && mkdir -p "\$OUTDIR"

    R1="\$OUTDIR/\${SAMPLE}.R1.fq.gz"
    R2="\$OUTDIR/\${SAMPLE}.R2.fq.gz"

    echo "Converting BAM to paired FASTQs (compressed)" >> "\$LOG"
    # samtools fastq doesn't auto-compress in this container, so output to temp files then compress
    R1_TMP="\$OUTDIR/\${SAMPLE}.R1.fq"
    R2_TMP="\$OUTDIR/\${SAMPLE}.R2.fq"
    samtools fastq -1 "\$R1_TMP" -2 "\$R2_TMP" -0 /dev/null -s /dev/null -n "\$BAM" >> "\$LOG" 2>&1
    gzip -f "\$R1_TMP" >> "\$LOG" 2>&1
    gzip -f "\$R2_TMP" >> "\$LOG" 2>&1

    if [[ ! -s "\$R1" || ! -s "\$R2" ]]; then
      echo "ERROR: failed to generate FASTQs for SpecHLA" >> "\$LOG"
      exit 1
    fi

    echo "Running SpecHLA..." >> "\$LOG"
    # Note: SpecHLA -r is MAF (float), not reference build; -g is G-group annotation (0/1), not genes
    # SpecHLA types all standard HLA genes by default
    bash "\$SPECHLA_SH" \\
      -n "\$SAMPLE" \\
      -1 "\$R1" -2 "\$R2" \\
      -o "\$OUTDIR" \\
      -u "\$EXON_ONLY" \\
      -j "\$CPUS" >> "\$LOG" 2>&1

    # SpecHLA outputs to OUTDIR/SAMPLE/ subdirectory
    RESULT_FILE="\$OUTDIR/\${SAMPLE}/hla.result.txt"
    if [[ ! -s "\$RESULT_FILE" ]]; then
      # Fallback to direct path
      if [[ -s "\$OUTDIR/hla.result.txt" ]]; then
        RESULT_FILE="\$OUTDIR/hla.result.txt"
      else
        echo "ERROR: SpecHLA produced no hla.result.txt (empty or missing)" >> "\$LOG"
        ls -lah "\$OUTDIR" >> "\$LOG" 2>&1 || true
        ls -lah "\$OUTDIR/\${SAMPLE}/" >> "\$LOG" 2>&1 || true
        exit 1
      fi
    fi

    cp "\$RESULT_FILE" "\${SAMPLE}_spechla.tsv"
    echo "OK" >> "\$LOG"
    """
}

process SPECHLA_FASTQ {
    tag "${sample_id}"
    memory { params.max_memory ?: '14.GB' }
    cpus { params.max_cpus ?: 4 }
    publishDir "${params.outdir}/spechla", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path("${sample_id}_spechla.tsv"), emit: results
    path("${sample_id}_spechla.log"), emit: logs

    script:
    def genes = (params.spechla_genes ?: params.hla_genes ?: 'A,B,C,DRB1,DQB1,DPB1').toString()
    def exon_only = (params.spechla_exon_only != null ? params.spechla_exon_only : 0) as Integer
    def ref_build = (params.reference_build ?: 'hg38').toString()
    def spechla_home   = (params.spechla_home   ?: '').toString()
    def spechla_script = (params.spechla_script ?: '').toString()
    def use_local = (params.spechla_use_local != null ? params.spechla_use_local : false)

    """
    set -euo pipefail

    SAMPLE="${sample_id}"
    R1="${r1}"
    R2="${r2}"
    GENES="${genes}"
    REF="${ref_build}"
    EXON_ONLY="${exon_only}"
    SPECHLA_HOME="${spechla_home}"
    SPECHLA_SCRIPT="${spechla_script}"
    USE_LOCAL="${use_local}"
    CPUS="${task.cpus}"

    LOG="\${SAMPLE}_spechla.log"
    echo "SpecHLA (FASTQ) sample=\${SAMPLE}" > "\$LOG"
    echo "R1=\${R1}" >> "\$LOG"
    echo "R2=\${R2}" >> "\$LOG"
    echo "GENES=\${GENES}" >> "\$LOG"
    echo "REF=\${REF}" >> "\$LOG"
    echo "EXON_ONLY=\${EXON_ONLY}" >> "\$LOG"
    echo "CPUS=\${CPUS}" >> "\$LOG"
    echo "USE_LOCAL=\${USE_LOCAL}" >> "\$LOG"
    echo "SPECHLA_HOME=\${SPECHLA_HOME}" >> "\$LOG"

    SPECHLA_SH="\$SPECHLA_SCRIPT"
    if [[ -z "\$SPECHLA_SH" && -n "\$SPECHLA_HOME" && -f "\$SPECHLA_HOME/script/whole/SpecHLA.sh" ]]; then
      SPECHLA_SH="\$SPECHLA_HOME/script/whole/SpecHLA.sh"
    fi
    if [[ -z "\$SPECHLA_SH" ]]; then
      for p in \\
        "/opt/SpecHLA/script/whole/SpecHLA.sh" \\
        "script/whole/SpecHLA.sh" \\
        "./script/whole/SpecHLA.sh" \\
        "./SpecHLA/script/whole/SpecHLA.sh" \\
        "/opt/SpecHLA/script/whole/SpecHLA.sh" \\
        "/usr/local/SpecHLA/script/whole/SpecHLA.sh" \\
        "/SpecHLA/script/whole/SpecHLA.sh"
      do
        [[ -f "\$p" ]] && SPECHLA_SH="\$p" && break
      done
    fi
    if [[ -z "\$SPECHLA_SH" ]] && command -v SpecHLA.sh &>/dev/null; then
      SPECHLA_SH="\$(command -v SpecHLA.sh)"
    fi

    if [[ -z "\$SPECHLA_SH" || ! -f "\$SPECHLA_SH" ]]; then
      echo "ERROR: SpecHLA.sh not found. Set params.spechla_home or params.spechla_script." >> "\$LOG"
      exit 127
    fi

    # Get SpecHLA base directory
    ORIG_SPECHLA_DIR="\$(dirname "\$(dirname "\$(dirname "\$SPECHLA_SH")")")"
    echo "Original SpecHLA directory: \$ORIG_SPECHLA_DIR" >> "\$LOG"

    # For local mode, use SpecHLA in place; for container mode, copy if read-only
    if [[ "\$USE_LOCAL" == "true" ]]; then
      # Local mode: use SpecHLA directory as-is
      SPECHLA_DIR="\$ORIG_SPECHLA_DIR"

      # Set PATH to include local tools
      if [[ -d "\$SPECHLA_DIR/local_bin" ]]; then
        export PATH="\$SPECHLA_DIR/local_bin:\$PATH"
        echo "Added local_bin to PATH" >> "\$LOG"
      fi

      # Set library paths for local tools
      if [[ -d "\$SPECHLA_DIR/local_tools/samtools-1.17/htslib-1.17" ]]; then
        export LD_LIBRARY_PATH="\$SPECHLA_DIR/local_tools/samtools-1.17/htslib-1.17:\${LD_LIBRARY_PATH:-}"
        echo "Set LD_LIBRARY_PATH for htslib" >> "\$LOG"
      fi
    else
      # Container mode: copy if read-only
      LOCAL_SPECHLA="\$PWD/SpecHLA_local"
      if [[ ! -w "\$ORIG_SPECHLA_DIR" ]]; then
        echo "Container is read-only, copying SpecHLA to local directory..." >> "\$LOG"
        rm -rf "\$LOCAL_SPECHLA"
        cp -r "\$ORIG_SPECHLA_DIR" "\$LOCAL_SPECHLA" 2>> "\$LOG"
        SPECHLA_DIR="\$LOCAL_SPECHLA"
        SPECHLA_SH="\$LOCAL_SPECHLA/script/whole/SpecHLA.sh"
      else
        SPECHLA_DIR="\$ORIG_SPECHLA_DIR"
      fi

      # Create spechla_env symlink to conda environment (container mode)
      if [[ -d "/opt/conda/envs/specHLA" ]]; then
        echo "Creating symlink for spechla_env to conda environment" >> "\$LOG"
        rm -f "\$SPECHLA_DIR/spechla_env" 2>/dev/null || true
        ln -sf /opt/conda/envs/specHLA "\$SPECHLA_DIR/spechla_env" 2>> "\$LOG"
        export PATH="/opt/conda/envs/specHLA/bin:\$PATH"
      fi
    fi

    echo "Working SpecHLA directory: \$SPECHLA_DIR" >> "\$LOG"
    echo "Using SpecHLA.sh=\${SPECHLA_SH}" >> "\$LOG"

    # Build bowtie2 indexes if they don't exist
    DB_REF="\$SPECHLA_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.fasta"
    DB_REF_V2="\$SPECHLA_DIR/db/ref/hla_gen.format.filter.extend.DRB.no26789.v2.fasta"

    if [[ -f "\$DB_REF" && ! -f "\${DB_REF}.1.bt2" ]]; then
      echo "Building bowtie2 index for HLA database..." >> "\$LOG"
      bowtie2-build "\$DB_REF" "\$DB_REF" >> "\$LOG" 2>&1 || echo "Warning: bowtie2-build failed for DB_REF" >> "\$LOG"
    fi

    if [[ -f "\$DB_REF_V2" && ! -f "\${DB_REF_V2}.1.bt2" ]]; then
      echo "Building bowtie2 index for HLA database v2..." >> "\$LOG"
      bowtie2-build "\$DB_REF_V2" "\$DB_REF_V2" >> "\$LOG" 2>&1 || echo "Warning: bowtie2-build failed for DB_REF_V2" >> "\$LOG"
    fi

    # Build SpecHap binary if missing (only in container mode)
    SPECHAP_BIN="\$SPECHLA_DIR/bin/SpecHap/build/SpecHap"
    if [[ ! -f "\$SPECHAP_BIN" && "\$USE_LOCAL" != "true" ]]; then
      echo "SpecHap binary missing, attempting to compile..." >> "\$LOG"
      SPECHAP_SRC="\$SPECHLA_DIR/bin/SpecHap"
      if [[ -f "\$SPECHAP_SRC/CMakeLists.txt" ]]; then
        # Set environment for htslib from conda
        export CPLUS_INCLUDE_PATH="/opt/conda/envs/specHLA/include:\${CPLUS_INCLUDE_PATH:-}"
        export C_INCLUDE_PATH="/opt/conda/envs/specHLA/include:\${C_INCLUDE_PATH:-}"
        export LIBRARY_PATH="/opt/conda/envs/specHLA/lib:\${LIBRARY_PATH:-}"
        export LD_LIBRARY_PATH="/opt/conda/envs/specHLA/lib:\${LD_LIBRARY_PATH:-}"

        mkdir -p "\$SPECHAP_SRC/build" 2>> "\$LOG"
        cd "\$SPECHAP_SRC/build"
        cmake .. \\
          -DCMAKE_PREFIX_PATH="/opt/conda/envs/specHLA" \\
          -DCMAKE_CXX_FLAGS="-I/opt/conda/envs/specHLA/include" \\
          -DCMAKE_EXE_LINKER_FLAGS="-L/opt/conda/envs/specHLA/lib" \\
          -DHTSlib_INCLUDE_DIR="/opt/conda/envs/specHLA/include" \\
          -DHTSlib_LIBRARY="/opt/conda/envs/specHLA/lib/libhts.so" \\
          >> "\$LOG" 2>&1 && make -j\${CPUS:-4} VERBOSE=1 >> "\$LOG" 2>&1 || echo "Warning: SpecHap compilation failed" >> "\$LOG"
        cd - > /dev/null
        if [[ -f "\$SPECHAP_BIN" ]]; then
          echo "SpecHap compiled successfully" >> "\$LOG"
        else
          echo "Warning: SpecHap binary not found after compilation" >> "\$LOG"
        fi
      else
        echo "Warning: SpecHap CMakeLists.txt not found" >> "\$LOG"
      fi
    fi

    OUTDIR="spechla_out"
    rm -rf "\$OUTDIR" && mkdir -p "\$OUTDIR"

    echo "Running SpecHLA..." >> "\$LOG"
    # Note: SpecHLA -r is MAF (float), not reference build; -g is G-group annotation (0/1), not genes
    # SpecHLA types all standard HLA genes by default
    bash "\$SPECHLA_SH" \\
      -n "\$SAMPLE" \\
      -1 "\$R1" -2 "\$R2" \\
      -o "\$OUTDIR" \\
      -u "\$EXON_ONLY" \\
      -j "\$CPUS" >> "\$LOG" 2>&1

    # SpecHLA outputs to OUTDIR/SAMPLE/ subdirectory
    RESULT_FILE="\$OUTDIR/\${SAMPLE}/hla.result.txt"
    if [[ ! -s "\$RESULT_FILE" ]]; then
      # Fallback to direct path
      if [[ -s "\$OUTDIR/hla.result.txt" ]]; then
        RESULT_FILE="\$OUTDIR/hla.result.txt"
      else
        echo "ERROR: SpecHLA produced no hla.result.txt (empty or missing)" >> "\$LOG"
        ls -lah "\$OUTDIR" >> "\$LOG" 2>&1 || true
        ls -lah "\$OUTDIR/\${SAMPLE}/" >> "\$LOG" 2>&1 || true
        exit 1
      fi
    fi

    cp "\$RESULT_FILE" "\${SAMPLE}_spechla.tsv"
    echo "OK" >> "\$LOG"
    """
}
