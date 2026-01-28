// modules/optitype.nf
// Final working OptiType module - handles actual output structure

process OPTITYPE {
    tag "$sample_id"
    label 'process_medium'
    container '/scratch/project_2008084/hla_references/singularity_cache/containers/optitype.sif'
    
    publishDir "${params.outdir}/${sample_id}/optitype", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_result.tsv"), emit: results
    path "${sample_id}_coverage_plot.pdf", emit: plot, optional: true
    path "${sample_id}.optitype.log", emit: log
    
    script:
    def seq_type = params.optitype_seq_type ?: 'dna'
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "==================================================" > ${sample_id}.optitype.log
    echo "OptiType HLA Typing - ${sample_id}" >> ${sample_id}.optitype.log
    echo "==================================================" >> ${sample_id}.optitype.log
    echo "Started: \$(date -Iseconds)" >> ${sample_id}.optitype.log
    echo "Sequence type: ${seq_type}" >> ${sample_id}.optitype.log
    echo "R1: ${read1}" >> ${sample_id}.optitype.log
    echo "R2: ${read2}" >> ${sample_id}.optitype.log
    echo "CPUs: ${task.cpus}" >> ${sample_id}.optitype.log
    echo "Memory: ${task.memory}" >> ${sample_id}.optitype.log
    echo "" >> ${sample_id}.optitype.log
    
    # Verify input files
    if [ ! -s ${read1} ]; then
        echo "ERROR: R1 file is missing or empty" >> ${sample_id}.optitype.log
        exit 1
    fi
    
    if [ ! -s ${read2} ]; then
        echo "ERROR: R2 file is missing or empty" >> ${sample_id}.optitype.log
        exit 1
    fi
    
    # Check read counts
    R1_READS=\$(zcat ${read1} 2>/dev/null | wc -l | awk '{print int(\$1/4)}' || echo "0")
    R2_READS=\$(zcat ${read2} 2>/dev/null | wc -l | awk '{print int(\$1/4)}' || echo "0")
    
    echo "Read counts:" >> ${sample_id}.optitype.log
    echo "  R1: \$R1_READS reads" >> ${sample_id}.optitype.log
    echo "  R2: \$R2_READS reads" >> ${sample_id}.optitype.log
    echo "" >> ${sample_id}.optitype.log
    
    if [ \$R1_READS -eq 0 ] || [ \$R2_READS -eq 0 ]; then
        echo "ERROR: No reads found in input files" >> ${sample_id}.optitype.log
        exit 1
    fi
    
    # Run OptiType
    echo "Running OptiType..." >> ${sample_id}.optitype.log
    
    # Don't exit on error immediately - we need to check output
    set +e
    OptiTypePipeline.py \\
        --input ${read1} ${read2} \\
        --${seq_type} \\
        --verbose \\
        --outdir optitype_out \\
        --prefix ${sample_id} \\
        2>&1 | tee -a ${sample_id}.optitype.log
    
    OPTITYPE_EXIT=\$?
    set -e
    
    echo "" >> ${sample_id}.optitype.log
    echo "OptiType process completed with exit code: \$OPTITYPE_EXIT" >> ${sample_id}.optitype.log
    echo "" >> ${sample_id}.optitype.log
    
    # Show directory structure for debugging
    echo "Output directory structure:" >> ${sample_id}.optitype.log
    if [ -d "optitype_out" ]; then
        find optitype_out -type f >> ${sample_id}.optitype.log
    else
        echo "  WARNING: optitype_out directory not found" >> ${sample_id}.optitype.log
    fi
    echo "" >> ${sample_id}.optitype.log
    
    # Find result file - OptiType creates it directly in optitype_out/
    # Pattern: optitype_out/SAMPLENAME_result.tsv
    RESULT_FILE=""
    
    # Check multiple possible locations
    if [ -f "optitype_out/${sample_id}_result.tsv" ]; then
        RESULT_FILE="optitype_out/${sample_id}_result.tsv"
    elif [ -f "optitype_out/${sample_id}/${sample_id}_result.tsv" ]; then
        RESULT_FILE="optitype_out/${sample_id}/${sample_id}_result.tsv"
    else
        # Search for any result.tsv file
        RESULT_FILE=\$(find optitype_out -name "*_result.tsv" -o -name "result.tsv" 2>/dev/null | head -1)
    fi
    
    if [ -n "\$RESULT_FILE" ] && [ -f "\$RESULT_FILE" ]; then
        echo "Found result file: \$RESULT_FILE" >> ${sample_id}.optitype.log
        cp "\$RESULT_FILE" ${sample_id}_result.tsv
        
        # Verify it has content
        if [ -s ${sample_id}_result.tsv ]; then
            LINES=\$(wc -l < ${sample_id}_result.tsv)
            echo "Result file copied successfully (\$LINES lines)" >> ${sample_id}.optitype.log
            echo "" >> ${sample_id}.optitype.log
            echo "Result preview:" >> ${sample_id}.optitype.log
            head -3 ${sample_id}_result.tsv >> ${sample_id}.optitype.log
            echo "" >> ${sample_id}.optitype.log
        else
            echo "ERROR: Result file is empty after copy" >> ${sample_id}.optitype.log
            exit 1
        fi
    else
        echo "ERROR: Could not find result file" >> ${sample_id}.optitype.log
        echo "Searched locations:" >> ${sample_id}.optitype.log
        echo "  - optitype_out/${sample_id}_result.tsv" >> ${sample_id}.optitype.log
        echo "  - optitype_out/${sample_id}/${sample_id}_result.tsv" >> ${sample_id}.optitype.log
        echo "  - optitype_out/*_result.tsv" >> ${sample_id}.optitype.log
        echo "" >> ${sample_id}.optitype.log
        echo "Available files:" >> ${sample_id}.optitype.log
        find . -name "*.tsv" >> ${sample_id}.optitype.log 2>&1
        exit 1
    fi
    
    # Try to find coverage plot
    PLOT_FILE=""
    if [ -f "optitype_out/${sample_id}_coverage_plot.pdf" ]; then
        PLOT_FILE="optitype_out/${sample_id}_coverage_plot.pdf"
    elif [ -f "optitype_out/${sample_id}/${sample_id}_coverage_plot.pdf" ]; then
        PLOT_FILE="optitype_out/${sample_id}/${sample_id}_coverage_plot.pdf"
    else
        PLOT_FILE=\$(find optitype_out -name "*_coverage_plot.pdf" 2>/dev/null | head -1)
    fi
    
    if [ -n "\$PLOT_FILE" ] && [ -f "\$PLOT_FILE" ]; then
        cp "\$PLOT_FILE" ${sample_id}_coverage_plot.pdf
        echo "Coverage plot copied: ${sample_id}_coverage_plot.pdf" >> ${sample_id}.optitype.log
    else
        echo "No coverage plot found (optional)" >> ${sample_id}.optitype.log
    fi
    
    echo "" >> ${sample_id}.optitype.log
    echo "==================================================" >> ${sample_id}.optitype.log
    echo "OptiType completed successfully!" >> ${sample_id}.optitype.log
    echo "Completed: \$(date -Iseconds)" >> ${sample_id}.optitype.log
    echo "==================================================" >> ${sample_id}.optitype.log
    """
}
