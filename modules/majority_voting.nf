// Nextflow module for HLA majority voting
// Add this to your pipeline's modules directory

process MAJORITY_VOTING {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/majority_voting", mode: 'copy'
    
    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.5.3 conda-forge::numpy=1.23.5" : null)
    
    input:
    tuple val(sample_id), 
          path(optitype_result, stageAs: 'optitype/*'), 
          path(arcashla_result, stageAs: 'arcashla/*'), 
          path(spechla_result, stageAs: 'spechla/*'),
          path(hlahd_result, stageAs: 'hlahd/*'),
          path(hlala_result, stageAs: 'hlala/*')
    
    output:
    tuple val(sample_id), path("${sample_id}.consensus.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}.consensus.json"), emit: json
    path("${sample_id}.majority_voting.log"), emit: log
    
    script:
    def optitype_arg = optitype_result.name != 'NO_FILE' ? "--optitype optitype/${optitype_result.name}" : ""
    def arcashla_arg = arcashla_result.name != 'NO_FILE' ? "--arcashla arcashla/${arcashla_result.name}" : ""
    def spechla_arg = spechla_result.name != 'NO_FILE' ? "--spechla spechla/${spechla_result.name}" : ""
    def hlahd_arg = hlahd_result.name != 'NO_FILE' ? "--hlahd hlahd/${hlahd_result.name}" : ""
    def hlala_arg = hlala_result.name != 'NO_FILE' ? "--hlala hlala/${hlala_result.name}" : ""
    
    """
    echo "Starting majority voting for sample: ${sample_id}" > ${sample_id}.majority_voting.log
    echo "Tools provided:" >> ${sample_id}.majority_voting.log
    ls -lh optitype/ arcashla/ spechla/ hlahd/ hlala/ >> ${sample_id}.majority_voting.log 2>&1
    
    python3 ${projectDir}/bin/majority_voting.py \\
        --sample-id ${sample_id} \\
        ${optitype_arg} \\
        ${arcashla_arg} \\
        ${spechla_arg} \\
        ${hlahd_arg} \\
        ${hlala_arg} \\
        --output ${sample_id}.consensus.tsv \\
        --json-output ${sample_id}.consensus.json \\
        --genes A,B,C,DQA1,DQB1,DRB1 \\
        2>&1 | tee -a ${sample_id}.majority_voting.log
    
    echo "" >> ${sample_id}.majority_voting.log
    echo "Majority voting completed successfully" >> ${sample_id}.majority_voting.log
    """
}

process AGGREGATE_MAJORITY_VOTING {
    publishDir "${params.outdir}/majority_voting", mode: 'copy'
    
    input:
    path(consensus_files)
    
    output:
    path("all_samples.consensus.tsv"), emit: combined_tsv
    path("consensus_summary.txt"), emit: summary
    
    script:
    """
    # Combine all consensus TSV files
    head -1 ${consensus_files[0]} > all_samples.consensus.tsv
    for file in ${consensus_files}; do
        tail -n +2 \$file >> all_samples.consensus.tsv
    done
    
    # Create summary
    echo "HLA Typing Majority Voting Summary" > consensus_summary.txt
    echo "===================================" >> consensus_summary.txt
    echo "" >> consensus_summary.txt
    echo "Total samples: \$(tail -n +2 all_samples.consensus.tsv | cut -f1 | sort -u | wc -l)" >> consensus_summary.txt
    echo "" >> consensus_summary.txt
    echo "Typing success rate by gene:" >> consensus_summary.txt
    for gene in A B C DQA1 DQB1 DRB1; do
        total=\$(grep -w "\$gene" all_samples.consensus.tsv | wc -l)
        typed=\$(grep -w "\$gene" all_samples.consensus.tsv | awk '\$3 != "None"' | wc -l)
        if [ \$total -gt 0 ]; then
            rate=\$(echo "scale=2; \$typed * 100 / \$total" | bc)
            echo "  \$gene: \$typed/\$total (\${rate}%)" >> consensus_summary.txt
        fi
    done
    
    echo "" >> consensus_summary.txt
    echo "Tool usage frequency:" >> consensus_summary.txt
    tail -n +2 all_samples.consensus.tsv | cut -f5 | sort | uniq -c | sort -rn >> consensus_summary.txt
    
    cat consensus_summary.txt
    """
}

// Workflow to integrate majority voting
workflow MAJORITY_VOTING_WORKFLOW {
    take:
    optitype_results   // channel: [ val(sample_id), path(result_file) ]
    arcashla_results   // channel: [ val(sample_id), path(result_file) ]
    spechla_results    // channel: [ val(sample_id), path(result_file) ]
    hlahd_results      // channel: [ val(sample_id), path(result_file) ]
    hlala_results      // channel: [ val(sample_id), path(result_file) ]
    
    main:
    // Combine all results by sample_id
    // Use join to merge channels, providing default empty files for missing tools
    
    def empty_file = file("${projectDir}/assets/NO_FILE")
    
    combined = optitype_results
        .map { sample, file -> tuple(sample, file) }
        .join(arcashla_results.map { sample, file -> tuple(sample, file) }, remainder: true)
        .join(spechla_results.map { sample, file -> tuple(sample, file) }, remainder: true)
        .join(hlahd_results.map { sample, file -> tuple(sample, file) }, remainder: true)
        .join(hlala_results.map { sample, file -> tuple(sample, file) }, remainder: true)
        .map { row ->
            def sample = row[0]
            def opti = row[1] ?: empty_file
            def arcas = row[2] ?: empty_file
            def spec = row[3] ?: empty_file
            def hd = row[4] ?: empty_file
            def la = row[5] ?: empty_file
            tuple(sample, opti, arcas, spec, hd, la)
        }
    
    // Run majority voting
    MAJORITY_VOTING(combined)
    
    // Aggregate results
    AGGREGATE_MAJORITY_VOTING(
        MAJORITY_VOTING.out.tsv.map { it[1] }.collect()
    )
    
    emit:
    consensus_tsv = MAJORITY_VOTING.out.tsv
    consensus_json = MAJORITY_VOTING.out.json
    combined_tsv = AGGREGATE_MAJORITY_VOTING.out.combined_tsv
    summary = AGGREGATE_MAJORITY_VOTING.out.summary
}
