process AGGREGATE_RESULTS {
    tag "$sample_id"
    label 'low_resources'
    publishDir "${params.outdir}/${sample_id}/consensus", mode: 'copy'
    
    input:
    tuple val(sample_id), path(results)
    
    output:
    tuple val(sample_id), path("${sample_id}.consensus.tsv"), emit: consensus
    path "${sample_id}.summary.txt", emit: summary
    path "${sample_id}.aggregation.log", emit: log
    
    script:
    """
    #!/usr/bin/env python3
    import json, os
    from collections import defaultdict, Counter
    
    hla_calls = defaultdict(list)
    tools_used = []
    log = open("${sample_id}.aggregation.log", "w")
    
    result_files = "${results}".split()
    for result_file in result_files:
        if not os.path.exists(result_file):
            continue
        
        if "optitype" in result_file:
            tool_name = "OptiType"
            try:
                with open(result_file, 'r') as f:
                    next(f)
                    for line in f:
                        parts = line.strip().split('\\t')
                        if len(parts) >= 6:
                            for i, gene in enumerate(['A', 'A', 'B', 'B', 'C', 'C']):
                                if i < len(parts) and parts[i]:
                                    hla_calls[gene].append((tool_name, parts[i]))
                tools_used.append(tool_name)
            except: pass
    
    consensus = {}
    for gene, calls in hla_calls.items():
        allele_counts = Counter([allele for tool, allele in calls])
        if allele_counts:
            top_allele, top_count = allele_counts.most_common()[0]
            consensus[gene] = {'allele': top_allele, 'support': top_count, 'total': len(calls)}
    
    with open("${sample_id}.consensus.tsv", 'w') as f:
        f.write("Gene\\tAllele\\tSupport\\tTotal_Calls\\n")
        for gene in sorted(consensus.keys()):
            info = consensus[gene]
            f.write(f"{gene}\\t{info['allele']}\\t{info['support']}\\t{info['total']}\\n")
    
    with open("${sample_id}.summary.txt", 'w') as f:
        f.write(f"Sample: ${sample_id}\\nTools: {', '.join(set(tools_used))}\\nGenes: {len(consensus)}\\n")
    
    log.close()
    """
    
    stub:
    """
    touch ${sample_id}.consensus.tsv ${sample_id}.summary.txt ${sample_id}.aggregation.log
    """
}
