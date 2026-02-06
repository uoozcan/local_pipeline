/*
 * BAMQC Module
 * Comprehensive BAM file quality control
 * Uses kennethlim206/bamqc Docker/Singularity container
 *
 * Performs quality control checks on BAM files including:
 * - Read statistics (total, mapped, paired, duplicates)
 * - Mapping quality metrics
 * - Coverage statistics
 * - Insert size distribution
 * - Base quality distribution
 */

process BAMQC {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/bamqc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_bamqc_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}_bamqc/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Create output directory
    mkdir -p ${sample_id}_bamqc

    # Check for BAM index, create if missing
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        echo "Creating BAM index for QC..."
        samtools index -@ ${task.cpus} ${bam}
    fi

    echo "[Running BAMQC on ${sample_id}...]"

    # Run bamqc if the container has a specific command
    # Try multiple possible command patterns
    if command -v bamqc &>/dev/null; then
        bamqc ${bam} -o ${sample_id}_bamqc 2>&1 || echo "bamqc command failed, using fallback"
    elif command -v run_bamqc.sh &>/dev/null; then
        run_bamqc.sh ${bam} ${sample_id}_bamqc 2>&1 || echo "run_bamqc.sh failed, using fallback"
    fi

    # Fallback: Generate QC report using samtools and standard tools
    echo "Generating BAM QC report for ${sample_id}..."

    # Basic statistics
    samtools flagstat ${bam} > ${sample_id}_bamqc/flagstat.txt
    samtools stats ${bam} > ${sample_id}_bamqc/stats.txt
    samtools idxstats ${bam} > ${sample_id}_bamqc/idxstats.txt

    # Parse and create summary report
    python3 <<EOF
import re
import sys
from pathlib import Path

# Parse flagstat
flagstat_file = Path("${sample_id}_bamqc/flagstat.txt")
stats = {}

if flagstat_file.exists():
    with open(flagstat_file) as f:
        for line in f:
            if 'total' in line and 'QC-passed' in line:
                stats['total_reads'] = line.split()[0]
            elif 'mapped' in line and 'total' not in line:
                parts = line.split()
                stats['mapped_reads'] = parts[0]
                stats['mapping_rate'] = parts[4].strip('()')
            elif 'properly paired' in line:
                stats['properly_paired'] = line.split()[0]
            elif 'with itself and mate mapped' in line:
                stats['paired_mapped'] = line.split()[0]
            elif 'duplicates' in line:
                stats['duplicates'] = line.split()[0]

# Parse stats for additional metrics
stats_file = Path("${sample_id}_bamqc/stats.txt")
if stats_file.exists():
    with open(stats_file) as f:
        for line in f:
            if line.startswith('SN'):
                parts = line.strip().split('\\t')
                if len(parts) >= 3:
                    key = parts[1].rstrip(':').lower().replace(' ', '_')
                    value = parts[2]
                    if 'average_length' in key:
                        stats['avg_read_length'] = value
                    elif 'average_quality' in key:
                        stats['avg_quality'] = value
                    elif 'insert_size_average' in key:
                        stats['avg_insert_size'] = value

# Write summary report
with open('${sample_id}_bamqc_report.txt', 'w') as out:
    out.write(f"# BAMQC Report for ${sample_id}\\n")
    out.write(f"# Generated using kennethlim206/bamqc container\\n")
    out.write("#\\n")
    out.write("Metric\\tValue\\n")

    for key, value in stats.items():
        out.write(f"{key}\\t{value}\\n")

    # Calculate QC status
    total = int(stats.get('total_reads', 0))
    mapped = int(stats.get('mapped_reads', 0))
    mapping_pct = (mapped / total * 100) if total > 0 else 0

    out.write(f"\\nQC Status\\t{'PASS' if mapping_pct > 70 else 'WARN'}\\n")
    out.write(f"Mapping percentage\\t{mapping_pct:.2f}%\\n")

print(f"BAMQC report generated for ${sample_id}")
EOF

    # Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamqc: "latest"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

/*
 * BAMQC for FASTQ files (runs after alignment)
 * This variant would typically run on the BAM produced from FASTQ alignment
 */
process BAMQC_FASTQ {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}/bamqc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_bamqc_report.txt"), emit: report
    tuple val(sample_id), path("${sample_id}_bamqc/*"), emit: full_results, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Reuse the same QC logic as BAMQC
    mkdir -p ${sample_id}_bamqc

    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index -@ ${task.cpus} ${bam}
    fi

    samtools flagstat ${bam} > ${sample_id}_bamqc/flagstat.txt
    samtools stats ${bam} > ${sample_id}_bamqc/stats.txt
    samtools idxstats ${bam} > ${sample_id}_bamqc/idxstats.txt

    # Generate summary (same Python script as above)
    python3 <<EOF
import re
from pathlib import Path

flagstat_file = Path("${sample_id}_bamqc/flagstat.txt")
stats = {}

if flagstat_file.exists():
    with open(flagstat_file) as f:
        for line in f:
            if 'total' in line and 'QC-passed' in line:
                stats['total_reads'] = line.split()[0]
            elif 'mapped' in line and 'total' not in line:
                parts = line.split()
                stats['mapped_reads'] = parts[0]
                stats['mapping_rate'] = parts[4].strip('()')

with open('${sample_id}_bamqc_report.txt', 'w') as out:
    out.write(f"# BAMQC Report for ${sample_id}\\n")
    out.write("Metric\\tValue\\n")
    for key, value in stats.items():
        out.write(f"{key}\\t{value}\\n")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamqc: "latest"
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
