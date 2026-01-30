#!/usr/bin/env nextflow

/*
 * HLA Typing Pipeline
 * Multi-tool HLA typing from BAM/FASTQ with weighted consensus voting
 * Supports both BAM and paired FASTQ inputs
 * Includes FastQC, visualizations, and MultiQC integration
 *
 * Authors: HLA Analysis Team
 * Version: 1.2.0
 */

nextflow.enable.dsl = 2

// Import modules - BAM versions
include { SPECHLA } from './modules/spechla'
include { HLAHD } from './modules/hlahd'
include { HLALA } from './modules/hlala'
include { ARCASHLA } from './modules/arcashla'
include { OPTITYPE } from './modules/optitype'

// Import modules - FASTQ versions
include { SPECHLA_FASTQ } from './modules/spechla'
include { HLAHD_FASTQ } from './modules/hlahd'
include { ARCASHLA_FASTQ } from './modules/arcashla'
include { OPTITYPE_FASTQ } from './modules/optitype'

// Import QC and consensus modules
include { QC_BAM; QC_FASTQ } from './modules/qc'
include { CONSENSUS } from './modules/consensus'

// Import FastQC modules
include { FASTQC_BAM; FASTQC_FASTQ } from './modules/fastqc'

// Import visualization modules
include { HLA_VISUALIZE; HLA_SUMMARY_REPORT } from './modules/visualize'

// Import MultiQC module
include { MULTIQC } from './modules/multiqc'

// Import LOH module
include { HLA_LOH; HLA_LOH_VISUALIZE; HLA_LOH_SUMMARY } from './modules/loh'

// Help message
def helpMessage() {
    log.info """
    ===========================================
    HLA Typing Pipeline  v${workflow.manifest.version}
    ===========================================

    Usage:
        nextflow run main.nf --input_bam sample.bam --outdir results
        nextflow run main.nf --input_fastq_1 R1.fq.gz --input_fastq_2 R2.fq.gz --outdir results

    Input options (choose one):
        --input_bam           Path to input BAM file (single sample)
        --input_fastq_1       Path to R1 FASTQ file (single sample)
        --input_fastq_2       Path to R2 FASTQ file (single sample)
        --input_samplesheet   Path to samplesheet CSV (multiple samples)

    Samplesheet format (CSV with header):
        For BAM:   sample_id,bam_path
        For FASTQ: sample_id,fastq_1,fastq_2

    Optional arguments:
        --outdir            Output directory (default: ./results)
        --reference         Reference genome: hg38 or hg19 (default: hg38)
        --tools             HLA typing tools to use (default: spechla,hlahd)
                            Options: spechla,hlahd,hlala,arcashla,optitype
                            Note: hlala only supports BAM input
        --seq_type          Sequence type for OptiType: dna or rna (default: dna)
        --hlala_graph       HLA*LA graph (default: PRG_MHC_GRCh38_withIMGT)
        --hla_genes         HLA genes to type (default: classical HLA genes)
        --resolution        Output resolution: 2-field or 4-field (default: 2-field)
        --min_tools         Minimum tools for consensus (default: 1)
        --max_cpus          Maximum CPUs per process (default: 8)
        --max_memory        Maximum memory per process (default: 32.GB)

    QC thresholds:
        --min_hla_reads     Minimum HLA reads for warning (default: 1000)
        --min_read_length   Minimum average read length (default: 50)

    Consensus weighting:
        --expected_reads    Expected reads per allele for confidence (default: 1000)
        --weighting         Weighting method: equal, read_confidence, tool_quality
                            (default: read_confidence)

    LOH (Loss of Heterozygosity) analysis:
        --run_loh           Enable LOH analysis (requires SpecHLA, default: false)
        --tumor_purity      Tumor purity estimate (0-1), required for LOH
        --tumor_ploidy      Tumor ploidy estimate, required for LOH
        --loh_het_cutoff    Minimum het SNPs for LOH call (default: 5)

    Output reports:
        - Per-sample: QC report, consensus HLA types, visualizations
        - Summary: Multi-sample report with allele frequencies
        - MultiQC: Aggregated FastQC and HLA statistics

    Profiles:
        -profile singularity    Use Singularity containers
        -profile docker         Use Docker containers
        -profile slurm          Run on SLURM cluster
        -profile test           Run with test parameters

    Examples:
        # Single BAM file
        nextflow run main.nf --input_bam sample.bam -profile singularity

        # Single paired FASTQ files
        nextflow run main.nf --input_fastq_1 R1.fq.gz --input_fastq_2 R2.fq.gz -profile singularity

        # Multiple samples with samplesheet (BAM)
        nextflow run main.nf --input_samplesheet samples_bam.csv -profile singularity

        # Multiple samples with samplesheet (FASTQ)
        nextflow run main.nf --input_samplesheet samples_fastq.csv -profile singularity

        # Use multiple tools
        nextflow run main.nf --input_bam sample.bam --tools hlahd,arcashla -profile singularity

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Determine input type
def input_type = null
if (params.input_bam) {
    input_type = 'bam'
} else if (params.input_fastq_1 && params.input_fastq_2) {
    input_type = 'fastq'
} else if (params.input_samplesheet) {
    // Determine from samplesheet header
    def samplesheet = file(params.input_samplesheet)
    def header = samplesheet.readLines()[0]
    if (header.contains('fastq_1') || header.contains('fastq1') || header.contains('fq1')) {
        input_type = 'fastq'
    } else {
        input_type = 'bam'
    }
} else {
    log.error "Please provide input: --input_bam, --input_fastq_1/--input_fastq_2, or --input_samplesheet"
    helpMessage()
    exit 1
}

// Parse tools
def tools_list = params.tools.tokenize(',')

// Warn if hlala is selected with FASTQ input
if (input_type == 'fastq' && 'hlala' in tools_list) {
    log.warn "HLA*LA (hlala) only supports BAM input. It will be skipped for FASTQ samples."
    tools_list = tools_list.findAll { it != 'hlala' }
}

// Log parameters
log.info """
===========================================
HLA Typing Pipeline  v${workflow.manifest.version}
===========================================
Input type      : ${input_type.toUpperCase()}
Input           : ${params.input_bam ?: params.input_fastq_1 ?: params.input_samplesheet}
Output          : ${params.outdir}
Reference       : ${params.reference}
Tools           : ${tools_list.join(', ')}
HLA genes       : ${params.hla_genes}
Resolution      : ${params.resolution}
Min HLA reads   : ${params.min_hla_reads}
Max CPUs        : ${params.max_cpus}
Max Memory      : ${params.max_memory}
===========================================
"""

// Create input channel for BAM
def create_bam_channel() {
    if (params.input_bam) {
        def bam_file = file(params.input_bam)
        def sample_id = bam_file.baseName.replaceAll(/\.sorted$|\.dedup$|\.bam$/, '')
        return Channel.of([sample_id, bam_file])
    } else {
        return Channel
            .fromPath(params.input_samplesheet)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, file(row.bam_path)] }
    }
}

// Create input channel for FASTQ
def create_fastq_channel() {
    if (params.input_fastq_1) {
        def fq1 = file(params.input_fastq_1)
        def fq2 = file(params.input_fastq_2)
        def sample_id = fq1.baseName.replaceAll(/\.fastq$|\.fq$|\.gz$|_R?1$|_1$/, '').replaceAll(/\.fastq$|\.fq$/, '')
        return Channel.of([sample_id, fq1, fq2])
    } else {
        return Channel
            .fromPath(params.input_samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def fq1_col = row.fastq_1 ?: row.fastq1 ?: row.fq1 ?: row.read1 ?: row.R1
                def fq2_col = row.fastq_2 ?: row.fastq2 ?: row.fq2 ?: row.read2 ?: row.R2
                [row.sample_id, file(fq1_col), file(fq2_col)]
            }
    }
}

// Main workflow
workflow {
    // Results channel
    ch_results = Channel.empty()

    // FastQC results channel
    ch_fastqc = Channel.empty()

    // QC reports channel
    ch_qc_reports = Channel.empty()

    if (input_type == 'bam') {
        // ===== BAM INPUT WORKFLOW =====
        ch_bam = create_bam_channel()

        // Run QC on BAM files
        QC_BAM(
            ch_bam,
            params.reference,
            params.min_hla_reads,
            params.min_read_length
        )

        // Use validated BAM for downstream
        ch_input = QC_BAM.out.validated_bam
        ch_qc_reports = ch_qc_reports.mix(QC_BAM.out.qc_report)

        // Run FastQC on BAM files
        FASTQC_BAM(ch_bam)
        ch_fastqc = ch_fastqc.mix(FASTQC_BAM.out.zip.map { sample_id, zip -> zip })

        // Run SpecHLA if requested
        if ('spechla' in tools_list) {
            SPECHLA(ch_input, params.reference)
            ch_results = ch_results.mix(SPECHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'spechla', result_file]
            })
        }

        // Run HLA-HD if requested
        if ('hlahd' in tools_list) {
            HLAHD(ch_input, params.reference, params.hla_genes)
            ch_results = ch_results.mix(HLAHD.out.results.map { sample_id, result_file ->
                [sample_id, 'hlahd', result_file]
            })
        }

        // Run HLA*LA if requested
        if ('hlala' in tools_list) {
            HLALA(ch_input, params.hlala_graph)
            ch_results = ch_results.mix(HLALA.out.results.map { sample_id, result_file ->
                [sample_id, 'hlala', result_file]
            })
        }

        // Run arcasHLA if requested
        if ('arcashla' in tools_list) {
            ARCASHLA(ch_input, params.reference)
            ch_results = ch_results.mix(ARCASHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'arcashla', result_file]
            })
        }

        // Run OptiType if requested
        if ('optitype' in tools_list) {
            OPTITYPE(ch_input)
            ch_results = ch_results.mix(OPTITYPE.out.results.map { sample_id, result_file ->
                [sample_id, 'optitype', result_file]
            })
        }

    } else {
        // ===== FASTQ INPUT WORKFLOW =====
        ch_fastq = create_fastq_channel()

        // Run QC on FASTQ files
        QC_FASTQ(
            ch_fastq,
            params.min_hla_reads,
            params.min_read_length
        )

        // Use validated FASTQ for downstream
        ch_input = QC_FASTQ.out.validated_fastq
        ch_qc_reports = ch_qc_reports.mix(QC_FASTQ.out.qc_report)

        // Run FastQC on FASTQ files
        FASTQC_FASTQ(ch_fastq)
        ch_fastqc = ch_fastqc.mix(FASTQC_FASTQ.out.zip.map { sample_id, zip -> zip })

        // Run SpecHLA if requested
        if ('spechla' in tools_list) {
            SPECHLA_FASTQ(ch_input)
            ch_results = ch_results.mix(SPECHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'spechla', result_file]
            })
        }

        // Run HLA-HD if requested
        if ('hlahd' in tools_list) {
            HLAHD_FASTQ(ch_input, params.hla_genes)
            ch_results = ch_results.mix(HLAHD_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'hlahd', result_file]
            })
        }

        // Run arcasHLA if requested
        if ('arcashla' in tools_list) {
            ARCASHLA_FASTQ(ch_input)
            ch_results = ch_results.mix(ARCASHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'arcashla', result_file]
            })
        }

        // Run OptiType if requested
        if ('optitype' in tools_list) {
            // Default seq_type parameter
            def seq_type = params.seq_type ?: 'dna'
            OPTITYPE_FASTQ(ch_input, seq_type)
            ch_results = ch_results.mix(OPTITYPE_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'optitype', result_file]
            })
        }
    }

    // Group results by sample
    ch_grouped = ch_results
        .groupTuple(by: 0)
        .map { sample_id, tools, files ->
            [sample_id, tools, files]
        }

    // Run consensus voting
    CONSENSUS(
        ch_grouped,
        params.resolution,
        params.min_tools
    )

    // ===== VISUALIZATION =====
    // Prepare inputs for visualization (join consensus with QC reports)
    ch_vis_input = CONSENSUS.out.consensus
        .join(CONSENSUS.out.comparison)
        .map { sample_id, consensus, comparison ->
            [sample_id, consensus, comparison]
        }

    ch_vis_qc = ch_qc_reports
        .map { sample_id, qc_report -> [sample_id, qc_report] }

    // Run per-sample visualization
    HLA_VISUALIZE(
        ch_vis_input,
        ch_vis_qc
    )

    // ===== SUMMARY REPORT =====
    // Collect all results for summary report
    ch_all_consensus = CONSENSUS.out.consensus.map { sample_id, file -> file }.collect()
    ch_all_comparison = CONSENSUS.out.comparison.map { sample_id, file -> file }.collect()
    ch_all_statistics = HLA_VISUALIZE.out.statistics.map { sample_id, file -> file }.collect()

    // Generate multi-sample summary report
    HLA_SUMMARY_REPORT(
        ch_all_consensus,
        ch_all_comparison,
        ch_all_statistics
    )

    // ===== MULTIQC =====
    // Collect FastQC and HLA reports for MultiQC
    ch_fastqc_collected = ch_fastqc.collect().ifEmpty([])
    ch_hla_reports = HLA_SUMMARY_REPORT.out.statistics

    // Run MultiQC
    MULTIQC(
        ch_fastqc_collected,
        ch_hla_reports,
        file("${projectDir}/assets/multiqc_config.yaml")
    )

    // ===== LOH ANALYSIS (optional) =====
    // Requires SpecHLA and tumor purity/ploidy estimates
    if (params.run_loh && 'spechla' in tools_list && params.tumor_purity && params.tumor_ploidy) {
        log.info "LOH analysis enabled with purity=${params.tumor_purity}, ploidy=${params.tumor_ploidy}"

        // Get SpecHLA output directories for LOH analysis
        if (input_type == 'bam') {
            ch_spechla_for_loh = SPECHLA.out.full_results
        } else {
            ch_spechla_for_loh = SPECHLA_FASTQ.out.full_results
        }

        // Run LOH detection
        HLA_LOH(
            ch_spechla_for_loh,
            params.tumor_purity,
            params.tumor_ploidy
        )

        // Visualize LOH results
        HLA_LOH_VISUALIZE(HLA_LOH.out.loh_results)

        // Generate LOH summary if multiple samples
        ch_all_loh = HLA_LOH.out.loh_results.map { sample_id, file -> file }.collect()
        HLA_LOH_SUMMARY(ch_all_loh)

    } else if (params.run_loh) {
        if (!('spechla' in tools_list)) {
            log.warn "LOH analysis requires SpecHLA. Add 'spechla' to --tools"
        }
        if (!params.tumor_purity || !params.tumor_ploidy) {
            log.warn "LOH analysis requires --tumor_purity and --tumor_ploidy parameters"
        }
    }
}

// Workflow completion handler
workflow.onComplete {
    log.info """
    ===========================================
    Pipeline completed!
    ===========================================
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Work dir    : ${workflow.workDir}
    Output dir  : ${params.outdir}
    ===========================================

    Output locations:
      - QC reports:          ${params.outdir}/<sample>/qc/
      - HLA consensus:       ${params.outdir}/<sample>/
      - Visualizations:      ${params.outdir}/<sample>/visualizations/
      - LOH analysis:        ${params.outdir}/<sample>/loh/ (if enabled)
      - Summary report:      ${params.outdir}/summary/
      - MultiQC report:      ${params.outdir}/multiqc/
    ===========================================
    """
}
