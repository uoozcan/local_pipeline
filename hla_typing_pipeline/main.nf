#!/usr/bin/env nextflow

/*
 * HLA Typing Pipeline
 * Multi-tool HLA typing from BAM/FASTQ with weighted consensus voting
 * Supports both BAM and paired FASTQ inputs
 * Includes FastQC, visualizations, and MultiQC integration
 *
 * Authors: HLA Analysis Team
 * Version: 1.3.0
 */

nextflow.enable.dsl = 2

// Import modules - BAM versions
include { SPECHLA } from './modules/spechla'
include { HLAHD } from './modules/hlahd'
include { HLALA } from './modules/hlala'
include { ARCASHLA } from './modules/arcashla'
include { OPTITYPE } from './modules/optitype'
include { XHLA } from './modules/xhla'
include { HLASCAN } from './modules/hlascan'
include { BAMQC } from './modules/bamqc'

// Import modules - FASTQ versions
include { SPECHLA_FASTQ } from './modules/spechla'
include { HLAHD_FASTQ } from './modules/hlahd'
include { ARCASHLA_FASTQ } from './modules/arcashla'
include { OPTITYPE_FASTQ } from './modules/optitype'
include { XHLA_FASTQ } from './modules/xhla'
include { HLASCAN_FASTQ } from './modules/hlascan'
include { T1K_FASTQ; T1K_LONGREADS } from './modules/t1k'
include { HIFIHLA } from './modules/hifihla'
include { BAMQC_FASTQ } from './modules/bamqc'

// Import QC and consensus modules
include { QC_BAM; QC_FASTQ } from './modules/qc'
include { CONSENSUS } from './modules/consensus'

// Import FastQC modules
include { FASTQC_BAM; FASTQC_FASTQ } from './modules/fastqc'

// Import visualization modules
include { HLA_VISUALIZE; HLA_SUMMARY_REPORT; HLA_PIPELINE_METRICS } from './modules/visualize'

// Import MultiQC module
include { MULTIQC } from './modules/multiqc'

// Import LOH module
include { HLA_LOH; HLA_LOH_VISUALIZE; HLA_LOH_SUMMARY } from './modules/loh'

// Import multi-source module
include { MULTISOURCE_CONSENSUS } from './modules/multisource'

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
        For multi-source (auto-detected when header contains patient_id and seq_type):
                   patient_id,sample_id,seq_type,bam_path,fastq_1,fastq_2
                   seq_type values: WGS, WES, RNAseq, targeted

    Optional arguments:
        --outdir            Output directory (default: ./results)
        --reference         Reference genome: hg38 or hg19 (default: hg38)
        --tools             HLA typing tools to use (default: spechla,hlahd)
                            Options: spechla,hlahd,hlala,arcashla,optitype,xhla,hlascan
                            Note: hlala and xhla work best with BAM input
                            Note: hlascan requires full genome BAM + hg19 reference (v2.1 limitation)
                            In multi-source mode, tools are also filtered by seq_type compatibility
        --seq_type          Sequence type: dna, rna, longreads_hifi, longreads_ont (default: dna)
                            dna / rna: short-read WGS/WES or RNA-seq (paired FASTQ or BAM)
                            longreads_hifi: PacBio HiFi/CCS BAM → T1K + HiFi-HLA (4-field)
                            longreads_ont:  Oxford Nanopore BAM  → T1K (long-read preset)
                            In multi-source mode, RNAseq samples auto-use --rna flag
        --run_bamqc         Enable BAMQC for comprehensive BAM quality control (default: false)
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

    Multi-source integration (activated automatically for multi-source samplesheets):
        --source_weights    Cross-source weights (default: WGS:1.0,WES:0.8,RNAseq:0.6,targeted:0.5)
        --auto_tools_by_seqtype  Auto-select compatible tools per seq_type (default: true)

    LOH (Loss of Heterozygosity) analysis:
        --run_loh           Enable LOH analysis (requires SpecHLA, default: false)
        --tumor_purity      Tumor purity estimate (0-1), required for LOH
        --tumor_ploidy      Tumor ploidy estimate, required for LOH
        --loh_het_cutoff    Minimum het SNPs for LOH call (default: 5)

    Output reports:
        - Per-sample: QC report, consensus HLA types, visualizations
        - Multi-source: Integrated consensus + source comparison per patient (if multi-source input)
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

        # Multi-source (WGS + WES + RNA-seq for same patient) — auto-detected
        nextflow run main.nf --input_samplesheet assets/samplesheet_multisource_template.csv -profile singularity

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
def is_multisource = false
if (params.input_bam) {
    input_type = 'bam'
} else if (params.input_fastq_1 && params.input_fastq_2) {
    input_type = 'fastq'
} else if (params.input_samplesheet) {
    // Determine from samplesheet header
    def samplesheet = file(params.input_samplesheet)
    def header = samplesheet.readLines()[0]
    if (header.contains('patient_id') && header.contains('seq_type')) {
        // Multi-source samplesheet: groups multiple seq types per patient
        input_type = 'multisource'
        is_multisource = true
    } else if (header.contains('fastq_1') || header.contains('fastq1') || header.contains('fq1')) {
        input_type = 'fastq'
    } else {
        input_type = 'bam'
    }
} else {
    log.error "Please provide input: --input_bam, --input_fastq_1/--input_fastq_2, or --input_samplesheet"
    helpMessage()
    exit 1
}

// Tool compatibility map for multi-source mode
// hlascan: WGS/WES only (requires full genome BAM + hg19); excluded from RNA/targeted
// longreads_hifi: PacBio HiFi BAM  → T1K (long-read mode) + HiFi-HLA (4-field)
// longreads_ont:  Oxford Nanopore  → T1K (long-read mode)
def SEQ_TYPE_TOOLS = [
    'WGS'            : ['spechla', 'hlahd', 'hlala', 'arcashla', 'optitype', 'xhla', 'hlascan', 't1k'],
    'WES'            : ['spechla', 'hlahd', 'optitype', 'xhla', 'hlascan', 't1k'],
    'RNAseq'         : ['arcashla', 'optitype', 't1k'],
    'targeted'       : ['optitype', 'hlahd'],
    'longreads_hifi' : ['t1k', 'hifihla'],
    'longreads_ont'  : ['t1k'],
]

// Resolve effective tools for a given seq_type and global tools_list
def resolveTools = { String seq_type, List global_tools ->
    def compatible = SEQ_TYPE_TOOLS.get(seq_type, SEQ_TYPE_TOOLS['WGS'])
    if (global_tools && !global_tools.isEmpty()) {
        def intersected = global_tools.findAll { it in compatible }
        if (intersected.isEmpty()) {
            log.warn "No compatible tools for seq_type=${seq_type} among --tools=${global_tools.join(',')}. Falling back to: ${compatible.join(',')}"
            return compatible
        }
        return intersected
    }
    return compatible
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
Input type      : ${input_type.toUpperCase()}${is_multisource ? ' (multi-source)' : ''}
Input           : ${params.input_bam ?: params.input_fastq_1 ?: params.input_samplesheet}
Output          : ${params.outdir}
Reference       : ${params.reference}
Tools           : ${tools_list.join(', ')}${is_multisource ? ' (filtered per seq_type)' : ''}
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

// Create input channel for multi-source samplesheet
// Returns: [patient_id, sample_id, seq_type, bam_or_null, fq1_or_null, fq2_or_null]
def create_multisource_channel() {
    return Channel
        .fromPath(params.input_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def patient_id = row.patient_id
            def sample_id  = row.sample_id
            def seq_type   = row.seq_type
            def bam_path   = (row.bam_path  && row.bam_path.trim())  ? file(row.bam_path)  : null
            def fq1_path   = (row.fastq_1   && row.fastq_1.trim())   ? file(row.fastq_1)   : null
            def fq2_path   = (row.fastq_2   && row.fastq_2.trim())   ? file(row.fastq_2)   : null
            [patient_id, sample_id, seq_type, bam_path, fq1_path, fq2_path]
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

        // Run BAMQC if requested
        if (params.run_bamqc) {
            BAMQC(ch_bam)
        }

        // Run xHLA if requested
        if ('xhla' in tools_list) {
            XHLA(ch_input)
            ch_results = ch_results.mix(XHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'xhla', result_file]
            })
        }

        // Run HLAscan if requested (BAM input; hg19 only in v2.1 container)
        if ('hlascan' in tools_list) {
            HLASCAN(ch_input)
            ch_results = ch_results.mix(HLASCAN.out.results.map { sample_id, result_file ->
                [sample_id, 'hlascan', result_file]
            })
        }

        // ── Long-read tools (activated via --seq_type longreads_hifi / longreads_ont) ──
        // Short-read tools above are implicitly skipped because resolveTools() returns
        // only ['t1k', 'hifihla'] or ['t1k'] for these seq_types; the `in tools_list`
        // guards above will still match if the user explicitly listed those tools, but
        // for typical use the long-read seq_types will not include short-read tools.

        // T1K long-read mode (HiFi or ONT)
        if ('t1k' in tools_list && params.seq_type in ['longreads_hifi', 'longreads_ont']) {
            def lr_platform = params.seq_type == 'longreads_hifi' ? 'hifi' : 'ont'
            T1K_LONGREADS(ch_input, lr_platform)
            ch_results = ch_results.mix(T1K_LONGREADS.out.results.map { sample_id, result_file ->
                [sample_id, 't1k', result_file]
            })
        }

        // HiFi-HLA (4-field resolution; PacBio HiFi only)
        if ('hifihla' in tools_list && params.seq_type == 'longreads_hifi') {
            HIFIHLA(ch_input)
            ch_results = ch_results.mix(HIFIHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'hifihla', result_file]
            })
        }

    } else if (input_type == 'fastq') {
        // ===== FASTQ INPUT WORKFLOW =====
        ch_fastq = create_fastq_channel()

        // Run QC on FASTQ files
        QC_FASTQ(
            ch_fastq,
            params.min_hla_reads,
            params.min_read_length
        )

        // Use raw FASTQ for HLA typing; QC runs independently for reporting only
        ch_input = ch_fastq
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

        // Run xHLA if requested (note: xHLA works best with BAM, FASTQ support is limited)
        if ('xhla' in tools_list) {
            XHLA_FASTQ(ch_input)
            ch_results = ch_results.mix(XHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'xhla', result_file]
            })
        }

        // Run HLAscan if requested (FASTQ mode; limited by gene-resolution issue with HLA-enriched reads)
        if ('hlascan' in tools_list) {
            HLASCAN_FASTQ(ch_input)
            ch_results = ch_results.mix(HLASCAN_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'hlascan', result_file]
            })
        }

        // Run T1K if requested (WGS/RNA-seq, Class I + II, fast ~5 min)
        if ('t1k' in tools_list) {
            T1K_FASTQ(ch_input)
            ch_results = ch_results.mix(T1K_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 't1k', result_file]
            })
        }

    } else if (input_type == 'multisource') {
        // ===== MULTI-SOURCE WORKFLOW =====
        // Accepts a samplesheet with patient_id + seq_type columns.
        // BAM and FASTQ rows are processed in parallel; each sample runs tools
        // appropriate for its seq_type. All per-sample results feed the shared
        // CONSENSUS process, then MULTISOURCE_CONSENSUS integrates them per patient.

        ch_ms = create_multisource_channel()

        // ---- BAM-based samples ----
        ch_ms_bam = ch_ms
            .filter { patient_id, sample_id, seq_type, bam, fq1, fq2 -> bam != null }
            .map    { patient_id, sample_id, seq_type, bam, fq1, fq2 ->
                [sample_id, bam, seq_type, patient_id]
            }

        if (ch_ms_bam) {
            QC_BAM(
                ch_ms_bam.map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                params.reference, params.min_hla_reads, params.min_read_length
            )
            ch_qc_reports = ch_qc_reports.mix(QC_BAM.out.qc_report)

            // Re-attach seq_type and patient_id to validated BAM
            ch_ms_bam_validated = QC_BAM.out.validated_bam
                .join(ch_ms_bam.map { sample_id, bam, seq_type, patient_id ->
                    [sample_id, seq_type, patient_id]
                })
                // emits: [sample_id, validated_bam, seq_type, patient_id]

            FASTQC_BAM(ch_ms_bam.map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] })
            ch_fastqc = ch_fastqc.mix(FASTQC_BAM.out.zip.map { sample_id, zip -> zip })

            if (params.run_bamqc) {
                BAMQC(ch_ms_bam.map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] })
            }

            // SpecHLA — WGS, WES
            SPECHLA(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'spechla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                params.reference
            )
            ch_results = ch_results.mix(SPECHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'spechla', result_file]
            })

            // HLA-HD — WGS, WES, targeted
            HLAHD(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'hlahd' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                params.reference, params.hla_genes
            )
            ch_results = ch_results.mix(HLAHD.out.results.map { sample_id, result_file ->
                [sample_id, 'hlahd', result_file]
            })

            // HLA*LA — WGS only
            HLALA(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'hlala' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                params.hlala_graph
            )
            ch_results = ch_results.mix(HLALA.out.results.map { sample_id, result_file ->
                [sample_id, 'hlala', result_file]
            })

            // arcasHLA (BAM mode) — WGS
            ARCASHLA(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'arcashla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                params.reference
            )
            ch_results = ch_results.mix(ARCASHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'arcashla', result_file]
            })

            // OptiType (BAM mode) — WGS, WES, targeted
            OPTITYPE(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'optitype' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] }
            )
            ch_results = ch_results.mix(OPTITYPE.out.results.map { sample_id, result_file ->
                [sample_id, 'optitype', result_file]
            })

            // xHLA — WGS, WES
            XHLA(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'xhla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] }
            )
            ch_results = ch_results.mix(XHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'xhla', result_file]
            })

            // T1K long-read mode — longreads_hifi / longreads_ont BAM samples
            ch_lr_bam = ch_ms_bam_validated
                .filter { sample_id, bam, seq_type, patient_id ->
                    't1k' in resolveTools(seq_type, tools_list) &&
                    seq_type in ['longreads_hifi', 'longreads_ont']
                }
            T1K_LONGREADS(
                ch_lr_bam.map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] },
                ch_lr_bam.map { sample_id, bam, seq_type, patient_id ->
                    seq_type == 'longreads_hifi' ? 'hifi' : 'ont'
                }
            )
            ch_results = ch_results.mix(T1K_LONGREADS.out.results.map { sample_id, result_file ->
                [sample_id, 't1k', result_file]
            })

            // HiFi-HLA — longreads_hifi BAM samples only
            HIFIHLA(
                ch_ms_bam_validated
                    .filter { sample_id, bam, seq_type, patient_id ->
                        'hifihla' in resolveTools(seq_type, tools_list) &&
                        seq_type == 'longreads_hifi'
                    }
                    .map { sample_id, bam, seq_type, patient_id -> [sample_id, bam] }
            )
            ch_results = ch_results.mix(HIFIHLA.out.results.map { sample_id, result_file ->
                [sample_id, 'hifihla', result_file]
            })
        }

        // ---- FASTQ-based samples (RNAseq, targeted with FASTQ) ----
        ch_ms_fastq = ch_ms
            .filter { patient_id, sample_id, seq_type, bam, fq1, fq2 ->
                bam == null && fq1 != null
            }
            .map { patient_id, sample_id, seq_type, bam, fq1, fq2 ->
                [sample_id, fq1, fq2, seq_type, patient_id]
            }

        if (ch_ms_fastq) {
            QC_FASTQ(
                ch_ms_fastq.map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] },
                params.min_hla_reads, params.min_read_length
            )
            ch_qc_reports = ch_qc_reports.mix(QC_FASTQ.out.qc_report)

            ch_ms_fastq_validated = QC_FASTQ.out.validated_fastq
                .join(ch_ms_fastq.map { sample_id, fq1, fq2, seq_type, patient_id ->
                    [sample_id, seq_type, patient_id]
                })
                // emits: [sample_id, fq1, fq2, seq_type, patient_id]

            FASTQC_FASTQ(ch_ms_fastq.map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] })
            ch_fastqc = ch_fastqc.mix(FASTQC_FASTQ.out.zip.map { sample_id, zip -> zip })

            // arcasHLA (FASTQ) — RNAseq
            ARCASHLA_FASTQ(
                ch_ms_fastq_validated
                    .filter { sample_id, fq1, fq2, seq_type, patient_id ->
                        'arcashla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] }
            )
            ch_results = ch_results.mix(ARCASHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'arcashla', result_file]
            })

            // HLA-HD (FASTQ) — targeted
            HLAHD_FASTQ(
                ch_ms_fastq_validated
                    .filter { sample_id, fq1, fq2, seq_type, patient_id ->
                        'hlahd' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] },
                params.hla_genes
            )
            ch_results = ch_results.mix(HLAHD_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'hlahd', result_file]
            })

            // SpecHLA (FASTQ) — targeted (if selected)
            SPECHLA_FASTQ(
                ch_ms_fastq_validated
                    .filter { sample_id, fq1, fq2, seq_type, patient_id ->
                        'spechla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] }
            )
            ch_results = ch_results.mix(SPECHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'spechla', result_file]
            })

            // OptiType (FASTQ) — RNAseq (auto --rna) and targeted (--dna)
            // Use multiMap to derive per-sample mode without ordering assumptions
            ch_ms_fastq_opti = ch_ms_fastq_validated
                .filter { sample_id, fq1, fq2, seq_type, patient_id ->
                    'optitype' in resolveTools(seq_type, tools_list)
                }
            ch_ms_fastq_opti
                .multiMap { sample_id, fq1, fq2, seq_type, patient_id ->
                    reads: [sample_id, fq1, fq2]
                    mode:  seq_type == 'RNAseq' ? 'rna' : 'dna'
                }
                .set { ch_opti_ms }

            OPTITYPE_FASTQ(ch_opti_ms.reads, ch_opti_ms.mode)
            ch_results = ch_results.mix(OPTITYPE_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'optitype', result_file]
            })

            // xHLA (FASTQ) — targeted (if selected)
            XHLA_FASTQ(
                ch_ms_fastq_validated
                    .filter { sample_id, fq1, fq2, seq_type, patient_id ->
                        'xhla' in resolveTools(seq_type, tools_list)
                    }
                    .map { sample_id, fq1, fq2, seq_type, patient_id -> [sample_id, fq1, fq2] }
            )
            ch_results = ch_results.mix(XHLA_FASTQ.out.results.map { sample_id, result_file ->
                [sample_id, 'xhla', result_file]
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

    // ===== MULTI-SOURCE INTEGRATION (per-patient cross-source consensus) =====
    if (is_multisource) {
        log.info "Running multi-source integration (MULTISOURCE_CONSENSUS) per patient..."

        // Re-read samplesheet to associate sample_id -> patient_id and seq_type
        ch_sample_metadata = Channel
            .fromPath(params.input_samplesheet)
            .splitCsv(header: true)
            .map { row -> [row.sample_id, row.patient_id, row.seq_type] }

        // Join per-sample consensus output with samplesheet metadata
        // CONSENSUS emits: [sample_id, consensus_file]
        // Result: [sample_id, consensus_file, patient_id, seq_type]
        ch_per_source = CONSENSUS.out.consensus
            .join(ch_sample_metadata, by: 0)
            .map { sample_id, consensus_file, patient_id, seq_type ->
                [patient_id, sample_id, seq_type, consensus_file]
            }

        // Group all sources per patient: [patient_id, [sample_ids], [seq_types], [files]]
        ch_patient_grouped = ch_per_source
            .groupTuple(by: 0)

        // Run cross-source integration per patient
        MULTISOURCE_CONSENSUS(
            ch_patient_grouped,
            params.source_weights,
            params.resolution
        )
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
      - Multi-source:        ${params.outdir}/<patient_id>/integrated/ (if multi-source input)
      - Summary report:      ${params.outdir}/summary/
      - MultiQC report:      ${params.outdir}/multiqc/
      - Pipeline metrics:    ${params.outdir}/pipeline_info/execution_metrics_report.html
    ===========================================
    """

    // Generate pipeline execution metrics after trace file is finalised
    if (workflow.success) {
        def pipelineInfoDir = file("${params.outdir}/pipeline_info")
        def traceFiles = pipelineInfoDir.listFiles()?.findAll { it.name.startsWith("trace_") && it.name.endsWith(".txt") }
        if (traceFiles) {
            def latestTrace = traceFiles.sort { it.lastModified() }.last()
            log.info "Generating pipeline execution metrics from ${latestTrace.name} ..."
            def cmd = "hla_pipeline_metrics.py --trace ${latestTrace} --outdir ${pipelineInfoDir}"
            def proc = ["bash", "-c", cmd].execute()
            proc.waitForOrKill(120_000)  // wait up to 2 min
            if (proc.exitValue() == 0) {
                log.info "Pipeline metrics saved to ${pipelineInfoDir}/execution_metrics_report.html"
            } else {
                log.warn "Pipeline metrics script exited with code ${proc.exitValue()} — check container has plotly installed"
            }
        } else {
            log.warn "No trace file found in ${pipelineInfoDir} — skipping pipeline metrics"
        }
    }
}
