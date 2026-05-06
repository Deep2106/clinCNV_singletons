#!/usr/bin/env nextflow

/*
========================================================================================
    CNV CALLING PIPELINE (ClinCNV WES) v1.0
========================================================================================
    Production pipeline:
    - GC-content BED annotation (once)
    - Per-sample coverage calculation (scatter)
    - Coverage merge (gather)
    - ClinCNV cohort-level CNV calling
    - TSV → VCF conversion (scatter, fixed CN logic)
    - bcftools sort
    - Optional: VEP annotation
    - Optional: R gene panel filter (requires VEP)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

import groovy.transform.Field

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================================================
     CNV Calling Pipeline v${workflow.manifest.version}
    =========================================================================

    Usage:
        nextflow run main.nf -profile slurm --input samplesheet.csv [options]

    Required:
        --input             Sample CSV (sample_id,bam,bai,sex)
        --fasta             Reference genome FASTA
        --cnv_bed           CNV target BED with gene names in col 4 (comma-separated)
                            e.g. S000059_cexome6_hg38_cnvs_targets.bed
                            Do NOT use SNV padded BED
        --ped_file          PED file for trio/duo aware calling (optional)
                            null = singleton/cohort mode only
                            Always add complete family units per batch

    Optional:
        --outdir            Output directory (default: ./results)

    VEP Annotation (--run_vep true):
        --vep_cache         VEP cache directory
        --vep_plugins       VEP plugins directory
        --vep_cache_version VEP cache version (default: 107)
        --sv_overlap_vcf    1000G SV VCF for StructuralVariantOverlap plugin

    Gene Panel Filter (--run_rfilter true, requires --run_vep true):
        --gene_panel        Gene panel XLSX with 'Gene Symbol' column

    ClinCNV tuning:
        --clincnv_score_threshold   Minimum score (default: 50)
        --clincnv_min_cluster       Minimum cluster size (default: 5)

    Sample CSV format:
        sample_id,bam,bai,sex
        S000021_S15423Nr1,/path/to/sample.bam,/path/to/sample.bai,1

        sex: 1=male(XY), 2=female(XX), 0/blank=unknown

    Outputs:
        results/1_clincnv_output/       Always produced
        results/2_vep_annotated/        If --run_vep true
        results/3_gene_panel_filtered/  If --run_vep true AND --run_rfilter true
    =========================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    INCLUDE SUBWORKFLOWS
========================================================================================
*/

include { COVERAGE_WORKFLOW      } from './subworkflows/local/coverage_workflow/coverage_workflow'
include { CNV_CALLING_WORKFLOW   } from './subworkflows/local/cnv_calling_workflow/cnv_calling_workflow'
include { ANNOTATION_WORKFLOW    } from './subworkflows/local/annotation_workflow/annotation_workflow'

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validateParameters() {
    def errors = []

    if (!params.input)
        errors << "ERROR: --input samplesheet required"
    else if (!file(params.input).exists())
        errors << "ERROR: Samplesheet not found: ${params.input}"

    if (!params.fasta)
        errors << "ERROR: --fasta required"

    if (!params.cnv_bed)
        errors << "ERROR: --cnv_bed required (CNV target BED with comma-separated gene names in col 4)"
    else if (!file(params.cnv_bed).exists())
        errors << "ERROR: cnv_bed not found: ${params.cnv_bed}"

    if (!params.fasta_fai)
        errors << "ERROR: --fasta_fai required"

    if (params.run_rfilter && !params.run_vep)
        errors << "ERROR: --run_rfilter requires --run_vep true. " +
                  "Gene panel filter depends on VEP-annotated output."

    if (params.run_rfilter && !params.gene_panel)
        errors << "ERROR: --gene_panel (XLSX) required when --run_rfilter true"

    if (params.run_vep) {
        if (!params.vep_cache)       errors << "ERROR: --vep_cache required when --run_vep true"
        if (!params.vep_plugins)     errors << "ERROR: --vep_plugins required when --run_vep true"
        if (!params.sv_overlap_vcf)  errors << "ERROR: --sv_overlap_vcf required when --run_vep true"
    }

    if (errors) {
        log.error "\nPARAMETER VALIDATION FAILED\n" + errors.join("\n")
        exit 1
    }
    log.info "Parameters validated successfully"
}

validateParameters()

/*
========================================================================================
    PRINT BANNER
========================================================================================
*/

log.info """
╔═══════════════════════════════════════════════════════════════════╗
║       ClinCNV WES CNV Calling Pipeline v${workflow.manifest.version}                  ║
╠═══════════════════════════════════════════════════════════════════╣
║  Input samplesheet  : ${params.input}
║  Reference FASTA    : ${params.fasta ?: 'NOT SET'}
║  CNV target BED     : ${params.cnv_bed ?: 'NOT SET'}
║  PED file           : ${params.ped_file ?: 'Not provided (singleton mode)'}
║  ClinCNV score ≥    : ${params.clincnv_score_threshold}
║  ClinCNV min cluster: ${params.clincnv_min_cluster}
║  Run VEP            : ${params.run_vep}
║  Run R filter       : ${params.run_rfilter}
║  Gene panel         : ${params.gene_panel ?: 'N/A'}
║  Output directory   : ${params.outdir}
╚═══════════════════════════════════════════════════════════════════╝
""".stripIndent()

/*
========================================================================================
    HELPER: PARSE SAMPLESHEET
========================================================================================
*/

def parseSamplesheet(csv_file) {
    Channel.fromPath(csv_file)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:  row.sample_id,
                sex: row.sex ? row.sex.toInteger() : 0
            ]
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            [ meta, bam, bai ]
        }
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    ch_versions      = Channel.empty()

    // Ensure assets directory exists
    file("${projectDir}/assets").mkdirs()
    if (!file("${projectDir}/assets/NO_FILE").exists()) {
        file("${projectDir}/assets/NO_FILE").text = ""
    }

    /*
    ============================================================================
        STAGE 1: PARSE SAMPLESHEET
    ============================================================================
    */

    ch_bam = parseSamplesheet(params.input)

    /*
    ============================================================================
        STAGE 2: COVERAGE (GC annotate → per-sample BedCoverage → merge)
    ============================================================================
    */

    ch_fasta      = Channel.value(file(params.fasta))
    ch_fasta_fai  = Channel.value(file(params.fasta_fai, checkIfExists: true))
    ch_vcf_header = Channel.value(file(params.vcf_header))
    ch_cnv_bed    = Channel.value(file(params.cnv_bed, checkIfExists: true))

    COVERAGE_WORKFLOW(
        ch_bam,
        ch_cnv_bed,
        ch_fasta,
        ch_fasta_fai
    )
    ch_annotated_bed = COVERAGE_WORKFLOW.out.annotated_bed
    ch_merged_cov    = COVERAGE_WORKFLOW.out.merged_cov
    ch_versions      = ch_versions.mix(COVERAGE_WORKFLOW.out.versions)

    /*
    ============================================================================
        STAGE 3: CNV CALLING (ClinCNV → TSV→VCF → sort)
    ============================================================================
    */

    // Use sex-resolved meta from COVERAGE_WORKFLOW output
    ch_sample_meta = COVERAGE_WORKFLOW.out.bam_with_sex.map { meta, bam, bai -> meta }
    ch_ped_file    = params.ped_file
        ? Channel.value(file(params.ped_file, checkIfExists: true))
        : Channel.value(file("${projectDir}/assets/NO_FILE"))

    CNV_CALLING_WORKFLOW(
        ch_merged_cov,
        ch_annotated_bed,
        ch_sample_meta,
        ch_vcf_header,
        ch_ped_file,
        ch_fasta_fai
    )

    ch_sorted_vcf = CNV_CALLING_WORKFLOW.out.vcf
    ch_versions   = ch_versions.mix(CNV_CALLING_WORKFLOW.out.versions)
    // ch_sample_meta now carries resolved sex (inferred if was 0 in samplesheet)

    /*
    ============================================================================
        STAGE 4: ANNOTATION (VEP + R filter  both optional)
    ============================================================================
    */

    ch_vep_cache      = params.run_vep ? Channel.value(file(params.vep_cache,      checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_vep_plugins    = params.run_vep ? Channel.value(file(params.vep_plugins,    checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_sv_overlap     = params.run_vep ? Channel.value(file(params.sv_overlap_vcf, checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_sv_overlap_tbi = params.run_vep ? Channel.value(file("${params.sv_overlap_vcf}.tbi", checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_gene_panel     = params.run_rfilter ? Channel.value(file(params.gene_panel, checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_rscript        = Channel.value(file("${projectDir}/bin/vep_filter_1000gene.R", checkIfExists: true))

    ANNOTATION_WORKFLOW(
        ch_sorted_vcf,
        ch_vep_cache,
        ch_vep_plugins,
        ch_sv_overlap,
        ch_sv_overlap_tbi,
        ch_gene_panel,
        ch_rscript
    )

    ch_versions = ch_versions.mix(ANNOTATION_WORKFLOW.out.versions)

    /*
    ============================================================================
        STAGE 5: SOFTWARE VERSIONS
    ============================================================================
    */

    ch_versions
        .flatten()
        .filter { it != null }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

/*
========================================================================================
    COMPLETION HANDLERS
========================================================================================
*/

workflow.onComplete {
    log.info """
╔═══════════════════════════════════════════════════════════════════╗
║                    PIPELINE COMPLETED                             ║
╠═══════════════════════════════════════════════════════════════════╣
║  Status   : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
║  Duration : ${workflow.duration}
║  Output   : ${params.outdir}
╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
╔═══════════════════════════════════════════════════════════════════╗
║                     PIPELINE FAILED                               ║
║  Error: ${workflow.errorMessage}
╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
