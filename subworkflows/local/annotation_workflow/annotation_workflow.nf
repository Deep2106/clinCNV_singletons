include { VEP_ANNOTATE   } from '../../../modules/local/vep_annotate/singleton_vep_annotate'
include { RSCRIPT_FILTER } from '../../../modules/local/rscript_filter/rscript_filter'

workflow ANNOTATION_WORKFLOW {

    take:
    ch_vcf            // channel: [ meta, vcf ]
    ch_vep_cache      // channel: value(path)
    ch_vep_sv_plugin  // channel: value(path to PATCHED StructuralVariantOverlap.pm)
    ch_sv_overlap     // channel: value(vcf.gz)
    ch_sv_overlap_tbi // channel: value(vcf.gz.tbi)
    ch_gene_panel     // channel: value(path) or NO_FILE
    ch_rscript        // channel: value(path to vep_filter_1000gene.R)

    main:
    ch_versions    = Channel.empty()
    ch_vep_vcf     = Channel.empty()
    ch_sv_qc       = Channel.empty()
    ch_filtered    = Channel.empty()

    if (params.run_vep) {
        VEP_ANNOTATE(
            ch_vcf,
            ch_vep_cache,
            ch_vep_sv_plugin,
            ch_sv_overlap,
            ch_sv_overlap_tbi
        )
        ch_vep_vcf  = VEP_ANNOTATE.out.vcf
        ch_sv_qc    = VEP_ANNOTATE.out.sv_qc
        ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions)

        // - SV-overlap QC: warn on suspiciously low gnomAD-SV hit rate -----
        // qc tsv columns: sample_id, total_cnvs, cnvs_with_overlap, fraction
        ch_sv_qc
            .map { meta, f -> [ meta, f.text.trim().split('\t') ] }
            .subscribe { meta, c ->
                def total = c[1] as int
                def frac  = c[3] as double
                if (total >= 10 && frac < (params.sv_overlap_min_frac as double)) {
                    log.warn "[${meta.id}] only ${c[2]}/${total} CNVs overlap gnomAD-SV " +
                             "(< ${params.sv_overlap_min_frac}) — check VEP SV annotation"
                }
            }

        if (params.run_rfilter) {
            RSCRIPT_FILTER(
                ch_vep_vcf,
                ch_gene_panel,
                ch_rscript
            )
            ch_filtered = RSCRIPT_FILTER.out.tsv
            ch_versions = ch_versions.mix(RSCRIPT_FILTER.out.versions)
        }
    }

    emit:
    vep_vcf  = ch_vep_vcf
    sv_qc    = ch_sv_qc
    filtered = ch_filtered
    versions = ch_versions
}
