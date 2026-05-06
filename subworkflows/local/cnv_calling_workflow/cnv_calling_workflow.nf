include { CLINCNV  } from '../../../modules/local/clincnv/clincnv'
include { MAKE_VCF } from '../../../modules/local/make_vcf/make_vcf'
include { SORT_VCF } from '../../../modules/local/sort_vcf/sort_vcf'

workflow CNV_CALLING_WORKFLOW {

    take:
    ch_merged_cov
    ch_annotated_bed
    ch_sample_meta
    ch_vcf_header
    ch_ped_file
    ch_fasta_fai

    main:
    ch_versions = Channel.empty()

    CLINCNV(
        ch_merged_cov,
        ch_annotated_bed,
        ch_ped_file
    )
    ch_versions = ch_versions.mix(CLINCNV.out.versions)

    // Use Channel.fromPath on emitted directory  evaluated lazily after CLINCNV
    ch_tsv = CLINCNV.out.cnv_output_dir
        .flatMap { dir ->
            def tsv_files = []
            new File("${dir}/normal").eachDir { sample_dir ->
                sample_dir.eachFileMatch(~/.*_cnvs\.tsv$/) { f ->
                    def sid = f.name.replaceAll('_cnvs\\.tsv$', '')
                    tsv_files << [sid, file(f.absolutePath)]
                }
            }
            return tsv_files
        }

    ch_cnv_matched = ch_tsv
        .join(
            ch_sample_meta.map { meta -> [meta.id, meta] },
            by: 0
        )
        .map { sid, tsv, meta -> [meta, tsv] }

    MAKE_VCF(
        ch_cnv_matched,
        ch_vcf_header,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(MAKE_VCF.out.versions)

    SORT_VCF(
        MAKE_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(SORT_VCF.out.versions)

    emit:
    vcf      = SORT_VCF.out.vcf
    versions = ch_versions
}
