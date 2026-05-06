include { GC_ANNOTATE_BED  } from '../../../modules/local/gc_annotate/gc_annotate_bed'
include { INFER_SEX        } from '../../../modules/local/infer_sex/infer_sex'
include { BED_COVERAGE     } from '../../../modules/local/bed_coverage/bed_coverage'
include { MERGE_COVERAGE   } from '../../../modules/local/merge_coverage/merge_coverage'
include { ARCHIVE_COVERAGE } from '../../../modules/local/archive_coverage/archive_coverage'

workflow COVERAGE_WORKFLOW {

    take:
    ch_bam         // channel: [ meta, bam, bai ]
    ch_bed         // channel: value(cnv_bed)  gene names in col 4 (comma-separated)
    ch_fasta       // channel: value(fasta)
    ch_fasta_fai   // channel: value(fasta.fai)

    main:
    ch_versions = Channel.empty()

    // Step 1: Annotate BED with GC content  runs once
    GC_ANNOTATE_BED(
        ch_bed,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(GC_ANNOTATE_BED.out.versions)

    // Step 2: Infer sex only for samples where sex is unknown (0)
    //   - Known sex (1/2): keep samplesheet value, skip SampleGender
    //   - Unknown sex (0): run SampleGender and replace with inferred value
    ch_known_sex   = ch_bam.filter { meta, bam, bai -> meta.sex != null && meta.sex != 0 }
    ch_unknown_sex = ch_bam.filter { meta, bam, bai -> meta.sex == null || meta.sex == 0  }

    INFER_SEX(ch_unknown_sex)
    ch_versions = ch_versions.mix(INFER_SEX.out.versions)

    // Parse SampleGender output and update meta.sex for unknown-sex samples
    ch_inferred = ch_unknown_sex
        .join(INFER_SEX.out.sex, by: 0)
        .map { meta, bam, bai, sex_txt ->
            def inferred_sex = 0
            sex_txt.eachLine { line ->
                if (!line.startsWith('#')) {
                    def cols = line.split('\t')
                    if (cols.size() >= 3) {
                        inferred_sex = cols[2].trim() == 'male'   ? 1 :
                                       cols[2].trim() == 'female' ? 2 : 0
                    }
                }
            }
            [ meta + [sex: inferred_sex], bam, bai ]
        }

    // Merge known-sex and inferred-sex channels
    ch_sex_updated = ch_known_sex.mix(ch_inferred)

    // Step 3: Per-sample coverage (scatter)  current batch only
    BED_COVERAGE(
        ch_sex_updated,
        GC_ANNOTATE_BED.out.annotated_bed
    )
    ch_versions = ch_versions.mix(BED_COVERAGE.out.versions)

    // Step 4: Build full cohort coverage channel
    ch_current_cov  = BED_COVERAGE.out.cov
    ch_previous_cov = params.previous_cov_dir
        ? Channel.fromPath("${params.previous_cov_dir}/*.cov", checkIfExists: false)
        : Channel.empty()
    ch_all_cov = ch_current_cov.mix(ch_previous_cov)

    // Step 5: Archive current batch .cov files to persistent store
    ARCHIVE_COVERAGE(
        ch_current_cov.collect()
    )

    // Step 6: Merge full cohort
    MERGE_COVERAGE(
        ch_all_cov.collect()
    )
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)

    emit:
    annotated_bed  = GC_ANNOTATE_BED.out.annotated_bed
    merged_cov     = MERGE_COVERAGE.out.merged_cov
    bam_with_sex   = ch_sex_updated   // updated meta with resolved sex
    versions       = ch_versions
}
