process SORT_VCF {
    tag "${meta.id}"
    label 'process_single'

    container params.containers.bcftools

    publishDir "${params.outdir}/1_clincnv_output/vcfs", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.sorted.vcf"), emit: vcf
    path "versions.yml",                            emit: versions

    script:
    """
    bcftools sort -T . ${vcf} -o ${meta.id}.sorted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
