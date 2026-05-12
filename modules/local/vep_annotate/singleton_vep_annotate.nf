process VEP_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    container params.containers.vep

    publishDir "${params.outdir}/2_vep_annotated", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path  vep_cache
    path  vep_plugins
    path  sv_overlap_vcf
    path  sv_overlap_tbi

    output:
    tuple val(meta), path("${meta.id}.vep.vcf"), emit: vcf
    path "versions.yml",                         emit: versions

    script:
    """
    vep \\
        --offline \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --assembly GRCh38 \\
        --species homo_sapiens \\
        --format vcf \\
        --vcf \\
        --force_overwrite \\
        --no_stats \\
        --mane_select \\
        --plugin StructuralVariantOverlap,file=${sv_overlap_vcf} \\
        --dir_plugins ${vep_plugins} \\
        -i ${vcf} \\
        -o ${meta.id}.vep.vcf \\
        --fork ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep "ensembl-vep" | awk '{print \$NF}')
    END_VERSIONS
    """
}
