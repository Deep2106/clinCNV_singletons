process RSCRIPT_FILTER {
    tag "${meta.id}"
    label 'process_single'

    container params.containers.clincnv

    publishDir "${params.outdir}/3_gene_panel_filtered", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vep_vcf)
    path  gene_panel
    path  rscript      // staged from projectDir/bin/

    output:
    tuple val(meta), path("${meta.id}.vep.filtered.tsv"), emit: tsv
    path "versions.yml",                              emit: versions

    script:
    """
    Rscript ${rscript} \\
        ${vep_vcf} \\
        . \\
        ${gene_panel}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | head -1)
        readxl: \$(Rscript -e "packageVersion('readxl')" 2>/dev/null | tr -d "[]'")
    END_VERSIONS
    """
}
