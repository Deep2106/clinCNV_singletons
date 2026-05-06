process INFER_SEX {
    tag "${meta.id}"
    label 'process_single'

    container params.containers.ngsbits

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.sex.txt"), emit: sex
    path "versions.yml",                         emit: versions

    script:
    """
    SampleGender \\
        -in     ${bam} \\
        -method xy \\
        -build  hg38 \\
        -out    ${meta.id}.sex.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(SampleGender --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
