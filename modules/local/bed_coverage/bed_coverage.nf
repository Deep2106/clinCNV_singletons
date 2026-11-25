process BED_COVERAGE {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.ngsbits

    input:
    tuple val(meta), path(bam), path(bai)
    path  annotated_bed

    output:
    path "${meta.id}.cov", emit: cov
    path "versions.yml",   emit: versions

    script:
    """
    BedCoverage \\
        -bam      ${bam} \\
        -in       ${annotated_bed} \\
        -min_mapq 5 \\
        -decimals 4 > ontarget_${meta.id}.cov

    # ClinCNV .cov = chr,start,end,coverage (GC/gene come from --bed).
    # Keep cols 1-3 + last column: width-independent of BED layout.
    awk -F'\\t' 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$NF}' ontarget_${meta.id}.cov > ${meta.id}.cov

    # Guard: 4 cols, numeric coverage, non-empty
    if ! awk -F'\\t' 'NF!=4 || (NR>1 && \$4!=\$4+0) {exit 1} END {if (NR<2) exit 1}' ${meta.id}.cov; then
        echo "ERROR: malformed ${meta.id}.cov" >&2
        head -3 ${meta.id}.cov >&2
        exit 1
    fi

    rm -f ontarget_${meta.id}.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(BedCoverage --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
