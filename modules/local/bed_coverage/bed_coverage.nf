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

    # Detect annotated BED column count to determine what to strip
    # 3-col original BED -> annotated = 4 cols -> body has chr,start,end,GC,cov
    #   -> strip col 4 (GC)              -> output: chr,start,end,cov
    # 4-col original BED -> annotated = 5 cols -> body has chr,start,end,gene,GC,cov
    #   -> strip col 4 (gene) + col 5 (GC) -> output: chr,start,end,cov
    BED_COLS=\$(awk 'NR==1 && !/^#/ {print NF; exit}' ${annotated_bed})

    head -n 1 ontarget_${meta.id}.cov > header.tmp

    if [ "\$BED_COLS" -ge 5 ]; then
        # 4-col original BED: strip gene (col4) and GC (col5)
        tail -n +2 ontarget_${meta.id}.cov | cut -f4,5 --complement > body.tmp
    else
        # 3-col original BED: strip GC only (col4)
        tail -n +2 ontarget_${meta.id}.cov | cut -f4 --complement > body.tmp
    fi

    cat header.tmp body.tmp > ${meta.id}.cov

    rm -f ontarget_${meta.id}.cov header.tmp body.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(BedCoverage --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
