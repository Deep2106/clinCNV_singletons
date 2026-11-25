process GC_ANNOTATE_BED {
    tag "gc_annotate"
    label 'process_low'

    container params.containers.ngsbits

    input:
    path bed
    path fasta
    path fasta_fai

    output:
    path "gcAnnotated.targets.bed", emit: annotated_bed
    path "versions.yml",            emit: versions

    script:
    """
    BedAnnotateGC \\
        -in  ${bed} \\
        -out gc_unsorted.bed \\
        -ref ${fasta}

    # BedCoverage 2025_12 requires sorted BED input
    BedSort -in gc_unsorted.bed -out gc_sorted.bed

    # Detect original BED column count
    # 3-col BED: BedAnnotateGC appends GC as col4 -> correct for ClinCNV (chr,start,end,GC)
    # 4-col BED: BedAnnotateGC appends GC as col5 -> must swap to (chr,start,end,GC,gene)
    #            ClinCNV expects GC at col4, gene at col5
    BED_COLS=\$(awk 'NR==1 && !/^#/ {print NF; exit}' ${bed})

    if [ "\$BED_COLS" -ge 4 ]; then
        # 4-col BED: reorder col5 (GC) before col4 (gene)
        # Input:  chr, start, end, gene, GC
        # Output: chr, start, end, GC,   gene
        awk 'BEGIN{OFS="\\t"} /^#/{print; next} {print \$1,\$2,\$3,\$5,\$4}' gc_sorted.bed > gcAnnotated.targets.bed
    else
        # 3-col BED: GC already at col4, nothing to reorder
        mv gc_sorted.bed gcAnnotated.targets.bed
    fi

    rm -f gc_unsorted.bed gc_sorted.bed 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(BedAnnotateGC --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
