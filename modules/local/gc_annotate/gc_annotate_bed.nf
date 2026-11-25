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

    # Locate GC from OUTPUT content, not input column count.
    # BedAnnotateGC column placement varies by ngs-bits version and input layout
    # (3-col, 4-col vendor, or pre-annotated 5-col chr/start/end/gc/gene).
    # is_gc_col N: true if column N is numeric in [0,1] on every non-header line
    is_gc_col() {
        awk -F'\\t' -v c="\$1" '
            !/^#/ { if (NF < c || !(\$c == \$c+0 && \$c >= 0 && \$c <= 1)) { bad=1; exit } }
            END   { exit bad }' gc_sorted.bed
    }

    if is_gc_col 4; then
        echo "INFO: GC at col4 -- no reorder"
        mv gc_sorted.bed gcAnnotated.targets.bed
    elif is_gc_col 5; then
        echo "INFO: GC at col5 -- reordering to chr,start,end,GC,gene"
        awk 'BEGIN{OFS="\\t"} /^#/{print; next} {print \$1,\$2,\$3,\$5,\$4}' gc_sorted.bed > gcAnnotated.targets.bed
    else
        echo "ERROR: no valid GC column (col4/col5) in BedAnnotateGC output" >&2
        head -3 gc_sorted.bed >&2
        exit 1
    fi

    # Final guard: ClinCNV requires GC at col4
    mv gcAnnotated.targets.bed gc_check.bed
    cp gc_check.bed gc_sorted.bed
    is_gc_col 4 || { echo "ERROR: col4 is not GC after processing" >&2; head -3 gc_check.bed >&2; exit 1; }
    mv gc_check.bed gcAnnotated.targets.bed

    rm -f gc_unsorted.bed gc_sorted.bed 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(BedAnnotateGC --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
