process MERGE_COVERAGE {
    tag "merge_coverage"
    label 'process_low'

    container params.containers.clincnv

    input:
    path cov_files   // collected list from all samples

    output:
    path "merged.cov", emit: merged_cov
    path "versions.yml", emit: versions

    script:
    """
    Rscript \${EBROOTCLINCNV}/mergeFilesFromFolder.R \\
        -i . \\
        -o merged.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clincnv: \$(Rscript -e "packageVersion('ClinCNV')" 2>/dev/null || echo "1.18.3")
    END_VERSIONS
    """
}
