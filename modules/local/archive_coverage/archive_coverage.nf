process ARCHIVE_COVERAGE {
    tag "archive_coverage"
    label 'process_single'

    // Shell-only, no container needed  runs on executor node
    // Uses cp not mv so Nextflow work dir is unaffected

    input:
    path cov_files   // collected current-batch .cov files

    output:
    path "archived.flag", emit: flag   // sentinel so downstream can depend on this

    script:
    def archive = params.previous_cov_dir
    """
    mkdir -p ${archive}

    for f in *.cov; do
        dest="${archive}/\${f}"
        if [ -f "\${dest}" ]; then
            echo "WARNING: \${f} already exists in archive  skipping (sample already in cohort?)" >&2
        else
            cp "\${f}" "${archive}/\${f}"
            echo "Archived: \${f}"
        fi
    done

    echo "Archived at: \$(date)" > archived.flag
    """
}
