process CLINCNV {
    tag "clincnv"
    label 'process_high'

    container params.containers.clincnv

    publishDir "${params.outdir}/1_clincnv_output", mode: params.publish_dir_mode,
        saveAs: { fn -> fn }

    input:
    path merged_cov
    path annotated_bed
    path ped_file

    output:
    path "cnv_output/",  emit: cnv_output_dir
    path "versions.yml", emit: versions

    script:
    """
    export TMPDIR=\$(pwd)/tmp
    export HOME=\$(pwd)          # host home is not bound in the container
    mkdir -p \${TMPDIR}
    mkdir -p \$(pwd)/cnv_output

    # - PED / trios detection
    # ClinCNV --triosFile: one trio per line, comma-separated, order CHILD,MOTHER,FATHER
    # PED columns: 1=family 2=sample 3=paternal 4=maternal -> print \$2,\$4,\$3
    TRIOS_ARG=""
    if [ -s ${ped_file} ]; then
        grep -v '^#' ${ped_file} \\
            | awk -F'\\t' '\$3!="0" && \$4!="0" {print \$2","\$4","\$3}' \\
            | tr '-' '.' > trios.txt      # ClinCNV converts '-' to '.' in sample names
        if [ -s trios.txt ]; then
            TRIOS_ARG="--triosFile \$(pwd)/trios.txt"
            echo "INFO: Trio-aware mode -- \$(wc -l < trios.txt) complete trios (CHILD,MOTHER,FATHER)"
        else
            echo "INFO: No complete trios found -- running GERMLINE_SINGLE"
        fi
    else
        echo "INFO: No PED file -- running GERMLINE_SINGLE"
    fi

    PAR_ARG=""
    if [ -n "${params.clincnv_par}" ]; then
        PAR_ARG="--par ${params.clincnv_par}"
    fi

    # - Locate ClinCNV inside the container (module var if set, else params.clincnv_dir)
    CLINCNV_DIR=\${EBROOTCLINCNV:-${params.clincnv_dir}}
    if [ ! -f "\${CLINCNV_DIR}/clinCNV.R" ]; then
        echo "ERROR: clinCNV.R not found in \${CLINCNV_DIR}" >&2
        exit 1
    fi
    mkdir -p \$(pwd)/clincnv_run
    cp -r "\${CLINCNV_DIR}/." \$(pwd)/clincnv_run/

    Rscript \$(pwd)/clincnv_run/clinCNV.R \\
        --normal  \$(pwd)/${merged_cov} \\
        --bed     \$(pwd)/${annotated_bed} \\
        --out     \$(pwd)/cnv_output \\
        --scoreG  ${params.clincnv_score_threshold} \\
        --lengthG ${params.clincnv_min_length} \\
        --maxNumGermCNVs ${params.clincnv_max_germ_cnvs} \\
        --maxNumIter 5 \\
        --numberOfThreads ${task.cpus} \\
        --minimumNumOfElemsInCluster ${params.clincnv_min_cluster} \\
        --hg38 \\
        \${PAR_ARG} \\
        --noPlot \\
        \${TRIOS_ARG}

    echo "=== ClinCNV output counts ==="
    echo "Singleton TSVs : \$(find \$(pwd)/cnv_output/normal -name '*_cnvs.tsv' 2>/dev/null | wc -l)"
    echo "Trio CNV files : \$(find \$(pwd)/cnv_output/trios  -name '*_CNVs.txt' 2>/dev/null | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clincnv: \$(git -C "\${CLINCNV_DIR}" log --oneline -1 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
