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
    mkdir -p \${TMPDIR}
    mkdir -p \$(pwd)/cnv_output

    # ── PED detection (bash only  empty file = singleton mode)
    # trios.txt must be comma-separated: ClinCNV reads with sep=","
    TRIOS_ARG=""
    if [ -s ${ped_file} ]; then
        grep -v '^#' ${ped_file} \\
            | awk -F'\\t' '\$3!="0" && \$4!="0" {print \$2","\$3","\$4}' > trios.txt
        if [ -s trios.txt ]; then
            TRIOS_ARG="--triosFile \$(pwd)/trios.txt"
            echo "INFO: Trio-aware mode -- \$(wc -l < trios.txt) complete trios"
        else
            echo "INFO: No complete trios found -- running GERMLINE_SINGLE"
        fi
    else
        echo "INFO: No PED file -- running GERMLINE_SINGLE"
    fi

    cp -r \${EBROOTCLINCNV}/. \$(pwd)/clincnv_run/

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
        --noPlot \\
        \${TRIOS_ARG}

    echo "=== ClinCNV TSV count ==="
    find \$(pwd)/cnv_output/normal -name "*_cnvs.tsv" 2>/dev/null | wc -l

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clincnv: \$(git -C /opt/clincnv log --oneline -1 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
