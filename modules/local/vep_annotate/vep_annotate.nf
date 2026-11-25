process VEP_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    container params.containers.vep

    publishDir "${params.outdir}/2_vep_annotated", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path  vep_cache
    path  vep_sv_plugin            // patched StructuralVariantOverlap.pm (assets/vep_plugins/)
    path  sv_overlap_vcf
    path  sv_overlap_tbi

    output:
    tuple val(meta), path("${meta.id}.vep.vcf"),            emit: vcf
    tuple val(meta), path("${meta.id}.sv_overlap_qc.tsv"),  emit: sv_qc
    path "versions.yml",                                    emit: versions

    script:
    def strict = params.sv_overlap_strict ? 'true' : 'false'
    """
    # -- Plugin: staged as a single file so -resume detects patch changes ----
    mkdir -p vep_plugins_local
    cp ${vep_sv_plugin} vep_plugins_local/StructuralVariantOverlap.pm
    if ! grep -q 'PATCH: avoid chrchr' vep_plugins_local/StructuralVariantOverlap.pm; then
        echo "ERROR: StructuralVariantOverlap.pm is not the patched version (chr-prefix bug -> zero SV overlaps)" >&2
        exit 1
    fi

    vep \\
        --offline \\
        --cache \\
        --dir_cache ${vep_cache} \\
        --assembly GRCh38 \\
        --species homo_sapiens \\
        --format vcf \\
        --vcf \\
        --force_overwrite \\
        --no_stats \\
        --mane_select \\
        --dir_plugins vep_plugins_local \\
        --plugin StructuralVariantOverlap,file=${sv_overlap_vcf},percentage=${params.sv_overlap_pct},same_type=${params.sv_overlap_same_type},cols=${params.sv_overlap_cols} \\
        -i ${vcf} \\
        -o ${meta.id}.vep.vcf \\
        --fork ${task.cpus}

    # -- QC: fraction of CNVs with >=1 gnomAD-SV overlap ----------------------
    # Column located by name from the CSQ header (robust to VEP field order)
    awk -F'\\t' -v id="${meta.id}" '
        /^##INFO=<ID=CSQ/ {
            s = \$0; sub(/.*Format: /, "", s); sub(/">.*/, "", s)
            n = split(s, h, "|"); for (i = 1; i <= n; i++) if (h[i] == "SV_overlap_name") col = i
            next
        }
        /^#/ { next }
        {
            tot++
            if (col && match(\$8, /CSQ=[^;]+/)) {
                csq = substr(\$8, RSTART + 4, RLENGTH - 4)
                m = split(csq, tx, ",")
                for (j = 1; j <= m; j++) { split(tx[j], f, "|"); if (f[col] != "") { hit++; break } }
            }
        }
        END {
            if (!col) { print "ERROR: SV_overlap_name not in CSQ header" > "/dev/stderr"; exit 1 }
            printf "%s\\t%d\\t%d\\t%.3f\\n", id, tot, hit, (tot ? hit / tot : 0)
        }' ${meta.id}.vep.vcf > ${meta.id}.sv_overlap_qc.tsv

    read -r _ TOT HIT FRAC < ${meta.id}.sv_overlap_qc.tsv
    if [ "${strict}" = "true" ] && [ "\$TOT" -ge 10 ] && \\
       awk -v f="\$FRAC" -v m="${params.sv_overlap_min_frac}" 'BEGIN{exit !(f < m)}'; then
        echo "ERROR: only \$HIT/\$TOT CNVs have a gnomAD-SV overlap (< ${params.sv_overlap_min_frac})" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep "ensembl-vep" | awk '{print \$NF}')
        StructuralVariantOverlap: "patched (chr-prefix fix)"
    END_VERSIONS
    """
}
