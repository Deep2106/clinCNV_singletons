process MAKE_VCF {
    tag "${meta.id}"
    label 'process_single'

    container params.containers.bcftools

    input:
    tuple val(meta), path(cnv_tsv)
    path  vcf_header
    path  fasta_fai   // hg38.fasta.fai  staged alongside process for contig lines

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path  "versions.yml",                    emit: versions

    script:
    def sex = meta.sex ?: 0
    """
    # Build VCF header: fileformat → contigs from FAI → meta lines → #CHROM
    sed "s/SAMPLE_ID/${meta.id}/g" ${vcf_header} > header.tmp
    grep "^##fileformat" header.tmp > ${meta.id}.vcf
    awk '{print "##contig=<ID="\$1",length="\$2">"}' ${fasta_fai} >> ${meta.id}.vcf
    grep "^##" header.tmp | grep -v "^##fileformat" >> ${meta.id}.vcf
    grep "^#CHROM" header.tmp >> ${meta.id}.vcf
    rm header.tmp

    # Convert ClinCNV TSV to VCF
    grep -v '^#' ${cnv_tsv} | awk -v sex=${sex} '
    BEGIN { OFS="\\t" }
    {
        chr=\$1; start=\$2; end=\$3; cn=int(\$4); ll=\$5+0
        qual=(ll>0) ? int(ll) : 0
        id="CNV_" chr "_" start "_" end
        info="END=" end
        alt=""; gt=""

        if      (cn==0) { alt="<DEL>"; gt="1/1" }
        else if (cn==1) { alt="<DEL>"; gt="0/1" }
        else if (cn==2) {
            if (chr=="chrX" || chr=="chrY") {
                if (sex==1) { alt="<DUP>"; gt="0/1" }
                else next
            } else next
        }
        else if (cn==3) { alt="<DUP>"; gt="0/1" }
        else if (cn==4) { alt="<DUP>"; gt="1/1" }
        else if (cn>4)  { alt="<DUP>"; gt="./." }

        if (alt!="") print chr, start, id, "N", alt, qual, "PASS", info, "GT:CN", gt ":" cn
    }
    ' >> ${meta.id}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
