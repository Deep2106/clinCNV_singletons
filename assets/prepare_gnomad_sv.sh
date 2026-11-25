#!/bin/bash
# =============================================================================
# prepare_gnomad_sv.sh - build the gnomAD-SV reference used by VEP
#                        StructuralVariantOverlap (params.sv_overlap_vcf)
# =============================================================================
# Input : gnomad.v4.1.sv.sites.vcf.gz (+ .tbi)
# Output: gnomad.v4.1.sv.sites.DELDUP.noME.vcf.gz (+ .tbi)
#
# What it does (and why):
#   1. Keeps plain <DEL> / <DUP> only
#      - drops CPX/BND/INV/INS (never comparable to exome CNVs)
#      - drops mobile-element deletions (<DEL:ME:*>): single-TE events, mostly
#        UNRESOLVED, that otherwise lend high AFs to unrelated exonic CNVs
#   2. Keeps ALL FILTER values (PASS-only removes most common segdup CNVs)
#   3. Copies FILTER into INFO/GNOMAD_FILTER (';' -> '+') so VEP can report it
#      via cols=AF:GNOMAD_FILTER and the report can show PASS-only vs all AF
#
# Usage: bash prepare_gnomad_sv.sh /path/to/gnomad.v4.1.sv.sites.vcf.gz [outdir]
# Needs: bcftools, tabix, bgzip
# =============================================================================
set -euo pipefail

SRC="${1:?Usage: prepare_gnomad_sv.sh SRC.vcf.gz [OUTDIR]}"
OUTDIR="${2:-$(dirname "${SRC}")}"
OUT="${OUTDIR}/gnomad.v4.1.sv.sites.DELDUP.noME.vcf.gz"
TMP="$(mktemp -d "${OUTDIR}/.prep_gnomad_sv.XXXXXX")"
trap 'rm -rf "${TMP}"' EXIT

for t in bcftools tabix bgzip; do command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done
[ -f "${SRC}.tbi" ] || { echo "ERROR: ${SRC}.tbi missing" >&2; exit 1; }

echo "[1/3] Filter plain DEL/DUP"
bcftools view -i '(INFO/SVTYPE="DEL" && ALT="<DEL>") || (INFO/SVTYPE="DUP" && ALT="<DUP>")' \
    "${SRC}" -Oz -o "${TMP}/deldup.vcf.gz"
tabix -p vcf "${TMP}/deldup.vcf.gz"

echo "[2/3] FILTER -> INFO/GNOMAD_FILTER"
bcftools query -f '%CHROM\t%POS\t%ID\t%FILTER\n' "${TMP}/deldup.vcf.gz" \
    | awk 'BEGIN{OFS="\t"}{gsub(/;/,"+",$4); print}' | bgzip > "${TMP}/filt.tsv.gz"
tabix -s1 -b2 -e2 "${TMP}/filt.tsv.gz"
echo '##INFO=<ID=GNOMAD_FILTER,Number=1,Type=String,Description="gnomAD-SV v4.1 FILTER value (; replaced by +)">' > "${TMP}/hdr.txt"

bcftools annotate -a "${TMP}/filt.tsv.gz" -c CHROM,POS,ID,INFO/GNOMAD_FILTER -h "${TMP}/hdr.txt" \
    "${TMP}/deldup.vcf.gz" -Oz -o "${OUT}"
tabix -f -p vcf "${OUT}"

echo "[3/3] Summary"
N_IN=$(bcftools view -H "${TMP}/deldup.vcf.gz" | wc -l)
N_TAG=$(bcftools view -H -i 'INFO/GNOMAD_FILTER!="."' "${OUT}" | wc -l)
echo "  DEL/DUP records : ${N_IN}"
echo "  with GNOMAD_FILTER: ${N_TAG}"
[ "${N_IN}" -eq "${N_TAG}" ] || echo "WARNING: not every record got GNOMAD_FILTER" >&2
bcftools query -f '%INFO/GNOMAD_FILTER\n' "${OUT}" | sort | uniq -c | sort -k1,1nr | head
echo "Done: ${OUT}"
