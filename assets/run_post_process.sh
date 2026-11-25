#!/bin/bash
# =============================================================================
# run_post_process.sh  Batch wrapper for post_process_cnv.py
# =============================================================================
# Processes VEP-annotated VCFs to Excel reports incrementally.
# Only processes NEW samples - skips samples that already have an Excel.
# Use --force to reprocess all samples regardless.
#
# Usage:
#   # Singleton / Duo batch:
#   bash run_post_process.sh --vcf_dir DIR --outdir DIR [--python /env/bin/python]
#
#   # Trio batch:
#   bash run_post_process.sh --vcf_dir DIR --cnvs_dir DIR --outdir DIR
#
#   # Force reprocess all (e.g. after VEP re-run):
#   bash run_post_process.sh --vcf_dir DIR --outdir DIR --force
#
# Arguments:
#   --vcf_dir    Directory containing *.vep.vcf files         [required]
#   --outdir     Output directory for Excel and BED files      [required]
#   --cnvs_dir   ClinCNV trios output dir (trios mode only)    [optional]
#   --script     Path to post_process_cnv.py                   [optional, default: same dir as this script]
#   --ped        PED file (proband sex for trio chrX/Y inheritance) [optional]
#   --common_af  gnomAD-SV AF threshold for 'common' flag       [optional, default 0.01]
#   --python     Python interpreter to use                     [optional, default: python3]
#   --force      Reprocess all samples including existing ones [optional flag]
# =============================================================================

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

# -- Defaults -----------------------------------------------------------------
VCF_DIR=""
CNVS_DIR=""
OUTDIR=""
PED=""
COMMON_AF=""
FORCE=false
PYTHON="python3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/post_process_cnv.py"

# -- Parse arguments ----------------------------------------------------------
need_val() { [ $# -ge 2 ] && [ -n "$2" ] || die "$1 requires a value"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf_dir)  need_val "$@"; VCF_DIR="$2";   shift 2 ;;
        --cnvs_dir) need_val "$@"; CNVS_DIR="$2";  shift 2 ;;
        --outdir)   need_val "$@"; OUTDIR="$2";    shift 2 ;;
        --script)   need_val "$@"; PY_SCRIPT="$2"; shift 2 ;;
        --python)   need_val "$@"; PYTHON="$2";    shift 2 ;;
        --ped)      need_val "$@"; PED="$2";       shift 2 ;;
        --common_af) need_val "$@"; COMMON_AF="$2"; shift 2 ;;
        --force)    FORCE=true;                    shift   ;;
        *)
            die "Unknown argument: $1
Usage: bash run_post_process.sh --vcf_dir DIR --outdir DIR [--cnvs_dir DIR] [--ped FILE] [--common_af F] [--python BIN] [--force]"
            ;;
    esac
done

# -- Validate -----------------------------------------------------------------
[ -n "${VCF_DIR}" ]   || die "--vcf_dir is required"
[ -n "${OUTDIR}" ]    || die "--outdir is required"
[ -d "${VCF_DIR}" ]   || die "VCF directory not found: ${VCF_DIR}"
[ -z "${PED}" ] || [ -f "${PED}" ] || die "PED file not found: ${PED}"
[ -f "${PY_SCRIPT}" ] || die "post_process_cnv.py not found at: ${PY_SCRIPT} (use --script)"

PYTHON_PATH="$(command -v "${PYTHON}" || true)"
[ -n "${PYTHON_PATH}" ] || die "Python interpreter not found: ${PYTHON}"

# Do NOT suppress stderr here - the real ImportError is the useful diagnostic
if ! "${PYTHON}" -c "import openpyxl"; then
    die "openpyxl not importable by ${PYTHON_PATH} (wrong interpreter/env?)"
fi

mkdir -p "${OUTDIR}"

# -- Mode detection -----------------------------------------------------------
if [ -n "${CNVS_DIR}" ]; then
    [ -d "${CNVS_DIR}" ] || die "cnvs_dir not found: ${CNVS_DIR}"
    MODE="trio"
else
    MODE="singleton"
fi

# -- Log file -----------------------------------------------------------------
LOG="${OUTDIR}/post_process_$(date +%Y%m%d_%H%M%S).log"

log() { echo "$*" | tee -a "${LOG}"; }

log "============================================================"
log "ClinCNV Post-Processing Batch Wrapper"
log "============================================================"
log "Date        : $(date)"
log "Mode        : ${MODE}"
log "Python      : ${PYTHON_PATH}"
log "VCF dir     : ${VCF_DIR}"
log "Output dir  : ${OUTDIR}"
[ -n "${CNVS_DIR}" ] && log "CNVs dir    : ${CNVS_DIR}"
log "PED         : ${PED:-none}"
log "Common AF   : ${COMMON_AF:-0.01 (default)}"
log "Force mode  : ${FORCE}"
log "Log file    : ${LOG}"
log "============================================================"

# -- Find all VCF files -------------------------------------------------------
mapfile -t VCF_FILES < <(find "${VCF_DIR}" -maxdepth 1 -name "*.vep.vcf" | sort)

TOTAL=${#VCF_FILES[@]}
[ "${TOTAL}" -gt 0 ] || { log "ERROR: No *.vep.vcf files found in ${VCF_DIR}"; exit 1; }

log "Found ${TOTAL} VCF file(s)"
log ""

# -- Trio helper: proband = first element of CHILD-MOTHER-FATHER dir name -----
find_cnvs_txt() {
    local sample="$1" fam_dir fam hit
    for fam_dir in "${CNVS_DIR}"/*/; do
        [ -d "${fam_dir}" ] || continue
        fam="$(basename "${fam_dir}")"
        [ "${fam%%-*}" = "${sample}" ] || continue
        hit="$(find "${fam_dir}" -maxdepth 1 -name "*_CNVs.txt" -print -quit)"
        if [ -n "${hit}" ]; then
            echo "${hit}"
            return 0
        fi
    done
    return 0
}

# -- Counters -----------------------------------------------------------------
COUNT_PROCESSED=0
COUNT_SKIPPED=0
COUNT_FAILED=0
FAILED_SAMPLES=()

# -- Process each sample ------------------------------------------------------
for VCF in "${VCF_FILES[@]}"; do

    SAMPLE="$(basename "${VCF}" .vep.vcf)"
    EXCEL="${OUTDIR}/${SAMPLE}_CNV_report.xlsx"

    if [ "${FORCE}" = false ] && [ -f "${EXCEL}" ]; then
        log "[SKIP]    ${SAMPLE}  Excel already exists"
        COUNT_SKIPPED=$((COUNT_SKIPPED + 1))
        continue
    fi

    CMD=("${PYTHON}" "${PY_SCRIPT}" --vcf "${VCF}" --outdir "${OUTDIR}")
    [ -n "${PED}" ]       && CMD+=(--ped "${PED}")
    [ -n "${COMMON_AF}" ] && CMD+=(--common_af "${COMMON_AF}")

    if [ "${MODE}" = "trio" ]; then
        CNVS_TXT="$(find_cnvs_txt "${SAMPLE}")"
        if [ -n "${CNVS_TXT}" ]; then
            CMD+=(--cnvs_txt "${CNVS_TXT}")
            log "[PROCESS] ${SAMPLE}  trio mode (${CNVS_TXT})"
        else
            log "[PROCESS] ${SAMPLE}  singleton mode (no _CNVs.txt found)"
        fi
    else
        log "[PROCESS] ${SAMPLE}  singleton mode"
    fi

    if "${CMD[@]}" >> "${LOG}" 2>&1; then
        COUNT_PROCESSED=$((COUNT_PROCESSED + 1))
    else
        log "[FAILED]  ${SAMPLE}  check log: ${LOG}"
        COUNT_FAILED=$((COUNT_FAILED + 1))
        FAILED_SAMPLES+=("${SAMPLE}")
    fi

done

# -- Summary ------------------------------------------------------------------
log ""
log "============================================================"
log "Summary"
log "============================================================"
log "Total found    : ${TOTAL}"
log "Processed      : ${COUNT_PROCESSED}"
log "Skipped        : ${COUNT_SKIPPED}  (already have Excel; use --force to reprocess)"
log "Failed         : ${COUNT_FAILED}"

if [ "${COUNT_FAILED}" -gt 0 ]; then
    log ""
    log "Failed samples:"
    for s in "${FAILED_SAMPLES[@]}"; do
        log "  - ${s}"
    done
fi

log "============================================================"
log "Output dir: ${OUTDIR}"
log "Log file  : ${LOG}"
log "============================================================"

[ "${COUNT_FAILED}" -eq 0 ] || exit 1
exit 0
