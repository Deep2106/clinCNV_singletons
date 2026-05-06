#!/bin/bash
#SBATCH --job-name=cnv_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --partition=<PARTITION>
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# ============================================================================
# ClinCNV singletons WES CNV Calling Pipeline v1.0 - RCSI HPC
# ============================================================================
# Usage:
#   FIRST RUN:  sbatch run_pipeline.slurm
#   NEXT BATCH: Edit BATCH_NAME, INPUT_CSV, PED_FILE, PREVIOUS_COV_DIR
# ============================================================================

# ============================================================================
# USER CONFIGURATION  EDIT THESE FOR EACH RUN
# ============================================================================
module purge
module --ignore-cache load Java/17

BATCH_NAME="batch1"
INPUT_CSV="samplesheet.csv"

RUN_VEP="true"    # set false if not required
RUN_RFILTER="true" # set false if data is not from renal-kidney project

# PED file for trio/duo aware calling
# First run: leave empty if PED not yet available
# Subsequent runs: set full path e.g. /path/to/cohort.ped
PED_FILE=""

# Incremental coverage store
# First run: leave empty  uses current batch only
# Subsequent runs: set full path e.g. ${PROJECT_DIR}/coverage_archive
PREVIOUS_COV_DIR=""

# ============================================================================
# PATHS  do not edit unless directory structure changes
# ============================================================================

PROJECT_DIR=/path/to/clinCNV_2026
PIPELINE_DIR=${PROJECT_DIR}/clinCNV_singletons

export NXF_WORK=${PROJECT_DIR}/work_cnv
mkdir -p ${NXF_WORK}

OUTDIR=${PROJECT_DIR}/results_cnv

# ============================================================================
# SINGULARITY CONFIGURATION
# ============================================================================

export NXF_SINGULARITY_LIBRARYDIR=/path/to/Docker_Images
export NXF_SINGULARITY_CACHEDIR=${PROJECT_DIR}/.singularity_cache
mkdir -p ${NXF_SINGULARITY_CACHEDIR}

export NXF_ANSI_LOG=false
export _JAVA_OPTIONS="-Djava.io.tmpdir=${TMPDIR:-/tmp}"

# ============================================================================
# PRE-RUN CHECKS
# ============================================================================

cd ${PIPELINE_DIR} || { echo "ERROR: Cannot cd to ${PIPELINE_DIR}"; exit 1; }

if [ ! -f "${PROJECT_DIR}/${INPUT_CSV}" ]; then
    echo "ERROR: Samplesheet not found: ${PROJECT_DIR}/${INPUT_CSV}"
    exit 1
fi

if [ -n "${PED_FILE}" ] && [ ! -f "${PED_FILE}" ]; then
    echo "ERROR: PED file not found: ${PED_FILE}"
    exit 1
fi

echo "=============================================================================="
echo "ClinCNV WES CNV Calling Pipeline v1.0"
echo "=============================================================================="
echo "Date:              $(date)"
echo "Node:              $(hostname)"
echo "SLURM Job ID:      ${SLURM_JOB_ID}"
echo "Batch:             ${BATCH_NAME}"
echo "Input CSV:         ${INPUT_CSV}"
echo "PED file:          ${PED_FILE:-NOT PROVIDED (singleton/cohort mode)}"
echo "Previous cov dir:  ${PREVIOUS_COV_DIR:-NOT SET (first run)}"
echo "Run VEP:           ${RUN_VEP}"
echo "Run R filter:      ${RUN_RFILTER}"
echo "Output:            ${OUTDIR}"
echo "NXF_WORK:          ${NXF_WORK}"
echo "=============================================================================="

# ============================================================================
# BUILD OPTIONAL ARGUMENTS
# Only pass flags when values are set  avoids passing empty strings to Nextflow
# ============================================================================

OPTIONAL_ARGS=""

if [ -n "${PED_FILE}" ]; then
    OPTIONAL_ARGS="${OPTIONAL_ARGS} --ped_file ${PED_FILE}"
fi

if [ -n "${PREVIOUS_COV_DIR}" ]; then
    OPTIONAL_ARGS="${OPTIONAL_ARGS} --previous_cov_dir ${PREVIOUS_COV_DIR}"
fi

# ============================================================================
# RUN PIPELINE
# ============================================================================

nextflow run main.nf \
    -profile slurm \
    -w ${NXF_WORK} \
    --input          ${PROJECT_DIR}/${INPUT_CSV} \
    --fasta          /path/to/reference/hg38.fasta \
    --fasta_fai      /path/to/reference/hg38.fasta.fai \
    --cnv_bed        /path/to/reference/hg38/S000021_hg38_targets.bed \
    --vep_plugins    /path/to/resources/VEP_plugins \
    --sv_overlap_vcf /path/to/resources/sv_cache/gnomad.v4.1.sv.sites.vcf.gz \
    --gene_panel     /path/to/Renal_superpanel_broad.xlsx \
    --outdir         ${OUTDIR} \
    --run_vep        ${RUN_VEP} \
    --run_rfilter    ${RUN_RFILTER} \
    ${OPTIONAL_ARGS} \
    --vep_cache /path/to/resources/vep_cache \
    -resume

EXIT_CODE=$?

echo "=============================================================================="
echo "Pipeline exit code: ${EXIT_CODE}"
echo "Date: $(date)"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS"
    echo "  1_clincnv_output:       ${OUTDIR}/1_clincnv_output/"
    echo "  2_vep_annotated:        ${OUTDIR}/2_vep_annotated/  (if run_vep=true)"
    echo "  3_gene_panel_filtered:  ${OUTDIR}/3_gene_panel_filtered/  (if run_rfilter=true)"
else
    echo "FAILED  check: nextflow_${SLURM_JOB_ID}.err and .nextflow.log"
fi
echo "=============================================================================="

exit ${EXIT_CODE}
