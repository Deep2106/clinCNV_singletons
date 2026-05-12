#!/bin/bash
#SBATCH --job-name=cnv_singleton
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --partition=<PARTITION>
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# ============================================================================
# ClinCNV WES CNV Calling Pipeline - Singletons + Duos
# ============================================================================
# Samplesheet format: sample_id,bam,bai,sex
# For trios use clinCNV_trios pipeline instead.
#
# Usage:
#   FIRST RUN:  sbatch run_pipeline.slurm
#   NEXT BATCH: Edit BATCH_NAME, INPUT_CSV, PREVIOUS_COV_DIR
# ============================================================================

# ============================================================================
# USER CONFIGURATION - EDIT THESE FOR EACH RUN
# ============================================================================
module purge
module --ignore-cache load Java/17

BATCH_NAME="renal_test_CNV"
INPUT_CSV="samplesheet.csv"

RUN_VEP="true"
RUN_RFILTER="true"   # set true only when --run_vep is also true

# Gene panel XLSX - required when RUN_RFILTER=true
GENE_PANEL="/path/to/Renal_superpanel_broad.xlsx"

# Incremental coverage store
# First run:       leave empty - uses current batch only
# Subsequent runs: set full path e.g. ${PROJECT_DIR}/results_singletons/coverage_archive
PREVIOUS_COV_DIR=""

# ============================================================================
# PATHS - do not edit unless directory structure changes
# ============================================================================

PROJECT_DIR=/path/to/clinCNV_2026
PIPELINE_DIR=${PROJECT_DIR}/clinCNV_singletons

export NXF_WORK=${PROJECT_DIR}/work_singletons
mkdir -p ${NXF_WORK}

OUTDIR=${PROJECT_DIR}/results_singletons

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

if [ "${RUN_RFILTER}" == "true" ] && [ "${RUN_VEP}" != "true" ]; then
    echo "ERROR: RUN_RFILTER=true requires RUN_VEP=true"
    exit 1
fi

if [ "${RUN_RFILTER}" == "true" ] && [ ! -f "${GENE_PANEL}" ]; then
    echo "ERROR: Gene panel not found: ${GENE_PANEL}"
    exit 1
fi

echo "=============================================================================="
echo "ClinCNV WES CNV Calling Pipeline - Singletons + Duos"
echo "=============================================================================="
echo "Date:              $(date)"
echo "Node:              $(hostname)"
echo "SLURM Job ID:      ${SLURM_JOB_ID}"
echo "Batch:             ${BATCH_NAME}"
echo "Input CSV:         ${INPUT_CSV}"
echo "Previous cov dir:  ${PREVIOUS_COV_DIR:-NOT SET (first run)}"
echo "Run VEP:           ${RUN_VEP}"
echo "Run R filter:      ${RUN_RFILTER}"
echo "Gene panel:        ${RUN_RFILTER:+${GENE_PANEL}}"
echo "Output:            ${OUTDIR}"
echo "NXF_WORK:          ${NXF_WORK}"
echo "NOTE: For trios use clinCNV_trios pipeline"
echo "=============================================================================="

# ============================================================================
# BUILD OPTIONAL ARGUMENTS
# ============================================================================

OPTIONAL_ARGS=""

if [ -n "${PREVIOUS_COV_DIR}" ]; then
    OPTIONAL_ARGS="${OPTIONAL_ARGS} --previous_cov_dir ${PREVIOUS_COV_DIR}"
fi

if [ "${RUN_RFILTER}" == "true" ]; then
    OPTIONAL_ARGS="${OPTIONAL_ARGS} --gene_panel ${GENE_PANEL}"
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
    --cnv_bed        /path/to/hg38/targets.bed \
    --vep_cache      /path/to/resources/vep_cache \
    --vep_plugins    /path/to/resources/VEP_plugins \
    --sv_overlap_vcf /path/to/resources/sv_cache/gnomad.v4.1.sv.sites.vcf.gz \
    --outdir         ${OUTDIR} \
    --run_vep        ${RUN_VEP} \
    --run_rfilter    ${RUN_RFILTER} \
    ${OPTIONAL_ARGS} \
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
    echo "FAILED - check: nextflow_${SLURM_JOB_ID}.err and .nextflow.log"
fi
echo "=============================================================================="

exit ${EXIT_CODE}
