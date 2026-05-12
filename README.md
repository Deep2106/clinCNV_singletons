# ClinCNV WES CNV Calling Pipeline  Singletons + Duos

A production-grade DSL2 Nextflow pipeline for germline CNV calling from WES data using [ClinCNV](https://github.com/imgag/ClinCNV), with optional VEP annotation and gene panel filtering.

> **For trios/duos with de novo calling, use the [clinCNV_trios](https://github.com/Deep2106/clinCNV_trios) repository instead.**

---

## Pipeline Overview

```
samplesheet.csv (sample_id, bam, bai, sex)
    │
    ├── SAMPLESHEET_CHECK   (validate CSV, auto-detect cohort type)
    ├── GC_ANNOTATE_BED     (ngs-bits BedAnnotateGC + BedSort  once per run)
    ├── INFER_SEX           (ngs-bits SampleGender -method xy  per sample)
    ├── BED_COVERAGE        (ngs-bits BedCoverage  per sample, scatter)
    ├── ARCHIVE_COVERAGE    (persist .cov files for incremental runs)
    ├── MERGE_COVERAGE      (ClinCNV mergeFilesFromFolder  full cohort gather)
    │
    ├── CLINCNV             (cohort-level germline CNV calling  GERMLINE_SINGLE)
    │
    ├── MAKE_VCF            (ClinCNV TSV → VCF with correct CN logic  per sample)
    ├── SORT_VCF            (bcftools sort)
    │
    ├── [optional] VEP_ANNOTATE    (StructuralVariantOverlap + MANE Select)
    └── [optional] RSCRIPT_FILTER  (gene panel filter)
```

### Outputs

| Directory | Contents | When |
|---|---|---|
| `results/1_clincnv_output/` | Coverage files, ClinCNV calls (`_cnvs.tsv`, `_cnvs.seg`, `_cov.seg`), sorted VCFs | Always |
| `results/2_vep_annotated/` | VEP-annotated VCFs | `--run_vep true` |
| `results/3_gene_panel_filtered/` | Gene panel filtered TSVs | `--run_vep true --run_rfilter true` |

---

## Requirements

### Software
- Nextflow ≥ 25.04.0
- Singularity ≥ 3.8
- Java 17+

### Containers (Singularity `.sif`)

| Tool | Container |
|---|---|
| Python 3.14 | `python_3.14.2.sif` |
| ngs-bits 2025_12 | `ngsbits.sif` |
| ClinCNV (master) + R 4.2 | `clincnv.sif` |
| BCFtools 1.23 | `bcftool_1.23.sif` |
| VEP 115 | `vep.sif` |

---

## Quick Start

### 1. Prepare samplesheet
```csv
sample_id,bam,bai,sex
SAMPLE_001,/data/bams/SAMPLE_001.markdup.recal.bam,/data/bams/SAMPLE_001.markdup.recal.bai,1
SAMPLE_002,/data/bams/SAMPLE_002.markdup.recal.bam,/data/bams/SAMPLE_002.markdup.recal.bai,2
SAMPLE_003,/data/bams/SAMPLE_003.markdup.recal.bam,/data/bams/SAMPLE_003.markdup.recal.bai,0
```

**Sex codes:** `1`=male, `2`=female, `0`=unknown (auto-inferred from BAM)

### 2. Edit and submit
```bash
# Edit paths in run_pipeline.slurm, then:
sbatch run_pipeline.slurm
```

---

## Samplesheet Format

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Unique sample identifier |
| `bam` | Yes | Absolute path to BAM file |
| `bai` | Yes | Absolute path to BAM index |
| `sex` | No | 1=male, 2=female, 0/blank=unknown |

---

## Parameters

### Required
```bash
--input          # samplesheet CSV
--fasta          # hg38 reference FASTA
--fasta_fai      # hg38 FASTA index
--cnv_bed        # CNV target BED (see note below)
--outdir         # output directory
```

### Optional  VEP Annotation
```bash
--run_vep        # true/false (default: false)
--vep_cache      # path to VEP cache (homo_sapiens_vep_115_GRCh38)
--vep_plugins    # path to VEP plugins directory
--sv_overlap_vcf # gnomAD-SV v4.1 VCF (tabix-indexed)
```

### Optional  Gene Panel Filter
```bash
--run_rfilter    # true/false (default: false, requires run_vep=true)
--gene_panel     # XLSX file with "Gene Symbol" column
```

### Optional  Incremental Runs
```bash
--previous_cov_dir  # NFS path to archived .cov files from previous batches
```

### ClinCNV Tuning
```bash
--clincnv_score_threshold    # default: 60
--clincnv_min_length         # default: 0 (allow single-exon CNVs)
--clincnv_min_cluster        # default: 10 (raise to 15 at >100 samples)
--clincnv_max_germ_cnvs      # default: 2000
```

---

## CNV Target BED

Use the **CNV-specific** target BED  not the SNV padded BED. Padding dilutes coverage signal and inflates false positive CNV calls.

---

## Incremental Cohort Design

**First run:** leave `PREVIOUS_COV_DIR=""` in `run_pipeline.slurm`  coverage files are archived to `results/coverage_archive/` automatically.

**Subsequent runs:** set `PREVIOUS_COV_DIR="${OUTDIR}/coverage_archive"`  new batch merges with all previous batches for ClinCNV normalisation.

> ClinCNV re-processes the full cohort every batch. Never split family units across batches.

---

## ClinCNV Parameters Guide

| Parameter | Default | Rationale |
|---|---|---|
| `--scoreG` | 60 | Clinical grade; lower (50) for discovery |
| `--lengthG` | 0 | Allow single-exon CNVs |
| `--minimumNumOfElemsInCluster` | 10 | Raise to 15 at >100 samples |
| `--maxNumGermCNVs` | 2000 | Permissive for affected cohorts |
| `--noPlot` | on | Prevents write failures in read-only containers |

---

## VEP Annotation

- **Mode:** `--offline --cache`
- **Cache:** homo_sapiens_vep_115_GRCh38
- **Assembly:** GRCh38 (explicit `--assembly GRCh38 --species homo_sapiens`)
- **Plugin:** `StructuralVariantOverlap` with gnomAD-SV v4.1

### Download VEP cache
```bash
singularity exec -B /data vep.sif \
    perl /opt/vep/src/ensembl-vep/INSTALL.pl \
    -a cf -s homo_sapiens -y GRCh38 \
    --CACHE_VERSION 115 \
    -c /data/db/vep_cache
```

### Download gnomAD-SV v4.1
```bash
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
```

---

## Pipeline Structure

```
clinCNV_singletons/
├── main.nf
├── nextflow.config
├── run_pipeline.slurm
├── assets/
│   ├── header.vcf
│   └── NO_FILE
├── bin/
│   └── vep_filter_1000gene.R
├── conf/
│   ├── base.config
│   ├── slurm.config
│   ├── modules.config
│   └── test.config
├── containers/
│   ├── clincnv.def
│   └── ngsbits.def
├── modules/local/
│   ├── samplesheet_check/
│   ├── gc_annotate_bed/
│   ├── bed_coverage/
│   ├── merge_coverage/
│   ├── archive_coverage/
│   ├── clincnv/
│   ├── infer_sex/
│   ├── make_vcf/
│   ├── sort_vcf/
│   ├── vep_annotate/
│   └── rscript_filter/
└── subworkflows/local/
    ├── input_validation/
    ├── coverage_workflow/
    ├── cnv_calling_workflow/
    └── annotation_workflow/
```

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| VEP cache not found | Wrong path or not extracted | Extract `homo_sapiens_vep_115_GRCh38.tar.gz` |
| `--cache_version null` VEP error | `vep_cache_version` param not set | Removed from `vep_annotate.nf`  VEP auto-detects |
| R `library(tidyverse)` not found | Not in container | Uses individual packages: `dplyr`, `tidyr`, `readr`, `stringr` |
| MAKE_VCF killed by SLURM | Transient node failure | Re-run with `-resume` |

---

## Citation

If you use this pipeline, please cite:

- **ClinCNV:** [Bessonov et al. 2023](https://github.com/imgag/ClinCNV)
- **ngs-bits:** [IMGAG](https://github.com/imgag/ngs-bits)
- **VEP:** [McLaren et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4)
- **gnomAD-SV:** [Collins et al. 2020](https://doi.org/10.1038/s41586-020-2287-8)
- **Nextflow:** [Di Tommaso et al. 2017](https://doi.org/10.1038/nbt.3820)

---

## Author

Deepak Bharti  Clinical Bioinformatician, RCSI Dublin

> CNV calling was performed using a DSL2 Nextflow pipeline based on ClinCNV (Bessonov et al. 2023), developed at RCSI (Bharti D, 2026). VEP annotation used StructuralVariantOverlap with gnomAD-SV v4.1. Sex was inferred from BAM files using ngs-bits SampleGender.
