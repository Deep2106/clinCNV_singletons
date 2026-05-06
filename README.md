# ClinCNV WES CNV Calling Pipeline (singletons)

A production-grade DSL2 Nextflow pipeline for germline CNV calling from whole-exome sequencing (WES) data using [ClinCNV](https://github.com/imgag/ClinCNV), with optional VEP annotation and gene panel filtering.

## Pipeline Overview

```
samplesheet.csv
    │
    ├── GC_ANNOTATE_BED   (ngs-bits BedAnnotateGC + BedSort  once per run)
    ├── INFER_SEX         (ngs-bits SampleGender -method xy  per sample)
    ├── BED_COVERAGE      (ngs-bits BedCoverage  per sample, scatter)
    ├── ARCHIVE_COVERAGE  (persist .cov files for incremental runs)
    ├── MERGE_COVERAGE    (ClinCNV mergeFilesFromFolder  full cohort gather)
    │
    ├── CLINCNV           (cohort-level germline CNV calling)
    │
    ├── MAKE_VCF          (ClinCNV TSV → VCF with correct CN logic  per sample)
    ├── SORT_VCF          (bcftools sort)
    │
    ├── [optional] VEP_ANNOTATE    (StructuralVariantOverlap + MANE Select)
    └── [optional] RSCRIPT_FILTER  (gene panel filter)
```

### Outputs

| Directory | Contents | When |
|---|---|---|
| `results/1_clincnv_output/` | Coverage files, ClinCNV calls, sorted VCFs | Always |
| `results/2_vep_annotated/` | VEP-annotated VCFs | `--run_vep true` |
| `results/3_gene_panel_filtered/` | Gene panel filtered TSVs | `--run_vep true --run_rfilter true` |

---

## Requirements

### Software
- [Nextflow](https://nextflow.io/) ≥ 25.04.0
- [Singularity](https://sylabs.io/singularity/) ≥ 3.8
- Java 17 or greater

### Containers (Singularity `.sif`)

| Tool | Container | Build |
|---|---|---|
| ngs-bits 2025_12 | `ngsbits.sif` | `singularity build ngsbits.sif containers/ngsbits.def` |
| ClinCNV (master) + R 4.2 | `clincnv.sif` | `singularity build clincnv.sif containers/clincnv.def` |
| BCFtools 1.23 | `bcftool_1.23.sif` | Pre-built |
| VEP 115 | `vep.sif` | Pre-built |

---

## Quick Start

### 1. Clone the pipeline
```bash
git clone https://github.com/Deep2106/clinCNV_singletons.git
cd cnv-pipeline
```

### 2. Prepare samplesheet
```csv
sample_id,bam,bai,sex
SAMPLE_001,/data/bams/SAMPLE_001.markdup.recal.bam,/data/bams/SAMPLE_001.markdup.recal.bai,1
SAMPLE_002,/data/bams/SAMPLE_002.markdup.recal.bam,/data/bams/SAMPLE_002.markdup.recal.bai,2
SAMPLE_003,/data/bams/SAMPLE_003.markdup.recal.bam,/data/bams/SAMPLE_003.markdup.recal.bai,0
```

**Sex codes:** `1` = male (XY), `2` = female (XX), `0` = unknown (inferred automatically from BAM via `SampleGender`)

### 3. Edit and submit
```bash
# Edit paths in run_pipeline.slurm, then:
sbatch run_pipeline.slurm
```

---

## Samplesheet Format

| Column | Required | Description |
|---|---|---|
| `sample_id` | Yes | Unique sample identifier. Must match BAM basename convention |
| `bam` | Yes | Absolute path to BAM file |
| `bai` | Yes | Absolute path to BAM index |
| `sex` | No | 1=male, 2=female, 0/blank=unknown (auto-inferred) |

---


## CNV Target BED

Use the **CNV-specific** target BED (not the SNV padded BED) :

```
chr1    65564   65573   
chr1    69036   70008   
chr1    358066  358183  
chr1    370137  370737 
```

> **Do NOT use SNV padded BED**  padding dilutes coverage signal and inflates false positive CNV calls.

---

## Parameters

### Required
```bash
--input          # samplesheet CSV
--fasta          # hg38 reference FASTA
--fasta_fai      # hg38 FASTA index
--cnv_bed        # CNV target BED with/without gene names in col 4
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
                    # null = first run (default)
```

### ClinCNV Tuning
```bash
--clincnv_score_threshold    # minimum score (default: 60)
--clincnv_min_length         # minimum targets per CNV (default: 0)
--clincnv_min_cluster        # minimum cluster size (default: 10, raise to 15 at >100 samples)
--clincnv_max_germ_cnvs      # max CNVs per sample (default: 2000)
```

---

## Incremental Cohort Design

The pipeline supports incremental addition of samples across batches:

**First run:**
```bash
PREVIOUS_COV_DIR=""   # leave empty
```

**Subsequent runs:**
```bash
PREVIOUS_COV_DIR="/data/project/coverage_archive"
```

Each run archives current batch `.cov` files and merges with all previous batches for ClinCNV normalisation. ClinCNV re-processes the full cohort every batch.

> **Note:** Always add complete family units per batch  never split trios/duos across batches.

---

## Example Run

### First batch (50 singletons)
```bash
nextflow run main.nf \
    -profile slurm \
    --input              /data/project/samplesheet_batch1.csv \
    --fasta              /data/reference/hg38.fasta \
    --fasta_fai          /data/reference/hg38.fasta.fai \
    --cnv_bed            /data/reference/cnv_targets.bed \
    --outdir             /data/project/results \
    --run_vep            true \
    --vep_cache          /data/db/vep_cache \
    --vep_plugins        /data/db/VEP_plugins \
    --sv_overlap_vcf     /data/db/sv_cache/gnomad.v4.1.sv.sites.vcf.gz \
    --run_rfilter        true \
    --gene_panel         /data/reference/gene_panel.xlsx \
    -resume
```

### Second batch (add 20 samples)
```bash
nextflow run main.nf \
    -profile slurm \
    --input              /data/project/samplesheet_batch2.csv \
    --fasta              /data/reference/hg38.fasta \
    --fasta_fai          /data/reference/hg38.fasta.fai \
    --cnv_bed            /data/reference/cnv_targets.bed \
    --outdir             /data/project/results \
    --previous_cov_dir   /data/project/coverage_archive \
    --run_vep            true \
    --vep_cache          /data/db/vep_cache \
    --vep_plugins        /data/db/VEP_plugins \
    --sv_overlap_vcf     /data/db/sv_cache/gnomad.v4.1.sv.sites.vcf.gz \
    --run_rfilter        true \
    --gene_panel         /data/reference/gene_panel.xlsx \
    -resume
```

---

## ClinCNV Parameters Guide

| Parameter | Value | Rationale |
|---|---|---|
| `--scoreG 60` | Clinical grade threshold | Lower (50) for discovery, higher (80) for clinical |
| `--lengthG 0` | Allow single-exon CNVs | Required for renal disease (e.g. PKD1, COL4A5) |
| `--minimumNumOfElemsInCluster 10` | Min cluster size | Raise to 15 at >100 samples |
| `--maxNumGermCNVs 2000` | Permissive for affected cohort | Lower (500) for healthy reference panels |
| `--noPlot` | Disable R plots | Prevents write failures in read-only containers |

---

## VEP Annotation

- **Mode:** `--offline --cache` (no database connection required)
- **Cache:** homo_sapiens_vep_115_GRCh38
- **Plugin:** `StructuralVariantOverlap` with gnomAD-SV v4.1 (GRCh38 native)
- **Flags:** `--mane_select` (canonical transcript per gene)
- **No cache download needed** if running VEP for SV-only  only `StructuralVariantOverlap` plugin is used

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
cnv-pipeline/
├── main.nf                          # Pipeline entry point
├── nextflow.config                  # Parameters and container paths
├── run_pipeline.slurm               # SLURM submission script
├── assets/
│   ├── header.vcf                   # VCF header template
│   └── NO_FILE                      # Placeholder for optional inputs
├── bin/
│   └── vep_filter_1000gene.R        # Gene panel filter script
├── conf/
│   ├── base.config                  # Resource labels
│   ├── slurm.config                 # SLURM executor settings
│   ├── modules.config               # Per-process publishDir
│   └── test.config                  # Test profile
├── containers/
│   ├── clincnv.def                  # Singularity definition for ClinCNV
│   └── ngsbits.def                  # Singularity definition for ngs-bits
├── modules/local/
│   ├── gc_annotate_bed/             # BedAnnotateGC + BedSort
│   ├── bed_coverage/                # BedCoverage per sample
│   ├── merge_coverage/              # Cohort coverage merge
│   ├── archive_coverage/            # Persistent .cov store
│   ├── clincnv/                     # ClinCNV cohort calling
│   ├── infer_sex/                   # SampleGender sex inference
│   ├── make_vcf/                    # TSV → VCF conversion
│   ├── sort_vcf/                    # bcftools sort
│   ├── vep_annotate/                # VEP annotation
│   └── rscript_filter/              # Gene panel filter
├── scripts/
│   ├── 1_make_vcf.slurm             # Standalone: TSV → VCF
│   ├── 2_vep_annotate.slurm         # Standalone: VEP annotation
│   └── 3_rfilter.slurm              # Standalone: R filter
└── subworkflows/local/
    ├── coverage_workflow/           # Steps 1-3: coverage calculation
    ├── cnv_calling_workflow/        # Steps 4-6: CNV calling + VCF
    └── annotation_workflow/         # Steps 7-8: VEP + filter
```

---

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `Input BED file has to be sorted` | ngs-bits 2025_12 requires sorted input | `BedSort` added to pipeline  ensure latest `gc_annotate_bed.nf` |
| `failed to map segment from shared object` | `/tmp` is noexec on compute nodes | `TMPDIR=$(pwd)/tmp` set in `clincnv.nf` |
| `cannot open the connection` in ClinCNV | Container filesystem read-only | ClinCNV copied to work dir  ensure latest `clincnv.nf` |
| `bcftools sort` mkdtemp fails | `/local/scratch` noexec | `-T .` flag in `sort_vcf.nf` |
| VEP cache not found | Cache path wrong or not extracted | Extract `homo_sapiens_vep_115_GRCh38.tar.gz` |
| R `library(tidyverse)` not found | Not installed in container | Use individual packages: `dplyr`, `tidyr`, `readr`, `stringr` |

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

Deepak Bharti 

Clinical Bioinformatician 

RCSI, Dublin

---

## License & Citation

This pipeline is freely available for use and modification under the [APACHE License](LICENSE).

**If you use or modify this pipeline in your work, please cite:**

> Bharti D. *ClinCNV WES CNV Calling Pipeline*. RCSI, FutureNeuro Group (2026).
> Available at: https://github.com/Deep2106/clinCNV_singletons.git

You may also include the following acknowledgement in your methods section:

> CNV calling was performed using a DSL2 Nextflow pipeline based on ClinCNV
> (Bessonov et al. 2023), developed at RCSI (Bharti D, 2026).
> VEP annotation used StructuralVariantOverlap with gnomAD-SV v4.1.
> Sex was inferred from BAM files using ngs-bits SampleGender.
