#!/usr/bin/env Rscript
# vep_filter_1000gene.R
# Args: [1] input VEP-annotated VCF, [2] output directory, [3] gene panel XLSX

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript vep_filter_1000gene.R <input.vcf> <outdir> <gene_panel.xlsx>")
}

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(readxl)

input_file  <- args[1]
out_dir     <- args[2]
gene_panel  <- args[3]

# Parse VEP-annotated VCF (skip meta lines, read from #CHROM header)
df <- fread(file = input_file, sep = '\t', header = TRUE, skip = '#CHROM')

new_df <- df %>%
    separate(INFO, c("END", "CSQ"), sep = ";") %>%
    mutate(CSQ = gsub("CSQ=", "", CSQ)) %>%
    separate_rows(CSQ, sep = ",") %>%
    separate(
        col   = CSQ,
        into  = c(
            "Allele", "Consequence", "IMPACT", "SYMBOL",
            "Gene", "Feature_type", "Feature", "Biotype",
            "Exon", "intron", "HGVSc", "HGVSp", "cDNA_pos",
            "CDS_position", "Protein_position", "Amino_acids",
            "Codons", "Existing_variation", "Distance",
            "STRAND", "Flags", "Symbol_source", "HGNC_ID",
            "MANE_SELECT", "SV_overlap_AF", "SV_overlap_PC",
            "SV_overlap_name"
        ),
        sep   = "\\|",
        extra = "drop",
        fill  = "right"
    )

# Load gene panel
panel_df  <- read_excel(gene_panel)
panel_vec <- panel_df[["Gene Symbol"]]

# Filter: gene in panel AND MANE_SELECT transcript present
filter_df <- new_df %>%
    filter(SYMBOL %in% panel_vec) %>%
    filter(!is.na(MANE_SELECT) & MANE_SELECT != "")

# Build output filename from input
base_name <- tools::file_path_sans_ext(basename(input_file))
out_file  <- file.path(out_dir, paste0(base_name, ".filtered.tsv"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_delim(x = filter_df, file = out_file, delim = "\t", col_names = TRUE)

message("Written: ", out_file)
message("Variants retained: ", nrow(filter_df))
