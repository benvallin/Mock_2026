# File description --------------------------------------------------------

### Goal
# add_count_metrics.R updates cell metadata with Salmon's count metrics.

### Execution steps
# x make_feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# ===>>> update_cell_metadata/add_count_metrics.R
# - update_cell_metadata/mark_barcode_quality.R
# - cell_clustering/build_adata.py
# - cell_clustering/integration.py
# - update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - dge_gsea/summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.csv available on GitHub repository 'benvallin/Mock_2026', at metadata/cell_metadata.csv
# - gene_metadata.csv file produced by make_feature_metadata.R
# - Gene-level estimated counts produced by summarise_to_gene_level.R
# - Salmon's quant.sf and meta_info.json files available from the Gene Expression Omnibus (GEO) under accession number GSE354120
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/update_cell_metadata/add_count_metrics"

# Path to cell metadata file (csv)
cell_metadata_path <- "publications/Mock_2026/metadata/cell_metadata.csv"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Mock_2026/data/make_feature_metadata/v44/gene_metadata.csv"

# Path to directory containing counts files
cnt_dir_path <- "publications/Mock_2026/data/count_processing/summarise_to_gene_level"

# Path to directory containing Salmon quant.sf and meta_info.json files
# salmon_quant_dir_path <- "publications/Mock_2026/seq_data_Mock_2026"
salmon_quant_dir_path <- "seq_data/ben/pff_single_cell/gencode_v46/04_salmon_quant"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(data.table)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE) 

# Load metadata -----------------------------------------------------------

# Load cell metadata
cell_metadata <- read_csv(cell_metadata_path)

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Load counts -------------------------------------------------------------

# Load raw counts
cnt <- fread(file.path(cnt_dir_path, "cnt_lengthScaledTPM.csv"))

# Load TPMs
tpm <- fread(file.path(cnt_dir_path, "tpm_lengthScaledTPM.csv"))

# Update cell metadata with Salmon metrics --------------------------------

cell_metadata <- cell_metadata %>%
  mutate(salmon_data = map(.x = barcode_spe_dir,
                           .f = ~ import_salmon_meta_info(file = file.path(normalizePath(salmon_quant_dir_path), .x, "aux_info/meta_info.json")) %>%
                             select(num_processed, num_mapped, percent_mapped, frag_length_mean))) %>%
  unnest(salmon_data)

# Update cell metadata with count type metrics ----------------------------

# Compute library size using raw counts
libsize <- colSums(select(cnt, all_of(cell_metadata$barcode)))

libsize <- tibble(barcode = names(libsize),
                  libsize = unname(libsize))

# Compute count type summary using TPMs
tpm_summary <- summarise_gene_types(input = tpm,
                                    gene_id = ensembl_gene_id_version,
                                    sample_id = barcode,
                                    gene_metadata = gene_metadata) %>%
  select(-libsize)

# Update cell metadata with count type summary
cell_metadata <- cell_metadata %>%
  left_join(libsize,
            by = join_by(barcode)) %>%
  left_join(tpm_summary,
            by = join_by(barcode))

# Tidy updated cell metadata
cell_metadata <- cell_metadata %>%
  select(DIV:n_cells, num_processed:pct_rb_genes, barcode_spe_dir)

# Save updated cell metadata ----------------------------------------------

# Write updated cell metadata to disk
write_csv(cell_metadata, file.path(out_dir_path, "cell_metadata.csv"))
