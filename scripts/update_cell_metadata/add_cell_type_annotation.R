# File description --------------------------------------------------------

### Goal
# add_cell_type_annotation.R updates cell metadata with cell type annotation results.

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# x feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# x cell_clustering/build_adata.py
# x cell_clustering/integration.py
# ===>>> update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.parquet file produced by mark_barcode_quality.R
# - cell_type_annotation.csv file produced by integration.py
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'
# - paramstream package available from GitHub with 'devtools::install_github("benvallin/paramstream")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/output/update_cell_metadata/add_cell_type_annotation"

# Path to directory containing cell metadata file
cell_metadata_dir_path <- "publications/Mock_2026/data/output/update_cell_metadata/mark_barcode_quality"

# Path to directory containing cell clustering results
clustering_dir_path <- "publications/Mock_2026/data/output/cell_clustering/integration"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(paramstream)
library(arrow)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE) 

# Define analysis parameters ----------------------------------------------

input1_params = list(min_libsize_incl = 20000,
                     min_gene_incl = 3000,
                     max_pct_mt_incl = 20)

input2_params <- c(input1_params,
                   list(n_top_genes = 2000,
                        cat_cov_key = "line_name_generic",
                        integration_method = "scvi"))

# Stream analysis parameter values ----------------------------------------

input1_param_log_id <- get_param_log_id(param_list = input1_params,
                                        in_dir_path = cell_metadata_dir_path)

input2_param_log_id <- get_param_log_id(param_list = input2_params,
                                        in_dir_path = clustering_dir_path)

write_param_log(param_list = input2_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = input2_params,
                                        in_dir_path = out_dir_path)

# Load metadata -----------------------------------------------------------

# Load cell metadata
cell_metadata <- read_parquet(file.path(cell_metadata_dir_path, 
                                        paste0("cell_metadata_", input1_param_log_id, ".parquet")))

# Load cell type annotation
cell_type_annotation <- read_csv(file.path(clustering_dir_path, 
                                           paste0("cell_type_annotation_", input2_param_log_id, ".csv"))) 

# Add cell type annotation to cell metadata -------------------------------

cell_metadata <- cell_metadata %>% 
  left_join(cell_type_annotation,
            by = join_by(barcode))

# Set key variables as factors and provide colors -------------------------

# cell_type
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "cell_type",
                                 levels = c("neuron", "astrocyte", "fibroblast-like", NA_character_),
                                 colors = c("#5d76cb", "#ca3767", "#29a655", "#b7b7b7"))

# Write updated cell metadata to disk
write_rds(cell_metadata, file.path(out_dir_path, paste0("cell_metadata_", output_param_log_id, ".rds")))

