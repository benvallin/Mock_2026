# File description --------------------------------------------------------

### Goal
# make_feature_metadata.R produces gene and transcript metadata from GENECODE annotation.

### Execution steps
# ===>>> make_feature_metadata/make_feature_metadata.R
# - count_processing/summarise_to_gene_level.R
# - update_cell_metadata/add_count_metrics.R
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
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/make_feature_metadata"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(arrow)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define annotation parameters --------------------------------------------

species <- "human"

release <- 44

# Update paths ------------------------------------------------------------

sr_dir <- ifelse(species == "mouse",
                 paste0("vM", release),
                 paste0("v", release))

out_dir_sub_path <- file.path(out_dir_path, sr_dir)

# Download GENCODE annotation ---------------------------------------------

download_gencode_annotation(species = species,
                            release = release,
                            out_dir_path = out_dir_sub_path)

# Import GENCODE annotation -----------------------------------------------

gencode_annotation <- import_gencode_annotation(file = list.files(path = out_dir_sub_path,
                                                                  pattern = "^.*\\.gtf\\.gz$",
                                                                  full.names = TRUE))

# Add transgene annotation ------------------------------------------------

transgenes <- c("mTagBFP2", "puromycinR")

transgene_annotation <- tibble(gene_id = transgenes,
                               transcript_id = transgenes,
                               gene_name = transgenes,
                               gene_type = "protein_coding")

gencode_annotation <- bind_rows(gencode_annotation, transgene_annotation)

# Construct gene metadata from GENCODE annotation -------------------------

gene_metadata <- make_gene_metadata(input = gencode_annotation)

write_csv(gene_metadata, file.path(out_dir_sub_path, "gene_metadata.csv"))

# Add Molecular Signatures Database (MSigDB) annotation -------------------

msigdb_species <- ifelse(species == "human", "Homo sapiens", "Mus musculus")

msigdb_db_species <- ifelse(species == "human", "Hs", "Mm")

msigdb_collection_table = get_msigdb_collections(species = msigdb_species,
                                                 db_species = msigdb_db_species)

gene_metadata <- gene_metadata %>%
  left_join(msigdb_collection_table, by = join_by(ensembl_gene_id))

write_rds(gene_metadata, file.path(out_dir_sub_path, "gene_metadata_msigdb.rds"))

write_parquet(gene_metadata, file.path(out_dir_sub_path, "gene_metadata_msigdb.parquet"))

# Construct transcript metadata from GENCODE annotation -------------------

transcript_metadata <- make_transcript_metadata(input = gencode_annotation)

write_csv(transcript_metadata, file.path(out_dir_sub_path, "transcript_metadata.csv"))
