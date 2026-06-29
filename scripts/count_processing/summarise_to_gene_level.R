# File description --------------------------------------------------------

### Goal
# summarise_to_gene_level.R produces gene count estimates from Salmon's transcript abundance estimates

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# x feature_metadata/make_feature_metadata.R
# ===>>> count_processing/summarise_to_gene_level.R
# - update_cell_metadata/add_count_metrics.R
# - update_cell_metadata/mark_barcode_quality.R
# - cell_clustering/build_adata.py
# - cell_clustering/integration.py
# - update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.csv file produced by make_cell_metadata.R
# - transcript_metadata.csv file produced by make_feature_metadata.R
# - Salmon's quant.sf files 

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/output/count_processing/summarise_to_gene_level"

# Path to cell metadata file (csv)
cell_metadata_path <- "publications/Mock_2026/data/output/make_cell_metadata/cell_metadata.csv"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "publications/Mock_2026/data/output/make_feature_metadata/v44/transcript_metadata.csv"

# Path to directory containing Salmon quant files
salmon_quant_dir_path <- "seq_data/ben/pff_single_cell/gencode_v46/04_salmon_quant"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(tximport)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE) 

# Import metadata --------------------------------------------------

# Import cell metadata
cell_metadata <- read_csv(cell_metadata_path) %>%
  mutate(files = setNames(object = file.path(normalizePath(salmon_quant_dir_path), barcode_spe_dir, "quant.sf"),
                          nm = barcode)) %>%
  select(files, everything())

# Import transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path)

# Summarise transcript data to gene level ---------------------------------

# Construct txi object
# - Using countsFromAbundance = "no"
# - Using countsFromAbundance = "scaledTPM"
# => Generate estimated counts using abundance estimates scaled up to library size
# - Using countsFromAbundance = "lengthScaledTPM"
# => Generate estimated counts using abundance estimates scaled using the average transcript length over samples and then the library size

for(i in c("no", "scaledTPM", "lengthScaledTPM")) {
  
  message("Processing countsFromAbundance: ", i, "\n")
  
  txi_nm <- paste0("txi", ifelse(i == "no", "", paste0("_", i)))
  
  assign(x = txi_nm,
         value = tximport(files = cell_metadata$files,
                          type = "salmon",
                          txOut = FALSE,
                          countsFromAbundance = i,
                          tx2gene = transcript_metadata))
  
  write_rds(x = eval(as.symbol(txi_nm)),
            file = file.path(out_dir_path, paste0(txi_nm, ".rds")))
  
}

# Load txi(s) object(s) ---------------------------------------------------

for(i in c("txi", "txi_scaledTPM", "txi_lengthScaledTPM")) {
  
  assign(x = i,
         value = read_rds(file.path(out_dir_path, paste0(i, ".rds"))))
  
}

# Extract estimated raw read counts and TPMs ------------------------------

for(i in c("counts", "abundance")) {
  
  for(j in c("txi", "txi_scaledTPM", "txi_lengthScaledTPM")) {
    
    nm_prefix <- switch(i,
                        counts = "cnt",
                        abundance = "tpm",
                        stop("Invalid count type value."))
    
    nm <- paste0(nm_prefix,
                 regmatches(x = j, m = regexpr(pattern = "_.*$", text = j)))
    
    message("Processing: ", nm, "\n")
    
    assign(x = nm,
           value = as_tibble(eval(as.symbol(j))[[i]],
                             rownames = "ensembl_gene_id_version"))
    
    fwrite(x = eval(as.symbol(nm)), file = file.path(out_dir_path, paste0(nm, ".csv")))
    
  }
  
}
