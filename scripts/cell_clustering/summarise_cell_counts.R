# File description --------------------------------------------------------

### Goal
# summarise_cell_counts.R produces a summary of cell counts by cell type.

### Execution steps
# x make_feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# x cell_clustering/build_adata.py
# x cell_clustering/integration.py
# x update_cell_metadata/add_cell_type_annotation.R
# ===>>> cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - dge_gsea/summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.rds produced by add_cell_type_annotation.R
# - paramstream package available from GitHub with 'devtools::install_github("benvallin/paramstream")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/cell_clustering/summarise_cell_counts"

# Path to directory containing cell metadata file
cell_metadata_dir_path <- "publications/Mock_2026/data/update_cell_metadata/add_cell_type_annotation"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(ggpubr)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define analysis parameters ----------------------------------------------

input_params = list(min_libsize_incl = 20000,
                    min_gene_incl = 3000,
                    max_pct_mt_incl = 20,
                    n_top_genes = 2000,
                    cat_cov_key = "line_name_generic",
                    integration_method = "scvi")

# Stream analysis parameter values ----------------------------------------

input_param_log_id <- get_param_log_id(param_list = input_params,
                                       in_dir_path = cell_metadata_dir_path)

write_param_log(param_list = input_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = input_params,
                                        in_dir_path = out_dir_path)

# Load cell metadata ------------------------------------------------------

# Load cell metadata
cell_metadata <- read_rds(file.path(cell_metadata_dir_path, 
                                    paste0("cell_metadata_", input_param_log_id, ".rds")))

# Subset cells passing barcode QC
cell_metadata <- cell_metadata %>% 
  filter(pass_qc == TRUE)

# Cell type summary -------------------------------------------------------

# Cell type summary - cell types by line ID 
cell_type_summary <- cell_metadata %>% 
  select(line_name_generic, cell_type) %>% 
  group_by(line_name_generic, cell_type) %>% 
  add_count() %>% 
  ungroup() %>% 
  unique() %>% 
  group_by(line_name_generic) %>% 
  mutate(n_all = sum(n)) %>% 
  ungroup() %>% 
  mutate(pct = round(n / n_all * 100, 2)) %>% 
  ungroup()

ggplot(data = cell_type_summary,
       mapping = aes(x = line_name_generic,
                     y = pct,
                     fill = fct_relevel(cell_type, rev))) +
  geom_col(color = "black") +
  geom_label(data = cell_type_summary %>% 
               filter(cell_type == "neuron"),
             mapping = aes(y = 0.9*pct,
                           label = paste0(pct, "%")),
             fill = "white", size = 3) +
  scale_fill_manual(limits = levels(cell_metadata$cell_type),
                    values = levels(cell_metadata$cell_type_color)) +
  scale_y_continuous(limits = c(0L, 101L), breaks = seq(0L, 100L, 10L)) +
  theme_pubr(legend = "bottom") +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0, vjust = 1),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12)) +
  labs(x = NULL,
       y = "Proportion (%)",
       fill = NULL) 

ggsave(filename = file.path(out_dir_path, paste0("cell_type_by_line.pdf")),
       width = 14,
       height = 6,
       device = cairo_pdf)




