# File description --------------------------------------------------------

### Goal
# mark_barcode_quality.R updates cell metadata with barcode quality scoring based on user-defined thresholds.

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# x feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# ===>>> update_cell_metadata/mark_barcode_quality.R
# - cell_clustering/build_adata.py
# - cell_clustering/integration.py
# - update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.csv file produced by add_count_metrics.R
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
out_dir_path <- "publications/Mock_2026/data/output/update_cell_metadata/mark_barcode_quality"

# Path to directory containing cell metadata file
cell_metadata_path <- "publications/Mock_2026/data/output/update_cell_metadata/add_count_metrics/cell_metadata.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(paramstream)
library(ggpubr)
library(arrow)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Define analysis parameters ----------------------------------------------

# Define barcode QC parameters 
output_params <- list(min_libsize_incl = 20000,
                      min_gene_incl = 3000,
                      max_pct_mt_incl = 20)

# Stream analysis parameter values ----------------------------------------

write_param_log(param_list = output_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = output_params,
                                        in_dir_path = out_dir_path)

# Load cell metadata ------------------------------------------------------

# Load cell metadata
cell_metadata <- read_csv(cell_metadata_path)

# Mark barcode quality ----------------------------------------------------

cell_metadata <- cell_metadata %>% 
  mutate(pass_libsize = ifelse(libsize >= output_params$min_libsize_incl, TRUE, FALSE),
         pass_n_genes = ifelse(n_genes >= output_params$min_gene_incl, TRUE, FALSE),
         pass_pct_mt_genes = ifelse(pct_mt_genes <= output_params$max_pct_mt_incl, TRUE, FALSE),
         pass_qc = ifelse(pass_libsize & pass_n_genes & pass_pct_mt_genes, TRUE, FALSE))

# Set key variables as factors and provide colors -------------------------

# icell8_chip_id 
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "icell8_chip_id",
                                 levels = c("A", "B"),
                                 colors = c("#ffbe0c","#f0396a"))

# line_name
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "line_name",
                                 levels = c("B856", "B156", 
                                            "B067", "831-1",
                                            "831-3", "831-5"),
                                 colors = c("#3966ce", "#ff6d00", 
                                            "#0b8c1f", "#642ca9", 
                                            "#ffbe0c", "#f0396a"))

# line_name_generic
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "line_name_generic",
                                 levels = c("CTRL1", "CTRL2", 
                                            "CTRL3", "TRIP1",
                                            "TRIP2", "TRIP3"),
                                 colors = c("#3966ce", "#ff6d00", 
                                            "#0b8c1f", "#f0396a", 
                                            "#642ca9", "#ffbe0c"))

# genotype
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "genotype",
                                 levels = c("CTRL", "TRIP"),
                                 colors = c("#0077b6", "#a4133c"))

# donor_sex 
cell_metadata <- set_level_color(input = cell_metadata,
                                 variable = "donor_sex",
                                 levels = c("F", "M"),
                                 colors = c("#5d76cb", "#ca3767"))

# Save marked cell metadata -----------------------------------------------

# Write marked cell metadata to disk
write_parquet(cell_metadata, 
              file.path(out_dir_path, 
                        paste0("cell_metadata_", output_param_log_id, ".parquet")))

# QC plots ----------------------------------------------------------------

# Metrics plots
metrics_plots <- map2(.x = c("percent_mapped", "num_mapped", "num_processed", "frag_length_mean",
                             "libsize", "n_genes", "n_pc_genes", "n_mt_genes",  "n_rb_genes",
                             "pct_pc_genes", "pct_mt_genes", "pct_rb_genes"),
                      .y = c("mapped fragments (%)", "mapped fragments (million)", "processed fragments (million)", "mean fragment length (nt)",
                             "raw counts (million)", "genes (n)", "protein coding genes (n)", "mitochondrial RNAs (n)", "ribosomal RNAs (n)",
                             "protein coding genes (%)",  "mitochondrial RNAs (%)", "ribosomal RNAs (%)"),
                      
                      .f = ~ ggplot(data = cell_metadata %>% 
                                      mutate(across(.cols = c(num_mapped, num_processed, libsize),
                                                    .fns = ~ .x / 10^6)),
                                    mapping = aes(x = eval(sym(.x)))) +
                        geom_density(linewidth = 1, alpha = 0.3, fill = "#3966ce") +
                        theme_pubr(legend = "none") +
                        theme(axis.title = element_text(size = 16),
                              axis.text.x = element_text(size = 14),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank()) +
                        labs(x = .y)) %>%
  setNames(nm = c("percent_mapped", "num_mapped", "num_processed", "frag_length_mean",
                  "libsize", "n_genes", "n_pc_genes", "n_mt_genes",  "n_rb_genes",
                  "pct_pc_genes", "pct_mt_genes", "pct_rb_genes"))

p_pct_mapped <- metrics_plots$percent_mapped +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 10)) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "inside",
        legend.justification = c(0.1, 1),
        legend.background = element_blank())

p_pct_mapped

p_libsize <- metrics_plots$libsize +
  scale_x_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, 0.5)) +
  rotate_x_text(angle = 45) +
  geom_vline(xintercept = output_params$min_libsize_incl / 10^6,
             linetype = "dashed",
             color = "orange",
             linewidth = 1) 

p_libsize

p_n_genes <- metrics_plots$n_genes +
  scale_x_continuous(limits = c(0, 16000),
                     breaks = seq(0, 16000, 2000)) +
  geom_vline(xintercept = output_params$min_gene_incl,
             linetype = "dashed",
             color = "orange",
             linewidth = 1)

p_n_genes

p_pct_mt_genes <- metrics_plots$pct_mt_genes +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 10)) +
  geom_vline(xintercept = output_params$max_pct_mt_incl,
             linetype = "dashed",
             color = "orange",
             linewidth = 1)

p_pct_mt_genes

p_n_pc_genes <- metrics_plots$n_pc_genes +
  scale_x_continuous(limits = c(0, 16000),
                     breaks = seq(0, 16000, 2000))

p_n_pc_genes

p_pct_pc_genes <- metrics_plots$pct_pc_genes +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 10))

p_pct_pc_genes

p_pct_rb_genes <- metrics_plots$pct_rb_genes +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 10))

p_pct_rb_genes

# Display cell numbers before / after filtering
summary <- cell_metadata %>%
  select(barcode, pass_qc) %>%
  nest(data = c(barcode, pass_qc)) %>%
  mutate(n_cells = map_int(.x = data, .f = ~ dim(.x)[[1]]),
         n_pass = map_int(.x = data, .f = ~ dim(.x[.x$pass_qc == TRUE,])[[1]]),
         pct_pass = n_pass / n_cells *100)

p_n_cells <- ggplot(data = summary,
                    mapping = aes(x = 1)) +
  geom_col(mapping = aes(y = n_cells),
           fill = "grey",
           color = "black", alpha = 0.8) +
  geom_col(mapping = aes(y = n_pass),
           color = "black", fill = "#3966ce", alpha = 0.8) +
  geom_label(mapping = aes(y = n_pass - 200,
                           label = paste0(n_pass, "\n(", round(pct_pass, 1), "%)")),
             fill = "white", size = 5) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, 400)) +
  theme_pubr(legend = "none") +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = "inside",
        legend.justification = c(0.5, 1),
        legend.background = element_blank()) +
  labs(x = NULL,
       y = "barcodes passing QC (n)")

p_n_cells

# Write plots to disk -----------------------------------------------------

walk(.x = syms(objects(pattern = "^p_.*$")),
     .f = ~ ggsave(plot = eval(.x),
                   filename = file.path(out_dir_path, paste0(.x, ".pdf")),
                   width = 7,
                   height = 7,
                   device = cairo_pdf))
