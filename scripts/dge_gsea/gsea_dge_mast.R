# File description --------------------------------------------------------

### Goal
# gsea_dge_mast.R performs gene set enrichment analyses.

### Execution steps
# x make_feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# x cell_clustering/build_adata.py
# x cell_clustering/integration.py
# x update_cell_metadata/add_cell_type_annotation.R
# x cell_clustering/summarise_cell_counts.R
# x dge_gsea/dge_mast.R
# ===>>> dge_gsea/gsea_dge_mast.R
# - dge_gsea/summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - gene_metadata_msigdb.rds file produced by make_feature_metadata.R
# - DGEA results produced by dge_mast.R
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
out_dir_path <- "publications/Mock_2026/data/dge_gsea/gsea_dge_mast"

# Path to directory containing MAST DGE results
dge_dir_path <- "publications/Mock_2026/data/dge_gsea/dge_mast"

# Path to gene metadata file (rds)
gene_metadata_path <- "publications/Mock_2026/data/make_feature_metadata/v44/gene_metadata_msigdb.rds"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(rnaseqtools)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(27)

# Define input analysis parameters ----------------------------------------

input_params = list(min_libsize_incl = 20000,
                    min_gene_incl = 3000,
                    max_pct_mt_incl = 20,
                    n_top_genes = 2000,
                    cat_cov_key = "line_name_generic",
                    integration_method = "scvi")

# Define DGE analysis parameters ------------------------------------------

# PFF experiment: TRIP vs CTRL DIV61 neurons
dge_params = list(
  input_count = "tpm_lengthScaledTPM",
  coef_var = "genotype",
  coef_ref_lvl = "CTRL",
  coef_test_lvl = "TRIP",
  model_formula = "~ n_gene_on + genotype + (1 | line_name)",
  # model_formula = "~ n_gene_on + donor_sex + genotype + (1 | line_name)",
  padj_dge = 0.05,
  sc_lvl_gene_filter = TRUE,
  min_cnt_excl = 0,
  min_freq_incl = 0.2,
  protein_coding_only = TRUE,
  cell_type = "neuron",
  treatment = "untreated"
)

# Define GSEA parameters --------------------------------------------------

gsea_params <- list(
  min_size = 15,
  max_size = 500, 
  padj_gsea = 0.05,
  collections_names = c("MSigDB_H",
                        "MSigDB_C2_CP:KEGG_LEGACY",
                        "MSigDB_C2_CP:WIKIPATHWAYS",
                        "MSigDB_C2_CP:REACTOME",
                        "MSigDB_C5_GO:BP",
                        "MSigDB_C5_GO:CC")
)

# Stream analysis parameter values ----------------------------------------

input_param_log_id <- get_param_log_id(param_list = c(input_params, dge_params),
                                       in_dir_path = dge_dir_path)

output_params <- c(input_params, dge_params, gsea_params)

write_param_log(param_list = output_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = output_params,
                                        in_dir_path = out_dir_path)

# Load metadata -----------------------------------------------------------

# Load gene metadata
gene_metadata <- read_rds(gene_metadata_path)

# Load DGE results --------------------------------------------------------

cell_metadata <- read_rds(file = file.path(dge_dir_path, paste0("cell_metadata_", input_param_log_id, ".rds")))

lfc <- read_csv(file = file.path(dge_dir_path, paste0("lfc_", input_param_log_id, ".csv")))

lfc_sign <- read_csv(file = file.path(dge_dir_path, paste0("lfc_sign_", input_param_log_id, ".csv")))

# Perform GSEA ------------------------------------------------------------

# fgsea::fgsea()
fgsea <- multi_fgsea(input = gene_metadata,
                     collections = gsea_params$collections_names,
                     stats = lfc %>%
                       filter(is.finite(sign_log2fc_times_minus_log10pvalue)) %>% 
                       arrange(desc(sign_log2fc_times_minus_log10pvalue)) %>%
                       pull(sign_log2fc_times_minus_log10pvalue, ensembl_gene_id_version) %>%
                       na.omit(),
                     gene_id = "ensembl_gene_id_version",
                     min_set_size = gsea_params$min_size,
                     max_set_size = gsea_params$max_size,
                     padj_threshold = Inf) %>%
  arrange(padj)

fgsea_sign_up <- fgsea %>%
  filter(NES > 0,
         padj < gsea_params$padj_gsea)

fgsea_sign_down <- fgsea %>%
  filter(NES < 0,
         padj < gsea_params$padj_gsea)

# fgsea::fora()
fora_up <- multi_fora(input = gene_metadata,
                      collections = gsea_params$collections_names,
                      genes = lfc %>%
                        filter(log2FoldChange > 0, padj < dge_params$padj_dge) %>%
                        pull(ensembl_gene_id_version),
                      universe = lfc %>%
                        filter(!is.na(log2FoldChange)) %>%
                        pull(ensembl_gene_id_version),
                      gene_id = "ensembl_gene_id_version",
                      min_set_size = gsea_params$min_size,
                      max_set_size = gsea_params$max_size,
                      padj_threshold = Inf) %>%
  arrange(padj)

fora_sign_up <- fora_up %>%
  filter(padj < gsea_params$padj_gsea) %>%
  arrange(padj)

fora_down <- multi_fora(input = gene_metadata,
                        collections = gsea_params$collections_names,
                        genes = lfc %>%
                          filter(log2FoldChange < 0, padj < dge_params$padj_dge) %>%
                          pull(ensembl_gene_id_version),
                        universe = lfc %>%
                          filter(!is.na(log2FoldChange)) %>%
                          pull(ensembl_gene_id_version),
                        gene_id = "ensembl_gene_id_version",
                        min_set_size = gsea_params$min_size,
                        max_set_size = gsea_params$max_size,
                        padj_threshold = Inf) %>%
  arrange(padj)

fora_sign_down <- fora_down %>%
  filter(padj < gsea_params$padj_gsea) %>%
  arrange(padj)

# Write results to file ---------------------------------------------------

walk(.x = c("fgsea", "fgsea_sign_up", "fgsea_sign_down",
            "fora_up", "fora_down", "fora_sign_up", "fora_sign_down"),
     .f = ~ write_rds(x = eval(as.symbol(.x)),
                      file = file.path(out_dir_path, paste0(.x, "_", output_param_log_id, ".rds"))))

# -------------------------------------------------------------------------


"HALLMARK_OXIDATIVE_PHOSPHORYLATION" %in% fgsea_sign_down$set_name
"HALLMARK_FATTY_ACID_METABOLISM" %in% fgsea_sign_down$set_name
"REACTOME_AUTOPHAGY" %in% fgsea_sign_down$set_name
"KEGG_PROTEASOME" %in% fgsea_sign_down$set_name
"WP_PROTEASOME_DEGRADATION" %in% fgsea_sign_down$set_name
"REACTOME_PROTEASOME_ASSEMBLY" %in% fgsea_sign_down$set_name
"GOBP_LIPID_CATABOLIC_PROCESS" %in% fgsea_sign_down$set_name
