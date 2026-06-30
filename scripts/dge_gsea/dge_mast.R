# File description --------------------------------------------------------

### Goal
# dge_mast.R performs differential gene expression analysis (DGEA) based on user-defined parameters.

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# x feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# x cell_clustering/build_adata.py
# x cell_clustering/integration.py
# x update_cell_metadata/add_cell_type_annotation.R
# x cell_clustering/summarise_cell_counts.R
# ===>>> dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.rds produced by add_cell_type_annotation.R
# - gene_metadata.csv file produced by make_feature_metadata.R
# - Gene-level estimated counts produced by summarise_to_gene_level.R
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
out_dir_path <- "publications/Mock_2026/data/output/dge_gsea/dge_mast"

# Path to directory containing cell metadata file
cell_metadata_dir_path <- "publications/Mock_2026/data/output/update_cell_metadata/add_cell_type_annotation"

# Path to gene metadata file (csv)
gene_metadata_path <- "publications/Mock_2026/data/output/make_feature_metadata/v44/gene_metadata.csv"

# Path to directory containing counts files
cnt_dir_path <- "publications/Mock_2026/data/output/count_processing/summarise_to_gene_level"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(MAST)
library(rnaseqtools)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Set n cores for parallel processing
options(mc.cores = 48)

# Set seed for reproducibility
set.seed(27)

# Define input analysis parameters ----------------------------------------

input_params = list(min_libsize_incl = 20000,
                    min_gene_incl = 3000,
                    max_pct_mt_incl = 20,
                    n_top_genes = 2000,
                    cat_cov_key = "line_name_generic",
                    integration_method = "scvi")

# Stream input analysis parameter values ----------------------------------

input_param_log_id <- get_param_log_id(param_list = input_params,
                                       in_dir_path = cell_metadata_dir_path)

# Load cell metadata ------------------------------------------------------

# Load cell metadata
cell_metadata <- read_rds(file.path(cell_metadata_dir_path, 
                                    paste0("cell_metadata_", input_param_log_id, ".rds")))

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

# Check DGE analysis parameters -------------------------------------------

# Check unique dge_params entries
n_dge_params <- length(names(dge_params))
n_unique_dge_params <- length(unique(names(dge_params)))

if(n_dge_params != n_unique_dge_params) {
  
  non_unique <- setNames(object = as.integer(table(names(dge_params))), nm = names(table(names(dge_params))))
  non_unique <- names(non_unique[non_unique > 1])
  
  stop("The following dge_params entries are not unique: ", paste0(non_unique, collapse = ", "))
  
}

# Define required dge_params entries
req_dge_params <- c("input_count", "coef_var", "coef_ref_lvl", 
                    "coef_test_lvl", "model_formula", "padj_dge")

# Define optional dge_params entries of type logical
opt_lgl_dge_params <- c("exclude_mixed_ctrls",
                        "exclude_doublets", 
                        "sc_lvl_gene_filter", 
                        "protein_coding_only",
                        "exclude_mt_genes")

# Define optional dge_params entries of type numeric
opt_num_dge_params <- c("max_prop_ambient_trip_incl",
                        "max_pct_mt_incl",
                        "min_cnt_excl", 
                        "min_freq_incl")

# Define optional dge_params entries of type character
opt_chr_dge_params <- c("singlet_filt_var")

# Define dge_params entries which should be names of logical columns in cell_metadata
met_lgl_dge_params <- c("singlet_filt_var")

# Define reserved dge_params entries
res_dge_params <- unique(c(req_dge_params, opt_lgl_dge_params, opt_num_dge_params, opt_chr_dge_params, met_lgl_dge_params))

# Define reserved entries (required + optional) included in user's dge_params
for(i in c("req_dge_params", "opt_lgl_dge_params", "opt_num_dge_params", 
           "opt_chr_dge_params", "met_lgl_dge_params", "res_dge_params")) {
  
  assign(x = paste0("incl_", i),
         value = names(dge_params)[names(dge_params) %in% eval(as.symbol(i))])
  
}

# Check mandatory dge_params entries
missing <- req_dge_params[!req_dge_params %in% incl_req_dge_params]

if(length(missing) > 0L) {
  
  stop("The following mandatory entries are missing from dge_params: ", paste0(missing, collapse = ", "))
  
}

# Check optional dge_params entries of type logical
not_lgl <- incl_opt_lgl_dge_params[map_lgl(.x = incl_opt_lgl_dge_params, 
                                           .f = ~ !is.logical(dge_params[[.x]]))]

if(length(not_lgl) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected logical type: ", paste0(not_lgl, collapse = ", "))
  
}

# Check optional dge_params entries of type numeric
not_num <- incl_opt_num_dge_params[map_lgl(.x = incl_opt_num_dge_params, 
                                           .f = ~ !is.numeric(dge_params[[.x]]))]

if(length(not_num) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected numeric type: ", paste0(not_num, collapse = ", "))
  
}

# Check optional dge_params entries of type character
not_chr <- incl_opt_chr_dge_params[map_lgl(.x = incl_opt_chr_dge_params, 
                                           .f = ~ !is.character(dge_params[[.x]]))]

if(length(not_chr) > 0L) {
  
  stop("The following optional dge_params entries are not of the expected character type: ", paste0(not_chr, collapse = ", "))
  
}

# Check dependencies for reserved dge_params entries 
if("exclude_doublets" %in% incl_res_dge_params &&
   dge_params$exclude_doublets &&
   !"singlet_filt_var" %in% incl_res_dge_params) {
  
  stop("The optional dge_params entry 'exclude_doublets' requires 'singlet_filt_var'")
  
}

if("sc_lvl_gene_filter" %in% incl_res_dge_params &&
   dge_params$sc_lvl_gene_filter &&
   !all(c("min_cnt_excl", "min_freq_incl") %in% incl_res_dge_params)) {
  
  stop("The optional dge_params entry 'sc_lvl_gene_filter' requires 'min_cnt_excl' and 'min_freq_incl'")
  
}

# Check dge_params entries which should be names of logical columns in cell_metadata
met_lgl_exp <- setNames(object = map_chr(.x = incl_met_lgl_dge_params, 
                                         .f = ~ as.character(dge_params[[.x]])), 
                        nm = incl_met_lgl_dge_params)

not_met_lgl <- met_lgl_exp[!map_lgl(.x = met_lgl_exp, 
                                    .f = ~ .x %in% colnames(cell_metadata) && is.logical(cell_metadata[[.x]]))]   

if(length(not_met_lgl) > 0L) {
  
  stop("The following dge_params entries are not names of logical columns in cell_metadata: ", 
       paste0(paste0(names(not_met_lgl), " - ", not_met_lgl), collapse = ", "))
  
}

# Define expected cell filter vars
# => dge_params entries that are not reserved
exp_cell_filt_vars <- setdiff(names(dge_params), incl_res_dge_params)

# Define actual cell filter vars
# => Expected cell filter vars which are column names of cell_metadata
cell_filt_vars <- exp_cell_filt_vars[exp_cell_filt_vars %in% colnames(cell_metadata)]

# Define problematic cell filter vars
# => Expected cell filter vars which are not column names of cell_metadata
prob_cell_filt_vars <- exp_cell_filt_vars[!exp_cell_filt_vars %in% cell_filt_vars]

if(length(prob_cell_filt_vars) > 0) {
  
  stop("The following dge_params entries are not column names of cell_metadata: ", paste0(prob_cell_filt_vars, collapse = ", "))
  
}

# Load gene metadata ------------------------------------------------------

# Load gene metadata
gene_metadata <- read_csv(gene_metadata_path)

# Filter cell metadata ----------------------------------------------------

# Subset to cells passing barcode QC 
cell_metadata <- cell_metadata %>% 
  filter(pass_qc == TRUE)

# Subset to cell filter var values defined in dge_params
if(length(cell_filt_vars) > 0) {
  
  for(i in cell_filt_vars) {
    
    cell_metadata <- cell_metadata[cell_metadata[[i]] %in% dge_params[[i]],]
    
  }
  
}

# Filter on max % mitochondrial TPMs
if("max_pct_mt_incl" %in% incl_res_dge_params) {
  
  cell_metadata <- cell_metadata %>%
    filter(pct_mt_genes <= dge_params$max_pct_mt_incl)
  
}

# Set coef_var levels
cell_metadata[[dge_params$coef_var]] = fct_relevel(cell_metadata[[dge_params$coef_var]],
                                                   c(dge_params$coef_ref_lvl, dge_params$coef_test_lvl))

# Drop unused levels in categorical variables
cell_metadata <- cell_metadata %>%
  mutate(across(.cols = where(is.factor),
                .fns = droplevels))

# Construct count matrix --------------------------------------------------

# Load input counts
cnt <- fread(file.path(cnt_dir_path, paste0(dge_params$input_count, ".csv"))) %>%
  column_to_rownames(var = "ensembl_gene_id_version") %>%
  as.matrix()

# Keep only the barcode present cell_metadata
cnt <- cnt[, cell_metadata$barcode]
identical(colnames(cnt), cell_metadata$barcode)

# Filter low freq genes ---------------------------------------------------

if(all(c("sc_lvl_gene_filter", "min_cnt_excl", "min_freq_incl") %in% incl_res_dge_params) && 
   dge_params$sc_lvl_gene_filter) {
  
  keep <- (rowSums(cnt > dge_params$min_cnt_excl) / dim(cnt)[[2]]) >= dge_params$min_freq_incl
  cnt <- cnt[keep,]
  
}

# Extra gene filtering ----------------------------------------------------

# Subset to protein-coding genes only
if("protein_coding_only" %in% incl_res_dge_params && 
   dge_params$protein_coding_only) {
  
  protein_coding_genes <- gene_metadata %>%
    filter(gene_type == "protein_coding") %>%
    pull(ensembl_gene_id_version)
  
  cnt <- cnt[rownames(cnt) %in% protein_coding_genes,]
  
}

# Exclude mitochondrial genes
if("exclude_mt_genes" %in% incl_res_dge_params && 
   dge_params$exclude_mt_genes) {
  
  mt_genes <- gene_metadata %>%
    filter(str_detect(gene_symbol, "^MT-")) %>%
    pull(ensembl_gene_id_version)
  
  cnt <- cnt[!rownames(cnt) %in% mt_genes,]
  
}

# MAST workflow -----------------------------------------------------------

# Produced log2 COUNT+1
log2_cnt <- log(x = cnt+1L, base = 2)

# Construct SingleCellAssay
sca <- FromMatrix(exprsArray = log2_cnt,
                  cData = as.data.frame(cell_metadata),
                  check_sanity = T)

# Compute and scale cellular detection rate (number of genes detected in a cell)
colData(sca)$n_gene_on <- scale(colSums(assay(sca) > 0))

# Define test
test <- ifelse(str_detect(dge_params$model_formula, "\\|"), "glmer", "bayesglm")

# Fit hurdle model
if(test == "glmer") {
  
  zlm <- zlm(formula = formula(dge_params$model_formula),
             sca = sca,
             method = "glmer",
             parallel = TRUE,
             ebayes = FALSE,
             strictConvergence = FALSE,
             fitArgsD = list(nAGQ = 1,
                             control = lme4::glmerControl(optimizer = "bobyqa",
                                                          optCtrl = list(maxfun = 1e7),
                                                          tolPwrss = 1e-7,
                                                          check.conv.grad = "warning",
                                                          check.conv.singular = "warning",
                                                          calc.derivs = TRUE)),
             fitArgsC = list(control = lme4::lmerControl(optimizer = "bobyqa",
                                                         optCtrl = list(maxfun = 1e7))))
} else if(test == "bayesglm") {
  
  zlm <- zlm(formula = formula(dge_params$model_formula),
             sca = sca,
             method = "bayesglm",
             ebayes = TRUE,
             parallel = TRUE)
  
}

# Perform likelihood ratio test
lrt <- summary(object = zlm,
               logFC = TRUE,
               doLRT = paste0(dge_params$coef_var, dge_params$coef_test_lvl),
               level = 1-dge_params$padj_dge,
               parallel = TRUE)

# Extract results
lfc <- tidy_summaryZlmFit(input = lrt,
                          coef = paste0(dge_params$coef_var, dge_params$coef_test_lvl),
                          cnt_scale = "log2",
                          level = 1-dge_params$padj_dge,
                          gene_id = "ensembl_gene_id_version")

lfc_fail <- lfc[is.na(lfc$log2FoldChange),]

if(test == "glmer" && nrow(lfc_fail) > 0L) {
  
  # For genes with missing estimates, run mixed effects model on continuous component only
  lfc_rescue <- fit_lmer_continuous(sca = sca,
                                    assay_name = "et",
                                    genes = lfc_fail$ensembl_gene_id_version,
                                    gene_id = "ensembl_gene_id_version",
                                    model_formula = dge_params$model_formula,
                                    coef_var = dge_params$coef_var,
                                    coef_test_lvl = dge_params$coef_test_lvl,
                                    cnt_scale = "log2",
                                    parallel = TRUE,
                                    control = lme4::lmerControl(optimizer = "bobyqa",
                                                                optCtrl = list(maxfun = 1e7)))
  
  lfc_rescue <- lfc_rescue[!is.na(lfc_rescue$log2FoldChange),]
  
  # Combine results from hurdle model and mixed effects model on continuous component only
  lfc <- lfc %>% 
    filter(!ensembl_gene_id_version %in% lfc_rescue$ensembl_gene_id_version) %>% 
    mutate(method = "zlm_glmer") %>% 
    bind_rows(lfc_rescue)
  
}

lfc <- lfc %>% 
  mutate(padj = stats::p.adjust(p = .data$pvalue, method = "BH"),
         sign_log2fc_times_minus_log10pvalue = sign(.data$log2FoldChange) * -log(x = .data$pvalue, base = 10)) %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")],
            by = join_by(ensembl_gene_id_version))

# Isolate DEGs
lfc_sign <- lfc[lfc$padj < dge_params$padj_dge & !is.na(lfc$pvalue),]

# Stream output analysis parameter values ---------------------------------

output_params <- c(input_params, dge_params)

write_param_log(param_list = output_params,
                out_dir_path = out_dir_path,
                param_log_table = TRUE)

output_param_log_id <- get_param_log_id(param_list = output_params,
                                        in_dir_path = out_dir_path)

# Write results to files --------------------------------------------------

write_rds(x = as_tibble(colData(sca), rownames = NA),
          file = file.path(out_dir_path, paste0("cell_metadata_", output_param_log_id, ".rds")))

write_csv(x = lfc,
          file = file.path(out_dir_path, paste0("lfc_", output_param_log_id, ".csv")))

write_csv(x = lfc_sign,
          file = file.path(out_dir_path, paste0("lfc_sign_", output_param_log_id, ".csv")))

fwrite(x = as_tibble(cnt, rownames = "ensembl_gene_id_version"),
       file = file.path(out_dir_path, paste0("cnt_", output_param_log_id, ".csv")))

fwrite(x = as_tibble(log2_cnt, rownames = "ensembl_gene_id_version"),
       file = file.path(out_dir_path, paste0("log2_cnt_1p_", output_param_log_id, ".csv")))

