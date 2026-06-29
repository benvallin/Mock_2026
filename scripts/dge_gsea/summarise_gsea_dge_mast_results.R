# File description --------------------------------------------------------

### Goal
# summarise_gsea_dge_mast_results.R summarises DGEA and gene set enrichment analysis results.

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
# x dge_gsea/dge_mast.R
# x dge_gsea/gsea_dge_mast.R
# ===>>> summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - gene_metadata_msigdb.rds file produced by make_feature_metadata.R
# - DGEA results produced by dge_mast.R
# - Gene set enrichment analysis results produced by gsea_dge_mast.R
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
out_dir_path <- "publications/Mock_2026/data/output/dge_gsea/summarise_gsea_dge_mast_results"

# Path to directory containing MAST DGE results
dge_dir_path <- "publications/Mock_2026/data/output/dge_gsea/dge_mast"

# Path to directory containing GSEA results (MAST)
gsea_dir_path <- "publications/Mock_2026/data/output/dge_gsea/gsea_dge_mast"

# Path to gene metadata file (rds)
gene_metadata_path <- "publications/Mock_2026/data/output/make_feature_metadata/v44/gene_metadata_msigdb.rds"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(ggpubr)
library(kableExtra)
library(rnaseqtools)
library(paramstream)

# Set working directory
setwd(work_dir_path)

# Create output directory if necessary
dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

# Stream analysis parameter values ----------------------------------------

input_params = list(min_libsize_incl = 20000,
                    min_gene_incl = 3000,
                    max_pct_mt_incl = 20,
                    n_top_genes = 2000,
                    cat_cov_key = "line_name_generic",
                    integration_method = "scvi")

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

gsea_params <- list(
  min_size = 15,
  max_size = 500, 
  padj_gsea = 0.05,
  collections_names = c("MSigDB_H",
                        "MSigDB_C2_CP:KEGG_LEGACY",
                        "MSigDB_C2_CP:WIKIPATHWAYS",
                        "MSigDB_C2_CP:REACTOME",
                        "MSigDB_C2_CP:BIOCARTA",
                        "MSigDB_C5_GO:BP")
)

dge_log_id <- get_param_log_id(param_list = c(input_params, dge_params),
                               in_dir_path = dge_dir_path)

gsea_log_id <- get_param_log_id(param_list = c(input_params, dge_params, gsea_params),
                                in_dir_path = gsea_dir_path)

# Load DGE/GSEA results ---------------------------------------------------

cell_metadata <- read_rds(file.path(dge_dir_path, paste0("cell_metadata_", dge_log_id, ".rds")))

lfc <- read_csv(file.path(dge_dir_path, paste0("lfc_", dge_log_id, ".csv")))

lfc_sign <- read_csv(file.path(dge_dir_path, paste0("lfc_sign_", dge_log_id, ".csv")))

tpm <- fread(file.path(dge_dir_path, paste0("cnt_", dge_log_id, ".csv")))

fgsea_sign_up <- read_rds(file.path(gsea_dir_path, paste0("fgsea_sign_up_", gsea_log_id, ".rds")))

fgsea_sign_down <- read_rds(file.path(gsea_dir_path, paste0("fgsea_sign_down_", gsea_log_id, ".rds")))

fora_sign_up <- read_rds(file.path(gsea_dir_path, paste0("fora_sign_up_", gsea_log_id, ".rds")))

fora_sign_down <- read_rds(file.path(gsea_dir_path, paste0("fora_sign_down_", gsea_log_id, ".rds")))

# Load metadata -----------------------------------------------------------

# Load gene metadata
gene_metadata <- read_rds(gene_metadata_path)

# Abundance / DGEA summary ------------------------------------------------

lfc <- lfc %>%
  left_join(gene_metadata %>% 
              select(ensembl_gene_id_version, starts_with("MSigDB")), 
            by = join_by(ensembl_gene_id_version))

lfc_summary <- lfc %>% 
  select(ensembl_gene_id, gene_symbol, 
         log2FoldChange, lfcSE, pvalue, padj,
         MSigDB_H, `MSigDB_C2_CP:REACTOME`, `MSigDB_C2_CP:KEGG_LEGACY`, `MSigDB_C5_GO:BP`, `MSigDB_C2_CP:WIKIPATHWAYS`) %>% 
  mutate(across(.cols = c(log2FoldChange, lfcSE),
                .fn = ~ ifelse(is.nan(.x), NA, .x)),
         across(.cols = starts_with("MSigDB"),
                .fn = ~ map_chr(.x, 
                                .f = function(x) paste0(x$set_name, collapse = ", "))),
         trip_vs_ctrl = case_when(log2FoldChange > 0 ~ "up",
                                  log2FoldChange < 0 ~ "down",
                                  log2FoldChange == 0 ~ "unchanged",
                                  T ~ NA_character_),
         adj_pval_below_0.05 = ifelse(padj < 0.05, T, F)) %>% 
  arrange(desc(log2FoldChange)) %>% 
  select(`Ensembl gene ID` = ensembl_gene_id,
         `gene symbol` = gene_symbol,
         log2FoldChange,
         lfcSE,
         pval = pvalue,
         `adj. pval` = padj, 
         `TRIP vs CTRL` = trip_vs_ctrl, 
         `adj. pval < 0.05` = adj_pval_below_0.05, 
         `MSigDB - Hallmark` = MSigDB_H,
         `MSigDB - Reactome` = `MSigDB_C2_CP:REACTOME`,
         `MSigDB - KEGG` = `MSigDB_C2_CP:KEGG_LEGACY`,
         `MSigDB - GO BP` = `MSigDB_C5_GO:BP`,
         `MSigDB - Wikipathways` = `MSigDB_C2_CP:WIKIPATHWAYS`) 

tpm_summary_line_lvl <- full_join(compute_mean_tpm(input = tpm,
                                                   gene_id = ensembl_gene_id_version,
                                                   sample_metadata = cell_metadata,
                                                   sample_id_var = barcode,
                                                   group_id_vars = line_name), 
                                  compute_prop_expr(input = tpm,
                                                    gene_id = ensembl_gene_id_version,
                                                    sample_metadata = cell_metadata,
                                                    sample_id_var = barcode,
                                                    group_id_vars = line_name,
                                                    min_cnt = 0),
                                  by = join_by(line_name, ensembl_gene_id_version)) %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  select(`Ensembl gene ID` = ensembl_gene_id,
         `gene symbol` = gene_symbol,
         line_name,
         mean_tpm, sd_tpm, sem_tpm, mean_tpm_quartile,
         n, n_expressing, pct_expressing)

tpm_summary_genotype_lvl <- full_join(compute_mean_tpm(input = tpm,
                                                       gene_id = ensembl_gene_id_version,
                                                       sample_metadata = cell_metadata,
                                                       sample_id_var = barcode,
                                                       group_id_vars = genotype), 
                                      compute_prop_expr(input = tpm,
                                                        gene_id = ensembl_gene_id_version,
                                                        sample_metadata = cell_metadata,
                                                        sample_id_var = barcode,
                                                        group_id_vars = genotype,
                                                        min_cnt = 0),
                                      by = join_by(genotype, ensembl_gene_id_version)) %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  select(`Ensembl gene ID` = ensembl_gene_id,
         `gene symbol` = gene_symbol,
         genotype,
         mean_tpm, sd_tpm, sem_tpm, mean_tpm_quartile,
         n, n_expressing, pct_expressing)

walk(.x = c("lfc_summary", "tpm_summary_line_lvl", "tpm_summary_genotype_lvl"),
     .f = ~ write_csv(x = eval(sym(.x)),
                      file = file.path(out_dir_path, paste0(.x, ".csv"))))

# Draw lists of top DEGs --------------------------------------------------

# Draw lists of topN upregulated genes
top_deg_up <- lapply(X = seq(from = 5, to = 50, by = 5),
                     FUN = function(x) {
                       
                       lfc %>%
                         filter(log2FoldChange > 0) %>%
                         arrange(padj) %>%
                         head(x) %>%
                         arrange(desc(log2FoldChange)) %>%
                         pull(padj, ensembl_gene_id_version)
                       
                     }) %>%
  setNames(nm = paste0("top", seq(from = 5, to = 50, by = 5), "_deg_up"))

# Draw lists of topN downregulated genes
top_deg_down <- lapply(X = seq(from = 5, to = 50, by = 5),
                       FUN = function(x) {
                         
                         lfc %>%
                           filter(log2FoldChange < 0) %>%
                           arrange(padj) %>%
                           head(x) %>%
                           arrange(log2FoldChange) %>%
                           pull(padj, ensembl_gene_id_version)
                         
                       }) %>%
  setNames(nm = paste0("top", seq(from = 5, to = 50, by = 5), "_deg_down"))

# Overview of DEGs --------------------------------------------------------

# Plot top 50 upregulated genes
plot_lfc(input = lfc,
         genes = names(top_deg_up$top50_deg_up),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         title = "Top 50 upregulated genes",
         desc_lfc = TRUE,
         xintercepts = seq(0.5, 3, 0.5),
         lab_nudge = 0.1)

ggsave(filename = file.path(out_dir_path, "top50_deg_up.pdf"),
       width = 9,
       height = 2+0.22*50,
       device = cairo_pdf)

# Plot top 50 downregulated genes
plot_lfc(input = lfc,
         genes = names(top_deg_down$top50_deg_down),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         title = "Top 50 downregulated genes",
         desc_lfc = FALSE,
         xintercepts = c(seq(-9, -1, 1), 1),
         lab_nudge = 0.1)

ggsave(filename = file.path(out_dir_path, "top50_deg_down.pdf"),
       width = 9,
       height = 2+0.22*50,
       device = cairo_pdf)

# Violin plot of -log10(padj) by log2(FC) for non-missing padj genes color-coded by fold change - padj status
plot_df <- lfc %>%
  filter(!is.na(padj)) %>%
  mutate(minus_log10_padj = -log(x = padj, base = 10),
         fc_sign_status = case_when(log2FoldChange < 0 & padj < 0.05 ~ "FC < 0 - padj < 0.05",
                                    log2FoldChange > 0 & padj < 0.05 ~ "FC > 0 - padj < 0.05",
                                    padj >= 0.05 ~ paste0("padj ", "\u2265", " 0.05"),
                                    T ~ NA_character_) %>%
           fct_relevel(c(paste0("padj ", "\u2265", " 0.05"),
                         "FC < 0 - padj < 0.05",
                         "FC > 0 - padj < 0.05")),
         label = ifelse(ensembl_gene_id_version %in% c(names(top_deg_up$top30_deg_up), names(top_deg_down$top30_deg_down)),
                        gene_symbol,
                        NA_character_))

ggplot(data = plot_df,
       mapping = aes(x = log2FoldChange,
                     y = minus_log10_padj)) +
  geom_vline(xintercept = seq(-9, 6, 1),
             linetype = "solid", color = "#403d39", linewidth = 0.1) +
  geom_hline(yintercept = c(-log(x = 0.05, base = 10),
                            -log(x = 0.01, base = 10),
                            -log(x = 0.001, base = 10)),
             linetype = "dashed", color = "#ff9e00", linewidth = 0.3) +
  annotate(geom = "label",
           x = -9,
           y = c(-log(x = 0.05, base = 10),
                 -log(x = 0.01, base = 10),
                 -log(x = 0.001, base = 10)),
           label = c("padj = 0.05",
                     "padj = 0.01",
                     "padj = 0.001"),
           color = "black", size = 3) +
  geom_point(mapping = aes(fill = fc_sign_status),
             shape = 21, color = "black", size = 2, stroke = 0.2, alpha = 0.6) +
  ggrepel::geom_text_repel(data = plot_df %>%
                             filter(!is.na(label)),
                           mapping = aes(label = label),
                           max.overlaps = 500,
                           force = 2,
                           size = 3,
                           box.padding = 0.4,
                           segment.color = "#c0c0c0",
                           segment.size = 0.2) +
  scale_x_continuous(breaks = seq(-9, 6, 1)) +
  scale_y_continuous(n.breaks = 10L) +
  scale_fill_manual(values = c("#c0c0c0", "#ca3767", "#5d76cb")) +
  theme_pubr(legend = "bottom") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  labs(x = paste0("log", "\u2082", " fold change (TRIP vs CTRL)"),
       y = paste0("-log", "\u2081", "\u2080", "(padj)"),
       fill = "fold change - padj status") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

ggsave(filename = file.path(out_dir_path, "volcano_fc_padj_status.pdf"),
       width = 10,
       height = 10,
       device = cairo_pdf)

# Proteasome --------------------------------------------------------------

temp <- fgsea_sign_down %>% 
  filter(str_detect(set_name, "PROTEASOME"))

### WP_PROTEASOME_DEGRADATION
# colection_name <- "MSigDB_C2_CP:WIKIPATHWAYS"
# gene_set_name <- "WP_PROTEASOME_DEGRADATION"

### KEGG_PROTEASOME
colection_name <- "MSigDB_C2_CP:KEGG_LEGACY"
gene_set_name <- "KEGG_PROTEASOME"


# Violin plot 
gene_set <- gene_metadata %>% 
  select(ensembl_gene_id_version, all_of(colection_name)) %>% 
  unnest(all_of(colection_name)) %>% 
  filter(set_name == gene_set_name) %>% 
  pull(ensembl_gene_id_version) 

leading_edge <- get_leading_edge(input = fgsea_sign_down,
                                 set_name = gene_set_name)

plot_df <- lfc %>%
  filter(!is.na(padj),
         ensembl_gene_id_version %in% gene_set) %>%
  mutate(minus_log10_padj = -log(x = padj, base = 10),
         leading_edge_status = ifelse(ensembl_gene_id_version %in% leading_edge, 
                                      "leading edge", 
                                      "other set members") %>% 
           fct_relevel(c("leading edge", 
                         "other set members")),
         fc_sign_status = case_when(log2FoldChange < 0 & padj < 0.1 ~ "FC < 0 - padj < 0.1",
                                    log2FoldChange > 0 & padj < 0.1 ~ "FC > 0 - padj < 0.1",
                                    padj >= 0.1 ~ paste0("padj ", "\u2265", " 0.1"),
                                    T ~ NA_character_) %>%
           fct_relevel(c(paste0("padj ", "\u2265", " 0.1"),
                         "FC < 0 - padj < 0.1",
                         "FC > 0 - padj < 0.1")),
         label = ifelse(padj < 0.1 & leading_edge_status == "leading edge",
                        gene_symbol,
                        NA_character_))

ggplot(data = plot_df,
       mapping = aes(x = log2FoldChange,
                     y = minus_log10_padj)) +
  geom_vline(xintercept = 0,
             linetype = "solid", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = seq(-1, 1, 0.5),
             linetype = "solid", color = "#403d39", linewidth = 0.1) +
  geom_hline(yintercept = c(-log(x = 0.1, base = 10),
                            -log(x = 0.05, base = 10),
                            -log(x = 0.01, base = 10),
                            -log(x = 0.001, base = 10)),
             linetype = "dashed", color = "#ff9e00", linewidth = 0.3) +
  annotate(geom = "label",
           x = 1,
           y = c(-log(x = 0.1, base = 10),
                 -log(x = 0.05, base = 10),
                 -log(x = 0.01, base = 10),
                 -log(x = 0.001, base = 10)),
           label = c("padj = 0.1",
                     "padj = 0.05",
                     "padj = 0.01",
                     "padj = 0.001"),
           color = "black", size = 3) +
  geom_point(mapping = aes(fill = fc_sign_status, 
                           shape = leading_edge_status, 
                           color = leading_edge_status),
             size = 4, stroke = 0.6, alpha = 0.6) +
  ggrepel::geom_text_repel(data = plot_df %>%
                             filter(!is.na(label)),
                           mapping = aes(label = label),
                           max.overlaps = 500,
                           force = 2,
                           size = 3,
                           box.padding = 0.4,
                           segment.color = "#c0c0c0",
                           segment.size = 0.2) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  scale_y_continuous(n.breaks = 10L) +
  scale_shape_manual(values = c(23, 21)) +
  scale_color_manual(values = c("#780000", "black")) +
  scale_fill_manual(values = c("#efd9ce", "#ca3767", "#5d76cb")) +
  theme_pubr(legend = "bottom") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.direction = "vertical",
        legend.box.just = "left") +
  labs(title = gene_set_name,
       x = paste0("log", "\u2082", " fold change (TRIP vs CTRL)"),
       y = paste0("-log", "\u2081", "\u2080", "(padj)"),
       shape = NULL,
       color = NULL,
       fill = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 1L),
         color = guide_legend(override.aes = list(size = 5), order = 1L),
         fill = guide_legend(override.aes = list(shape = 22, size = 5), order = 2L))

ggsave(filename = file.path(out_dir_path, paste0("volcano_fc_", gene_set_name, ".pdf")),
       width = 9,
       height = 9,
       device = cairo_pdf)

# Log fold change plot for leading edge
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = c(seq(-1, -0.5, 0.5), 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, ".pdf")),
       width = 12,
       height = 2+0.22*dim(lfc %>% 
                             filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Log fold change plot for leading edge with padj < 0.1 
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(padj < 0.1,
                  ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = seq(-1, -0.5, 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, "_sign.pdf")),
       width = 6,
       height = 2+0.22*dim(lfc %>% 
                             filter(padj < 0.1,
                                    ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Proportion plot for leading edge
plot_pct_exp(input = lfc,
             genes = get_leading_edge(input = fgsea_sign_down,
                                      set_name = gene_set_name),
             gene_id = "ensembl_gene_id_version",
             pretty_gene_id = "gene_symbol",
             ref_cond = "CTRL",
             test_cond = "TRIP",
             title = gene_set_name,
             fill_colors = c("#5d76cb", "#ca3767"))

ggsave(filename = file.path(out_dir_path, paste0("pctexp_", gene_set_name, ".pdf")),
       height = 5,
       width = 2+0.3*dim(lfc %>% 
                           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Mitochondria ------------------------------------------------------------

temp <- fgsea_sign_down %>% 
  filter(str_detect(set_name, "OXIDATIVE_PHOSPHORYLATION") | 
           str_detect(set_name, "RESPIRAT"))

### HALLMARK_OXIDATIVE_PHOSPHORYLATION
colection_name <- "MSigDB_H"
gene_set_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

# Violin plot 
gene_set <- gene_metadata %>% 
  select(ensembl_gene_id_version, all_of(colection_name)) %>% 
  unnest(all_of(colection_name)) %>% 
  filter(set_name == gene_set_name) %>% 
  pull(ensembl_gene_id_version) 

leading_edge <- get_leading_edge(input = fgsea_sign_down,
                                 set_name = gene_set_name)

plot_df <- lfc %>%
  filter(!is.na(padj),
         ensembl_gene_id_version %in% gene_set) %>%
  mutate(minus_log10_padj = -log(x = padj, base = 10),
         leading_edge_status = ifelse(ensembl_gene_id_version %in% leading_edge, 
                                      "leading edge", 
                                      "other set members") %>% 
           fct_relevel(c("leading edge", 
                         "other set members")),
         fc_sign_status = case_when(log2FoldChange < 0 & padj < 0.1 ~ "FC < 0 - padj < 0.1",
                                    log2FoldChange > 0 & padj < 0.1 ~ "FC > 0 - padj < 0.1",
                                    padj >= 0.1 ~ paste0("padj ", "\u2265", " 0.1"),
                                    T ~ NA_character_) %>%
           fct_relevel(c(paste0("padj ", "\u2265", " 0.1"),
                         "FC < 0 - padj < 0.1",
                         "FC > 0 - padj < 0.1")),
         label = ifelse(padj < 0.1 & leading_edge_status == "leading edge",
                        gene_symbol,
                        NA_character_))

ggplot(data = plot_df,
       mapping = aes(x = log2FoldChange,
                     y = minus_log10_padj)) +
  geom_vline(xintercept = 0,
             linetype = "solid", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = seq(-2, 1, 0.5),
             linetype = "solid", color = "#403d39", linewidth = 0.1) +
  geom_hline(yintercept = c(-log(x = 0.1, base = 10),
                            -log(x = 0.05, base = 10),
                            -log(x = 0.01, base = 10),
                            -log(x = 0.001, base = 10)),
             linetype = "dashed", color = "#ff9e00", linewidth = 0.3) +
  annotate(geom = "label",
           x = 1,
           y = c(-log(x = 0.1, base = 10),
                 -log(x = 0.05, base = 10),
                 -log(x = 0.01, base = 10),
                 -log(x = 0.001, base = 10)),
           label = c("padj = 0.1",
                     "padj = 0.05",
                     "padj = 0.01",
                     "padj = 0.001"),
           color = "black", size = 3) +
  geom_point(mapping = aes(fill = fc_sign_status, 
                           shape = leading_edge_status, 
                           color = leading_edge_status),
             size = 4, stroke = 0.6, alpha = 0.6) +
  ggrepel::geom_text_repel(data = plot_df %>%
                             filter(!is.na(label)),
                           mapping = aes(label = label),
                           max.overlaps = 500,
                           force = 2,
                           size = 3,
                           box.padding = 0.4,
                           segment.color = "#c0c0c0",
                           segment.size = 0.2) +
  scale_x_continuous(breaks = seq(-2, 1, 0.5)) +
  scale_y_continuous(n.breaks = 10L) +
  scale_shape_manual(values = c(23, 21)) +
  scale_color_manual(values = c("#780000", "black")) +
  scale_fill_manual(values = c("#efd9ce", "#ca3767", "#5d76cb")) +
  theme_pubr(legend = "bottom") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.direction = "vertical",
        legend.box.just = "left") +
  labs(title = gene_set_name,
       x = paste0("log", "\u2082", " fold change (TRIP vs CTRL)"),
       y = paste0("-log", "\u2081", "\u2080", "(padj)"),
       shape = NULL,
       color = NULL,
       fill = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 1L),
         color = guide_legend(override.aes = list(size = 5), order = 1L),
         fill = guide_legend(override.aes = list(shape = 22, size = 5), order = 2L))

ggsave(filename = file.path(out_dir_path, paste0("volcano_fc_", gene_set_name, ".pdf")),
       width = 9,
       height = 9,
       device = cairo_pdf)

# Log fold change plot for leading edge
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = c(seq(-0.5, -2, -0.5), 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, ".pdf")),
       width = 12,
       height = 2+0.22*dim(lfc %>% 
                             filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Log fold change plot for leading edge with padj < 0.1 
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(padj < 0.1,
                  ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = c(seq(-0.5, -2, -0.5), 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, "_sign.pdf")),
       width = 6,
       height = 2+0.22*dim(lfc %>% 
                             filter(padj < 0.1,
                                    ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Proportion plot for leading edge
plot_pct_exp(input = lfc,
             genes = get_leading_edge(input = fgsea_sign_down,
                                      set_name = gene_set_name),
             gene_id = "ensembl_gene_id_version",
             pretty_gene_id = "gene_symbol",
             ref_cond = "CTRL",
             test_cond = "TRIP",
             title = gene_set_name,
             fill_colors = c("#5d76cb", "#ca3767"))

ggsave(filename = file.path(out_dir_path, paste0("pctexp_", gene_set_name, ".pdf")),
       height = 5,
       width = 2+0.3*dim(lfc %>% 
                           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Lipid metabolism --------------------------------------------------------

temp <- fgsea_sign_down %>% 
  filter(str_detect(set_name, "LIPID") | 
           str_detect(set_name, "FATTY_ACID") |
           str_detect(set_name, "CHOLESTEROL"))

# Find gene sets matching the words "LIPID", "FATTY ACID", or "CHOLESTEROL" 
lipid_gene_sets <- fgsea_sign_down %>% 
  filter(str_detect(set_name, "LIPID") | 
           str_detect(set_name, "FATTY_ACID") |
           str_detect(set_name, "CHOLESTEROL")) %>% 
  mutate(matched_words = map_chr(.x = set_name,
                                 .f = function(x) {
                                   map_chr(.x = c("LIPID", "FATTY_ACID", "CHOLESTEROL"),
                                           .f = function(y) {
                                             str_extract(string = x, pattern = y)}) %>% na.omit() %>% paste0(collapse = ", ")
                                 }),
         collection = case_when(collection == "MSigDB_H" ~ "HALLMARK",
                                collection == "MSigDB_C2_CP:KEGG_LEGACY" ~ "KEGG",
                                collection == "MSigDB_C2_CP:WIKIPATHWAYS" ~ "Wikipathways",
                                collection == "MSigDB_C2_CP:REACTOME" ~ "REACTOME",
                                collection == "MSigDB_C2_CP:BIOCARTA" ~ "Biocarta",
                                collection == "MSigDB_C5_GO:BP" ~ "GO BP",
                                T ~ collection) %>% 
           fct_relevel(c("HALLMARK")),
         padj = as.character(format(padj, scientific = T, digits = 2)),
         NES = as.character(round(NES, 2)),
         pct_leadingEdge = as.character(round(pct_leadingEdge, 2))) %>% 
  arrange(collection, padj) %>% 
  select(collection, 
         `set name` = set_name, 
         `adj. p` = padj, 
         NES, 
         `leading edge size / gene set size` = n_leadingEdge_over_total, 
         `% leading edge` = pct_leadingEdge,
         `matched words` = matched_words) 

lipid_gene_sets %>% 
  kable() %>%
  kable_classic(full_width = T) %>%
  kable_styling(bootstrap_options = "striped", font_size = 25, position = "left", full_width = F) %>% 
  save_kable(file = file.path(out_dir_path, "lipid_gene_sets_pretty_table.html"))

### HALLMARK_FATTY_ACID_METABOLISM
colection_name <- "MSigDB_H"
gene_set_name <- "HALLMARK_FATTY_ACID_METABOLISM"

# Violin plot 
gene_set <- gene_metadata %>% 
  select(ensembl_gene_id_version, all_of(colection_name)) %>% 
  unnest(all_of(colection_name)) %>% 
  filter(set_name == gene_set_name) %>% 
  pull(ensembl_gene_id_version) 

leading_edge <- get_leading_edge(input = fgsea_sign_down,
                                 set_name = gene_set_name)

plot_df <- lfc %>%
  filter(!is.na(padj),
         ensembl_gene_id_version %in% gene_set) %>%
  mutate(minus_log10_padj = -log(x = padj, base = 10),
         leading_edge_status = ifelse(ensembl_gene_id_version %in% leading_edge, 
                                      "leading edge", 
                                      "other set members") %>% 
           fct_relevel(c("leading edge", 
                         "other set members")),
         fc_sign_status = case_when(log2FoldChange < 0 & padj < 0.1 ~ "FC < 0 - padj < 0.1",
                                    log2FoldChange > 0 & padj < 0.1 ~ "FC > 0 - padj < 0.1",
                                    padj >= 0.1 ~ paste0("padj ", "\u2265", " 0.1"),
                                    T ~ NA_character_) %>%
           fct_relevel(c(paste0("padj ", "\u2265", " 0.1"),
                         "FC < 0 - padj < 0.1",
                         "FC > 0 - padj < 0.1")),
         label = ifelse(padj < 0.1 & leading_edge_status == "leading edge",
                        gene_symbol,
                        NA_character_))

ggplot(data = plot_df,
       mapping = aes(x = log2FoldChange,
                     y = minus_log10_padj)) +
  geom_vline(xintercept = 0,
             linetype = "solid", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = seq(-2, 1.5, 0.5),
             linetype = "solid", color = "#403d39", linewidth = 0.1) +
  geom_hline(yintercept = c(-log(x = 0.1, base = 10),
                            -log(x = 0.05, base = 10),
                            -log(x = 0.01, base = 10),
                            -log(x = 0.001, base = 10)),
             linetype = "dashed", color = "#ff9e00", linewidth = 0.3) +
  annotate(geom = "label",
           x = 1.5,
           y = c(-log(x = 0.1, base = 10),
                 -log(x = 0.05, base = 10),
                 -log(x = 0.01, base = 10),
                 -log(x = 0.001, base = 10)),
           label = c("padj = 0.1",
                     "padj = 0.05",
                     "padj = 0.01",
                     "padj = 0.001"),
           color = "black", size = 3) +
  geom_point(mapping = aes(fill = fc_sign_status, 
                           shape = leading_edge_status, 
                           color = leading_edge_status),
             size = 4, stroke = 0.6, alpha = 0.6) +
  ggrepel::geom_text_repel(data = plot_df %>%
                             filter(!is.na(label)),
                           mapping = aes(label = label),
                           max.overlaps = 500,
                           force = 2,
                           size = 3,
                           box.padding = 0.4,
                           segment.color = "#c0c0c0",
                           segment.size = 0.2) +
  scale_x_continuous(breaks = seq(-2, 1.5, 0.5)) +
  scale_y_continuous(n.breaks = 10L) +
  scale_shape_manual(values = c(23, 21)) +
  scale_color_manual(values = c("#780000", "black")) +
  scale_fill_manual(values = c("#efd9ce", "#ca3767", "#5d76cb")) +
  theme_pubr(legend = "bottom") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.direction = "vertical",
        legend.box.just = "left") +
  labs(title = gene_set_name,
       x = paste0("log", "\u2082", " fold change (TRIP vs CTRL)"),
       y = paste0("-log", "\u2081", "\u2080", "(padj)"),
       shape = NULL,
       color = NULL,
       fill = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 5), order = 1L),
         color = guide_legend(override.aes = list(size = 5), order = 1L),
         fill = guide_legend(override.aes = list(shape = 22, size = 5), order = 2L))

ggsave(filename = file.path(out_dir_path, paste0("volcano_fc_", gene_set_name, ".pdf")),
       width = 9,
       height = 9,
       device = cairo_pdf)

# Log fold change plot for leading edge
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = c(seq(-2, -0.5, 0.5), 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, ".pdf")),
       width = 12,
       height = 2+0.22*dim(lfc %>% 
                             filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Log fold change plot for leading edge with padj < 0.1 
plot_lfc(input = lfc,
         genes = lfc %>% 
           filter(padj < 0.1,
                  ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                set_name = gene_set_name)) %>% 
           pull(ensembl_gene_id_version),
         gene_id = ensembl_gene_id_version,
         pretty_gene_id = gene_symbol,
         ref_cond = "CTRL",
         test_cond = "TRIP",
         key_genes = NULL,
         title = gene_set_name,
         desc_lfc = FALSE,
         xintercepts = c(seq(-2, -0.5, 0.5), 0.5),
         padj_lab = TRUE,
         lab_nudge = 0.02)

ggsave(filename = file.path(out_dir_path, paste0("lfc_", gene_set_name, "_sign.pdf")),
       width = 8,
       height = 2+0.22*dim(lfc %>% 
                             filter(padj < 0.1,
                                    ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                  set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Proportion plot for leading edge
plot_pct_exp(input = lfc,
             genes = get_leading_edge(input = fgsea_sign_down,
                                      set_name = gene_set_name),
             gene_id = "ensembl_gene_id_version",
             pretty_gene_id = "gene_symbol",
             ref_cond = "CTRL",
             test_cond = "TRIP",
             title = gene_set_name,
             fill_colors = c("#5d76cb", "#ca3767"))

ggsave(filename = file.path(out_dir_path, paste0("pctexp_", gene_set_name, ".pdf")),
       height = 5,
       width = 2+0.3*dim(lfc %>% 
                           filter(ensembl_gene_id_version %in% get_leading_edge(input = fgsea_sign_down,
                                                                                set_name = gene_set_name)))[[1]],
       device = cairo_pdf)

# Key gene sets -----------------------------------------------------------

key_gene_sets <- fgsea_sign_down %>% 
  filter(set_name %in% c(
    "WP_PROTEASOME_DEGRADATION",
    "KEGG_PROTEASOME",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "GOBP_LIPID_CATABOLIC_PROCESS"
  )) %>% 
  mutate(collection = case_when(collection == "MSigDB_H" ~ "HALLMARK",
                                collection == "MSigDB_C2_CP:KEGG_LEGACY" ~ "KEGG",
                                collection == "MSigDB_C2_CP:WIKIPATHWAYS" ~ "Wikipathways",
                                collection == "MSigDB_C2_CP:REACTOME" ~ "REACTOME",
                                collection == "MSigDB_C2_CP:BIOCARTA" ~ "Biocarta",
                                collection == "MSigDB_C5_GO:BP" ~ "GO BP",
                                T ~ collection) %>% 
           fct_relevel(c("HALLMARK")),
         padj = as.character(format(padj, scientific = T, digits = 2)),
         NES = as.character(round(NES, 2)),
         pct_leadingEdge = as.character(round(pct_leadingEdge, 2))) %>% 
  arrange(collection, padj) %>% 
  select(collection, 
         `set name` = set_name, 
         `adj. p` = padj, 
         NES, 
         `leading edge size / gene set size` = n_leadingEdge_over_total, 
         `% leading edge` = pct_leadingEdge) 

key_gene_sets %>% 
  kable() %>%
  kable_classic(full_width = T) %>%
  kable_styling(bootstrap_options = "striped", font_size = 25, position = "left", full_width = F) %>% 
  save_kable(file = file.path(out_dir_path, paste0("key_gene_sets_pretty_table.html")))

# single cell level TPMs --------------------------------------------------

cell_ids <- cell_metadata %>% 
  select(line_name, barcode) %>% 
  arrange(line_name) %>% 
  group_by(line_name) %>% 
  mutate(line_name_cell_id = paste0(line_name, "_cell", seq_along(.data$line_name))) %>% 
  ungroup()

tpm_pretty <- tpm %>% 
  column_to_rownames(var = "ensembl_gene_id_version")

tpm_pretty <- tpm_pretty[, cell_ids$barcode]

identical(tpm_pretty %>% colnames(), cell_ids$barcode)

colnames(tpm_pretty) <- cell_ids$line_name_cell_id

tpm_pretty <- rownames_to_column(tpm_pretty,
                                 var = "ensembl_gene_id_version") %>% 
  left_join(gene_metadata[, c("ensembl_gene_id_version", "gene_symbol")],
            by = join_by(ensembl_gene_id_version)) %>% 
  select(ensembl_gene_id_version, gene_symbol, everything())

write_csv(tpm_pretty,
          file.path(out_dir_path, "tpm_pretty.csv"))

# Significant GSEA --------------------------------------------------------

fgsea_sign_down_res <- fgsea_sign_down %>% 
  select(collection, set_name, pval, padj, NES, set_size = size, n_contrib_genes = n_leadingEdge) %>% 
  arrange(collection, padj)

fgsea_sign_up_res <- fgsea_sign_up %>% 
  select(collection, set_name, pval, padj, NES, set_size = size, n_contrib_genes = n_leadingEdge) %>% 
  arrange(collection, padj)

write_csv(fgsea_sign_down_res,
          file.path(out_dir_path, "fgsea_sign_down_res.csv"))

write_csv(fgsea_sign_up_res,
          file.path(out_dir_path, "fgsea_sign_up_res.csv"))

