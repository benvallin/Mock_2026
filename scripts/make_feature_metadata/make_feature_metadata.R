# File description --------------------------------------------------------

### Goal
# make_feature_metadata.R produces gene and transcript metadata from GENECODE annotation.

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# ===>>> feature_metadata/make_feature_metadata.R
# - count_processing/summarise_to_gene_level.R
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
# - rnaseqtools package available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
work_dir_path <- "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path <- "publications/Mock_2026/data/output/make_feature_metadata"

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

# Add cell types from Jerber et al. ---------------------------------------

# List cell types from Jerber et al. (PMID: 33664506)
proliferating = c("MKI67", "TOP2A", "CIP2A") # CIP2A == KIAA1524 in original publication

neural_prog = c("LIX1", "RAX", "NR2F1", "NES")

fp_prog = c("ZEB2", "DMBX1","HMGA1", "HMGB2")

glia = c("TNC", "SOX2", "CDH2", "HES1")

fp = c("LMX1A", "FOXA2")

astrocyte = c("S100B", "AQP4", "GFAP", "SLC1A3", "SOX9")

ependymal = c("STOML3", "CCDC153", "CDHR4", "FOXJ1", "DNAH11", "TTR", "MLF1")

serotonergic = c("CHGB", "DDC","FEV","GATA2", "GATA3", "GCH1", "GCHFR", "HTR1A",
                 "HTR1B","MAOA","MAOB","SLC18A2","SLC29A4","SLC6A4","TPH2")

neuron = c("SYT1", "SNAP25")

neuroblast = c("NEUROG1", "NEUROD1", "NEUROG2", "NHLH1", "SIM1")

dopaminergic = c("ABCC8","ACOT7","ALDH1A1","AMER3","ARG2","ASB4","BNC2","CADPS2","CALB1","CALB2","CAMK2N1",
                 "CCK","CDK14","CDKN1C","CHL1", "CHRNA4","CPEB3","CXCR4","DCC","DKK3","DRD2",
                 "EBF2","EN1","EN2","EPHA5","ERC2","FGF13","FOXA1","FOXA2","GDAP1","GFRA1","GRIA3","GRIK3",
                 "GRP","ICA1L","IGF1","KCNIP4","KCNJ6","KLHL1","KLHL13","LGI1","LMO3","LMX1A","LMX1B",
                 "LRRC3B","LRRTM2","LSAMP","LY6H","NETO2","NR4A2","NTSR1","OTX2","PBX1","PRICKLE2","PRKCA",
                 "PRL","PRRT4","PITX3","PTPN5","PTPRO","RET","SCG2","SLC10A4","SLC18A1","SLC18A2","SLC18A3",
                 "SLC6A3","SNCA","SOX6","TH","TMCC3","TMEFF2","TMEM255A","TUB","VGF")

cortical_hem <- c("DKK3","EOMES","RELN","LHX1","LHX5", "TP73", "CALB2","EBF3","SAMD3","SOX2")

# Check that all gene names are present in gene metadata
walk(.x = c("proliferating", "neural_prog", "fp_prog", "glia", "fp", "astrocyte", "ependymal",
            "serotonergic", "neuron", "neuroblast", "dopaminergic", "cortical_hem"),
     .f = ~ all(print(paste0(.x, ": ", all(eval(sym(.x)) %in% gene_metadata$gene_symbol)))))

# Add cell type annotations to gene metadata
gene_metadata <- gene_metadata %>%
  mutate(jerber_marker = ifelse(gene_symbol %in% unique(c(proliferating, neural_prog, fp_prog, glia, fp, astrocyte, ependymal,
                                                          serotonergic, neuron, neuroblast, dopaminergic, cortical_hem)),
                                T, F),
         proliferating_marker = ifelse(gene_symbol %in% proliferating, T, F),
         neural_prog_marker = ifelse(gene_symbol %in% neural_prog, T, F),
         fp_prog_marker = ifelse(gene_symbol %in% fp_prog, T, F),
         glia_marker = ifelse(gene_symbol %in% glia, T, F),
         fp_marker = ifelse(gene_symbol %in% fp, T, F),
         astrocyte_marker = ifelse(gene_symbol %in% astrocyte, T, F),
         ependymal_marker = ifelse(gene_symbol %in% ependymal, T, F),
         serotonergic_marker = ifelse(gene_symbol %in% serotonergic, T, F),
         neuron_marker = ifelse(gene_symbol %in% neuron, T, F),
         neuroblast_marker = ifelse(gene_symbol %in% neuroblast, T, F),
         dopaminergic_marker = ifelse(gene_symbol %in% dopaminergic, T, F),
         cortical_hem_marker = ifelse(gene_symbol %in% cortical_hem, T, F))

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
