# File description --------------------------------------------------------

### Goal
# build_adata.py constructs the AnnData object for the scRNA-seq datataset.

### Execution steps
# x make_feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# ===>>> cell_clustering/build_adata.py
# - cell_clustering/integration.py
# - update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - dge_gsea/summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - cell_metadata.parquet file produced by mark_barcode_quality.R
# - gene_metadata.csv file produced by make_feature_metadata.R
# - Gene-level estimated counts produced by summarise_to_gene_level.R
# - py_functions.py defining helper python functions
# - paramstream package available from PyPI with 'pip install py-paramstream'

# => Provided the above requirement(s) are met, one only needs to adjust the paths in the section below to run this script <=

### Author
# Benjamin Vallin (benjamin.vallin@dpag.ox.ac.uk)

# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
work_dir_path = "/scratch/ben/rnaseq"

# Path to output directory
out_dir_path = "publications/Mock_2026/data/cell_clustering/build_adata"

# Path to directory containing Ben's python functions file (py_functions.py)
py_funs_dir_path = "tools/rnaseq_functions"

# Path to directory containing counts files
cnt_dir_path = (
    "publications/Mock_2026/data/count_processing/summarise_to_gene_level"
)

# Path to directory containing cell metadata file
cell_metadata_dir_path = (
    "publications/Mock_2026/data/update_cell_metadata/mark_barcode_quality"
)

gene_metadata_path = (
    "publications/Mock_2026/data/make_feature_metadata/v44/gene_metadata.csv"
)

# %% Set up ----

# Import required libraries
import os, sys
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import paramstream as ps

# Set working directory
os.chdir(work_dir_path)

# Load Ben's python functions
sys.path.append(py_funs_dir_path)
sf = __import__("py_functions")

# Create output directory if necessary
Path(out_dir_path).mkdir(parents=True, exist_ok=True)

# %% Define analysis parameters ----

# Define input analysis parameters
input_params = {
    "min_libsize_incl": 20000,
    "min_gene_incl": 3000,
    "max_pct_mt_incl": 20,
}

# Define current analysis parameters
current_params = {"n_cnt_threshold": 5, "pct_cell_threshold": 0.01}

# %% Stream analysis parameter values ----


input_param_log_id = ps.get_param_log_id(
    param_dict=input_params, in_dir_path=cell_metadata_dir_path
)

output_params = {**input_params, **current_params}

ps.write_param_log(
    param_dict=output_params, out_dir_path=out_dir_path, param_log_table=True
)

output_param_log_id = ps.get_param_log_id(
    param_dict=output_params, in_dir_path=out_dir_path
)

# %% Build adata ----

# Construct adata with raw count matrix
adata = sc.read_csv(filename=Path(cnt_dir_path, "cnt_lengthScaledTPM.csv")).transpose()

# Convert raw count matrix to sparse
adata.X = sparse.csr_matrix(adata.X)

# Set raw count layer
adata.layers["cnt"] = adata.X.copy()

# Add TPM layer
adata.layers["tpm"] = sparse.csr_matrix(
    sc.read_csv(filename=Path(cnt_dir_path, "tpm_lengthScaledTPM.csv")).transpose().X
)

# Compute log2 TPM+1
adata.layers["log2_tpm1p"] = sc.pp.log1p(adata.layers["tpm"].copy(), base=2)

# %% Add cell metadata ----

# Load cell metadata
cell_metadata = pd.read_parquet(
    Path(cell_metadata_dir_path, f"cell_metadata_{input_param_log_id}.parquet")
)

bool_cols = cell_metadata.select_dtypes(include="object").columns
bool_cols = [
    i for i in bool_cols if cell_metadata[i].dropna().isin([True, False]).all()
]
cell_metadata[bool_cols] = cell_metadata[bool_cols].apply(lambda x: x.astype("bool"))

cell_metadata = cell_metadata.set_index(keys="barcode", drop=False)

# Subset to cells passing barcode QC
cell_metadata = cell_metadata.loc[cell_metadata["pass_qc"] == True, :]

# Subset adata to barcodes in cell metadata
adata = adata[cell_metadata["barcode"], :]

# Set adata.obs to cell metadata
adata.obs = cell_metadata

# %% Add gene metadata ----

# Load gene metadata
gene_metadata = pd.read_csv(gene_metadata_path)
gene_metadata = gene_metadata.set_index(keys="ensembl_gene_id_version", drop=False)

# Filter/reorder rows to match adata.var_names
gene_metadata = gene_metadata.loc[adata.var_names]

# Set adata.var to gene metadata
adata.var = gene_metadata

# %% Filter adata.var ----

# Subset to genes with at least n_cnt_threshold raw counts in at least pct_cell_threshold of the cells
adata = adata[
    :,
    (adata.X >= current_params["n_cnt_threshold"]).sum(axis=0)
    > (current_params["pct_cell_threshold"] * adata.n_obs),
].copy()

# %% Tidy adata.obs ----

# Set key obs variables as categories and supply colors for plotting

# icell8_chip_id
icell8_chip_id_colors = {
    "A": "#ffbe0c",
    "B": "#f0396a",
}

sf.set_pretty_cat(adata=adata, obsnm="icell8_chip_id", coldict=icell8_chip_id_colors)

# line_name
line_name_colors = {
    "B856": "#3966ce",
    "B156": "#ff6d00",
    "B067": "#0b8c1f",
    "831-5": "#f0396a",
    "831-1": "#642ca9",
    "831-3": "#ffbe0c",
}

sf.set_pretty_cat(adata=adata, obsnm="line_name", coldict=line_name_colors)

# line_name_generic
line_name_generic_colors = {
    "CTRL1": "#3966ce",
    "CTRL2": "#ff6d00",
    "CTRL3": "#0b8c1f",
    "TRIP1": "#f0396a",
    "TRIP2": "#642ca9",
    "TRIP3": "#ffbe0c",
}

sf.set_pretty_cat(
    adata=adata, obsnm="line_name_generic", coldict=line_name_generic_colors
)

# genotype
genotype_colors = {
    "CTRL": "#0077b6",
    "TRIP": "#a4133c",
}

sf.set_pretty_cat(adata=adata, obsnm="genotype", coldict=genotype_colors)

# donor_sex
donor_sex_colors = {"F": "#5d76cb", "M": "#ca3767"}

sf.set_pretty_cat(adata=adata, obsnm="donor_sex", coldict=donor_sex_colors)


# %% Write adata to disk ----

# Write filtered adata to disk
list_cols = [
    i
    for i in adata.obs.columns
    if adata.obs[i].apply(lambda x: isinstance(x, (list, tuple, np.ndarray))).any()
]

for i in list_cols:
    adata.obs[i] = adata.obs[i].apply(
        lambda x: (
            ",".join(map(str, x)) if isinstance(x, (list, tuple, np.ndarray)) else x
        )
    )

adata.write_h5ad(filename=Path(out_dir_path, f"adata_{output_param_log_id}.h5ad"))

# %%
