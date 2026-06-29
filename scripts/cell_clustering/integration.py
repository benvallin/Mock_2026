# File description --------------------------------------------------------

### Goal
# integration.py performs data integration, cell clustering and cell type annotation.

### Execution steps
# x make_cell_metadata/make_cell_metadata.R
# x feature_metadata/make_feature_metadata.R
# x count_processing/summarise_to_gene_level.R
# x update_cell_metadata/add_count_metrics.R
# x update_cell_metadata/mark_barcode_quality.R
# x cell_clustering/build_adata.py
# ===>>> cell_clustering/integration.py
# - update_cell_metadata/add_cell_type_annotation.R
# - cell_clustering/summarise_cell_counts.R
# - dge_gsea/dge_mast.R
# - dge_gsea/gsea_dge_mast.R
# - summarise_gsea_dge_mast_results.R

# => Any script listed before the current file must be executed beforehand

### Requirements
# - adata.h5ad produced by build_adata.py
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
out_dir_path = "publications/Mock_2026/data/output/cell_clustering/integration"

# Path to directory containing Ben's python functions file (py_functions.py)
py_funs_dir_path = "tools/rnaseq_functions"

# Path to directory containing adata file
adata_dir_path = "publications/Mock_2026/data/output/cell_clustering/build_adata"

# %% Set up ----

# Import required libraries
import os, sys, scvi
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import paramstream as ps
import torch

# Set working directory
os.chdir(work_dir_path)

# Load Ben's python functions
sys.path.append(py_funs_dir_path)
sf = __import__("py_functions")

# Create output directory if necessary
Path(out_dir_path).mkdir(parents=True, exist_ok=True)

# Set scanpy output directory
sc.settings.figdir = out_dir_path

# Set number of workers for scVI train_dataloader
scvi.settings.dl_num_workers = 8

# Set scVI seed for reproducibility
scvi.settings.seed = 27

# Set seed for other stochastic steps
random_state = 27

# Force CUDA to use deterministic algorithms
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

# %% Define analysis parameters ----

# Define donor annotation + barcode QC parameters (input)
input_params = {
    "min_libsize_incl": 20000,
    "min_gene_incl": 3000,
    "max_pct_mt_incl": 20,
}

# Define current analysis parameters (current)
current_params = {
    "n_top_genes": 2000,
    "cat_cov_key": "line_name_generic",
    "integration_method": "scvi",
}

# Define input + current parameters (output)
output_params = {**input_params, **current_params}

# %% Stream analysis parameter values ----

input_param_log_id = ps.get_param_log_id(
    param_dict=input_params, in_dir_path=adata_dir_path
)

ps.write_param_log(
    param_dict=output_params, out_dir_path=out_dir_path, param_log_table=True
)

output_param_log_id = ps.get_param_log_id(
    param_dict=output_params, in_dir_path=out_dir_path
)

# %% Load and filter adata ----

# Load filtered adata
adata = sc.read_h5ad(filename=Path(adata_dir_path, f"adata_{input_param_log_id}.h5ad"))

# Provide cell_type colors for future annotations
for i in [
    "cell_type_unint",
    "cell_type_scvi",
    "cell_type",
]:

    adata.uns[f"{i}_colors"] = ["#5d76cb", "#ca3767", "#29a655"]

for i in [
    "cell_type_bbknn",
    "cell_type_harmony",
]:

    adata.uns[f"{i}_colors"] = ["#5d76cb", "#ca3767"]

# %% Compute PCA on HVGs ----

# Find highly variable genes (HVG) using flavor 'seurat_v3' on raw counts
sc.pp.highly_variable_genes(
    adata=adata,
    flavor="seurat_v3",
    layer="cnt",
    n_top_genes=current_params["n_top_genes"],
    batch_key=current_params["cat_cov_key"],
)

# Set X to log2 TPM+1
adata.X = adata.layers["log2_tpm1p"].copy()

# Z-score transform log2(TPM+1) values (scaling to unit variance and zero mean)
sc.pp.scale(adata)

# Run PCA on Z-score transformed log2(TPM+1) values
sc.pp.pca(data=adata, mask_var="highly_variable", random_state=random_state)

# Reset X to raw counts
adata.X = adata.layers["cnt"].copy()

# %% Find the 'top X% marker gene-expressing cells' ----

# Define marker_genes
marker_genes = [
    "SNAP25",
    "SYT1",
    "TUBB3",
    "TH",
    "SNCA",
    "AQP4",
    "mTagBFP2",
    "puromycinR",
]

# Define top_pct
top_pct = 25

# Find the top marker gene-expressing cells based on log2 TPM+1 values
for i in marker_genes:
    sf.find_top_expressing_cells(
        adata=adata, layer="log2_tpm1p", gene_symbol=i, top_pct=top_pct
    )

# %% Cluster unintegrated data ----

# Compute neighbors
sc.pp.neighbors(adata=adata, use_rep="X_pca", key_added="neig_unint")

# Compute UMAP
sc.tl.umap(
    adata=adata,
    neighbors_key="neig_unint",
    key_added="umap_unint",
    random_state=random_state,
)

sc.pl.embedding(
    adata=adata,
    basis="umap_unint",
    neighbors_key="neig_unint",
    color=[
        "line_name",
        "donor_sex",
        "icell8_chip_id",
        "libsize",
        "pct_mt_genes",
        "pct_rb_genes",
    ],
    title=[
        "line ID",
        "donor sex",
        "iCELL8 chip ID",
        "total counts",
        "% mitochondrial counts",
        "% ribosomal counts",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=3,
    legend_fontsize="medium",
)

sc.pl.embedding(
    adata=adata,
    basis="umap_unint",
    neighbors_key="neig_unint",
    color=[
        "line_name",
        "genotype",
        "donor_sex",
        "icell8_chip_id",
    ],
    title=[
        "line ID",
        "genotype",
        "donor sex",
        "iCELL8 chip ID",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    legend_fontsize="medium",
)

# Find clusters
leiden_key = "leid_unint_"
for i in np.round(np.arange(0.1, 1.1, 0.1), 1):
    sc.tl.leiden(
        adata=adata,
        neighbors_key="neig_unint",
        resolution=i,
        key_added=f"{leiden_key}{str(i)}",
        random_state=random_state,
    )

# sc.pl.embedding(
#     adata=adata,
#     basis="umap_unint",
#     neighbors_key="neig_unint",
#     color=[
#         "line_name",
#         f"{leiden_key}0.1",
#         "SNAP25",
#         "top25pct_SNAP25_expressing",
#         "TH",
#         "top25pct_TH_expressing",
#         "AQP4",
#         "top25pct_AQP4_expressing",
#     ],
#     title=[
#         "line ID",
#         "leiden unint. 0.1",
#         "SNAP25",
#         "top 25% SNAP25-expressing",
#         "TH",
#         "top 25% TH-expressing",
#         "AQP4",
#         "top 25% AQP4-expressing",
#     ],
#     layer="tpm",
#     sort_order=True,
#     frameon=False,
#     edges=False,
#     cmap="Reds",
#     ncols=4,
#     gene_symbols="gene_symbol",
#     vmin=0,
#     vmax="p99",
#     legend_fontsize="medium",
# )

# Define cell types
cell_type = "cell_type_unint"
reso = 0.1
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["0", "2"]), "neuron", "other"
)
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["1"]), "astrocyte", adata.obs[cell_type]
)
adata.obs[cell_type] = adata.obs[cell_type].astype("category")
adata.obs[cell_type] = adata.obs[cell_type].cat.reorder_categories(
    ["neuron", "astrocyte", "other"]
)

sc.pl.embedding(
    adata=adata,
    basis="umap_unint",
    neighbors_key="neig_unint",
    color=[
        f"{leiden_key}{reso}",
        cell_type,
        "SNAP25",
        "top25pct_SNAP25_expressing",
        "TH",
        "top25pct_TH_expressing",
        "AQP4",
        "top25pct_AQP4_expressing",
    ],
    title=[
        f"leiden unint. {reso}",
        "cell type unint.",
        "SNAP25",
        "top 25% SNAP25-expressing",
        "TH",
        "top 25% TH-expressing",
        "AQP4",
        "top 25% AQP4-expressing",
    ],
    layer="tpm",
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    legend_fontsize="medium",
)

# %% Cluster BBKNN-integrated data ----

# Perform BBKNN integration
temp = sc.external.pp.bbknn(
    adata, batch_key=current_params["cat_cov_key"], use_rep="X_pca", copy=True
)

# Transfer integration results to adata
adata.uns["neig_bbknn"] = temp.uns["neighbors"].copy()

for i in ["distances_key", "connectivities_key"]:
    adata.uns["neig_bbknn"][i] = f"neig_bbknn_{adata.uns['neig_bbknn'][i]}"

adata.obsp["neig_bbknn_distances"] = temp.obsp["distances"].copy()

adata.obsp["neig_bbknn_connectivities"] = temp.obsp["connectivities"].copy()

del temp

# Compute UMAP
sc.tl.umap(
    adata=adata,
    neighbors_key="neig_bbknn",
    key_added="umap_bbknn",
    random_state=random_state,
)

sc.pl.embedding(
    adata=adata,
    basis="umap_bbknn",
    neighbors_key="neig_bbknn",
    color=[
        "line_name",
        "donor_sex",
        "icell8_chip_id",
        "libsize",
        "pct_mt_genes",
        "pct_rb_genes",
    ],
    title=[
        "line ID",
        "donor sex",
        "iCELL8 chip ID",
        "total counts",
        "% mitochondrial counts",
        "% ribosomal counts",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=3,
    legend_fontsize="medium",
)

sc.pl.embedding(
    adata=adata,
    basis="umap_bbknn",
    neighbors_key="neig_bbknn",
    color=[
        "line_name",
        "genotype",
        "donor_sex",
        "icell8_chip_id",
    ],
    title=[
        "line ID",
        "genotype",
        "donor sex",
        "iCELL8 chip ID",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    legend_fontsize="medium",
)

# Find clusters
leiden_key = "leid_bbknn_"
for i in np.round(np.arange(0.1, 2.1, 0.1), 1):
    sc.tl.leiden(
        adata=adata,
        neighbors_key="neig_bbknn",
        resolution=i,
        key_added=f"{leiden_key}{str(i)}",
        random_state=random_state,
    )

# sc.pl.embedding(
#     adata=adata,
#     basis="umap_bbknn",
#     neighbors_key="neig_bbknn",
#     color=[
#         "line_name",
#         f"{leiden_key}0.1",
#         "SNAP25",
#         "top25pct_SNAP25_expressing",
#         "TH",
#         "top25pct_TH_expressing",
#         "AQP4",
#         "top25pct_AQP4_expressing",
#     ],
#     title=[
#         "line ID",
#         "leiden BBKNN 0.1",
#         "SNAP25",
#         "top 25% SNAP25-expressing",
#         "TH",
#         "top 25% TH-expressing",
#         "AQP4",
#         "top 25% AQP4-expressing",
#     ],
#     layer="tpm",
#     sort_order=True,
#     frameon=False,
#     edges=False,
#     cmap="Reds",
#     ncols=4,
#     gene_symbols="gene_symbol",
#     vmin=0,
#     vmax="p99",
#     legend_fontsize="medium",
# )

# Define cell types
cell_type = "cell_type_bbknn"
reso = 0.1
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["0", "2"]), "neuron", "astrocyte / other"
)
adata.obs[cell_type] = adata.obs[cell_type].astype("category")
adata.obs[cell_type] = adata.obs[cell_type].cat.reorder_categories(
    ["neuron", "astrocyte / other"]
)

sc.pl.embedding(
    adata=adata,
    basis="umap_bbknn",
    neighbors_key="neig_bbknn",
    color=[
        f"{leiden_key}{reso}",
        cell_type,
        "SNAP25",
        "top25pct_SNAP25_expressing",
        "TH",
        "top25pct_TH_expressing",
        "AQP4",
        "top25pct_AQP4_expressing",
    ],
    title=[
        f"leiden BBKNN {reso}",
        "cell type BBKNN",
        "SNAP25",
        "top 25% SNAP25-expressing",
        "TH",
        "top 25% TH-expressing",
        "AQP4",
        "top 25% AQP4-expressing",
    ],
    layer="tpm",
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    legend_fontsize="medium",
)

# %% Cluster Harmony-integrated data ----

# Perform Harmony integration
sc.external.pp.harmony_integrate(
    adata,
    key=current_params["cat_cov_key"],
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
    random_state=random_state,
)

# Compute neighbors using the Harmony-corrected embedding
sc.pp.neighbors(adata=adata, use_rep="X_pca_harmony", key_added="neig_harmony")

# Compute UMAP
sc.tl.umap(
    adata=adata,
    neighbors_key="neig_harmony",
    key_added="umap_harmony",
    random_state=random_state,
)

sc.pl.embedding(
    adata=adata,
    basis="umap_harmony",
    neighbors_key="neig_harmony",
    color=[
        "line_name",
        "donor_sex",
        "icell8_chip_id",
        "libsize",
        "pct_mt_genes",
        "pct_rb_genes",
    ],
    title=[
        "line ID",
        "donor sex",
        "iCELL8 chip ID",
        "total counts",
        "% mitochondrial counts",
        "% ribosomal counts",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=3,
    legend_fontsize="medium",
)

sc.pl.embedding(
    adata=adata,
    basis="umap_harmony",
    neighbors_key="neig_harmony",
    color=[
        "line_name",
        "genotype",
        "donor_sex",
        "icell8_chip_id",
    ],
    title=[
        "line ID",
        "genotype",
        "donor sex",
        "iCELL8 chip ID",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    legend_fontsize="medium",
)

# Find clusters
leiden_key = "leid_harmony_"
for i in np.round(np.arange(0.1, 1.1, 0.1), 1):
    sc.tl.leiden(
        adata=adata,
        neighbors_key="neig_harmony",
        resolution=i,
        key_added=f"{leiden_key}{str(i)}",
        random_state=random_state,
    )

# sc.pl.embedding(
#     adata=adata,
#     basis="umap_harmony",
#     neighbors_key="neig_harmony",
#     color=[
#         "line_name",
#         f"{leiden_key}0.1",
#         "SNAP25",
#         "top25pct_SNAP25_expressing",
#         "TH",
#         "top25pct_TH_expressing",
#         "AQP4",
#         "top25pct_AQP4_expressing",
#     ],
#     title=[
#         "line ID",
#         "leiden Harmony 0.1",
#         "SNAP25",
#         "top 25% SNAP25-expressing",
#         "TH",
#         "top 25% TH-expressing",
#         "AQP4",
#         "top 25% AQP4-expressing",
#     ],
#     layer="tpm",
#     sort_order=True,
#     frameon=False,
#     edges=False,
#     cmap="Reds",
#     ncols=4,
#     gene_symbols="gene_symbol",
#     vmin=0,
#     vmax="p99",
#     legend_fontsize="medium",
# )

# Define cell types
cell_type = "cell_type_harmony"
reso = 0.1
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["0", "2", "3"]),
    "neuron",
    "astrocyte / other",
)
adata.obs[cell_type] = adata.obs[cell_type].astype("category")
adata.obs[cell_type] = adata.obs[cell_type].cat.reorder_categories(
    ["neuron", "astrocyte / other"]
)

sc.pl.embedding(
    adata=adata,
    basis="umap_harmony",
    neighbors_key="neig_harmony",
    color=[
        f"{leiden_key}{reso}",
        cell_type,
        "SNAP25",
        "top25pct_SNAP25_expressing",
        "TH",
        "top25pct_TH_expressing",
        "AQP4",
        "top25pct_AQP4_expressing",
    ],
    title=[
        f"leiden Harmony {reso}",
        "cell type Harmony",
        "SNAP25",
        "top 25% SNAP25-expressing",
        "TH",
        "top 25% TH-expressing",
        "AQP4",
        "top 25% AQP4-expressing",
    ],
    layer="tpm",
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    legend_fontsize="medium",
)

# %% Cluster scVI-integrated data ----

# Prepare adata for scVI integration
# => Create a copy of adata containing only the HVG
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# => Set up the HVG-only adata
scvi.model.SCVI.setup_anndata(
    adata=adata_hvg,
    layer="cnt",
    categorical_covariate_keys=[current_params["cat_cov_key"]],
)

# Build scVI model
model_scvi = scvi.model.SCVI(adata=adata_hvg)

model_scvi.view_anndata_setup()

# Train scVI model
model_scvi.train(accelerator="cuda", enable_progress_bar=False)

# Write scVI model to disk
model_scvi.save(dir_path=Path(out_dir_path, "scvi_model"), overwrite=True)

# Load scVI model from disk
model_scvi = scvi.model.SCVI.load(
    dir_path=Path(out_dir_path, "scvi_model"), adata=adata_hvg
)

# Pass the scVI-corrected embedding to the full-length adata
adata.obsm["X_scvi"] = model_scvi.get_latent_representation()

# Compute neighbors using the scVI-corrected embedding
sc.pp.neighbors(adata, use_rep="X_scvi", key_added="neig_scvi")

# Compute UMAP
sc.tl.umap(
    adata=adata,
    neighbors_key="neig_scvi",
    key_added="umap_scvi",
    random_state=random_state,
)

sc.pl.embedding(
    adata=adata,
    basis="umap_scvi",
    neighbors_key="neig_scvi",
    color=[
        "line_name",
        "donor_sex",
        "icell8_chip_id",
        "libsize",
        "pct_mt_genes",
        "pct_rb_genes",
    ],
    title=[
        "line ID",
        "donor sex",
        "iCELL8 chip ID",
        "total counts",
        "% mitochondrial counts",
        "% ribosomal counts",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=3,
    legend_fontsize="medium",
)

sc.pl.embedding(
    adata=adata,
    basis="umap_scvi",
    neighbors_key="neig_scvi",
    color=[
        "line_name",
        "genotype",
        "donor_sex",
        "icell8_chip_id",
    ],
    title=[
        "line ID",
        "genotype",
        "donor sex",
        "iCELL8 chip ID",
    ],
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    legend_fontsize="medium",
)

# Find clusters
leiden_key = "leid_scvi_"
for i in np.round(np.arange(0.1, 1.1, 0.1), 1):
    sc.tl.leiden(
        adata=adata,
        neighbors_key="neig_scvi",
        resolution=i,
        key_added=f"{leiden_key}{str(i)}",
        random_state=random_state,
    )

# sc.pl.embedding(
#     adata=adata,
#     basis="umap_scvi",
#     neighbors_key="neig_scvi",
#     color=[
#         "line_name",
#         f"{leiden_key}0.6",
#         "SNAP25",
#         "top25pct_SNAP25_expressing",
#         "TH",
#         "top25pct_TH_expressing",
#         "AQP4",
#         "top25pct_AQP4_expressing",
#     ],
#     title=[
#         "line ID",
#         "leiden scVI 0.6",
#         "SNAP25",
#         "top 25% SNAP25-expressing",
#         "TH",
#         "top 25% TH-expressing",
#         "AQP4",
#         "top 25% AQP4-expressing",
#     ],
#     layer="tpm",
#     sort_order=True,
#     frameon=False,
#     edges=False,
#     cmap="Reds",
#     ncols=4,
#     gene_symbols="gene_symbol",
#     vmin=0,
#     vmax="p99",
#     legend_fontsize="medium",
# )

# Define neuronal cell type
cell_type = "cell_type_scvi"
reso = 0.6
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["5"]),
    "astrocyte",
    "neuron",
)
adata.obs[cell_type] = np.where(
    adata.obs[f"{leiden_key}{reso}"].isin(["6"]), "other", adata.obs[cell_type]
)
adata.obs[cell_type] = adata.obs[cell_type].astype("category")
adata.obs[cell_type] = adata.obs[cell_type].cat.reorder_categories(
    ["neuron", "astrocyte", "other"]
)

sc.pl.embedding(
    adata=adata,
    basis="umap_scvi",
    neighbors_key="neig_scvi",
    color=[
        f"{leiden_key}{reso}",
        cell_type,
        "SNAP25",
        "top25pct_SNAP25_expressing",
        "TH",
        "top25pct_TH_expressing",
        "AQP4",
        "top25pct_AQP4_expressing",
    ],
    title=[
        f"leiden scVI {reso}",
        "cell type scVI",
        "SNAP25",
        "top 25% SNAP25-expressing",
        "TH",
        "top 25% TH-expressing",
        "AQP4",
        "top 25% AQP4-expressing",
    ],
    layer="tpm",
    sort_order=True,
    frameon=False,
    edges=False,
    cmap="Reds",
    ncols=4,
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    legend_fontsize="medium",
)

# %% Choose cell type annotation for next analyses ----

# Define definitive integration method based on user's parameters
integration_method = current_params["integration_method"]

# Choose cell type annotation from chosen integration method
adata.obs["cell_type"] = adata.obs[f"cell_type_{integration_method}"].copy()

# %% Write adata to disk ----

# Write annotated adata to disk
adata.write_h5ad(filename=Path(out_dir_path, f"adata_{output_param_log_id}.h5ad"))

# Write cell type annotation to disk
cell_type_annotation = adata.obs.loc[:, ["barcode", "cell_type"]].copy()

cell_type_annotation.to_csv(
    Path(out_dir_path, f"cell_type_annotation_{output_param_log_id}.csv"), index=False
)

# %% Pretty figures ----

# Unintegrated
ax = sc.pl.embedding(
    adata=adata,
    basis="umap_unint",
    neighbors_key="neig_unint",
    color=["line_name_generic"],
    frameon=False,
    legend_fontsize=14,
    legend_fontoutline=2,
    show=False,
)

ax.set_title(
    label="unintegrated", fontdict={"fontsize": 14, "fontweight": "bold"}, loc="center"
)
plt.savefig(
    Path(out_dir_path, "line_name_unintegrated.png"), bbox_inches="tight", pad_inches=0
)

# BBKNN integration
ax = sc.pl.embedding(
    adata=adata,
    basis="umap_bbknn",
    neighbors_key="neig_bbknn",
    color=["line_name_generic"],
    frameon=False,
    legend_fontsize=14,
    legend_fontoutline=2,
    show=False,
)

ax.set_title(
    label="BBKNN integration",
    fontdict={"fontsize": 14, "fontweight": "bold"},
    loc="center",
)
plt.savefig(
    Path(out_dir_path, "line_name_bbknn.png"), bbox_inches="tight", pad_inches=0
)

# Harmony integration
ax = sc.pl.embedding(
    adata=adata,
    basis="umap_harmony",
    neighbors_key="neig_harmony",
    color=["line_name_generic"],
    frameon=False,
    legend_fontsize=14,
    legend_fontoutline=2,
    show=False,
)

ax.set_title(
    label="Harmony integration",
    fontdict={"fontsize": 14, "fontweight": "bold"},
    loc="center",
)
plt.savefig(
    Path(out_dir_path, "line_name_harmony.png"), bbox_inches="tight", pad_inches=0
)

# scVI integration
ax = sc.pl.embedding(
    adata=adata,
    basis="umap_scvi",
    neighbors_key="neig_scvi",
    color=["line_name_generic"],
    frameon=False,
    legend_fontsize=14,
    legend_fontoutline=2,
    show=False,
)

ax.set_title(
    label="scVI integration",
    fontdict={"fontsize": 14, "fontweight": "bold"},
    loc="center",
)
plt.savefig(Path(out_dir_path, "line_name_scvi.png"), bbox_inches="tight", pad_inches=0)

# Chosen integration method
reso = 0.6
ax = sc.pl.embedding(
    adata=adata,
    basis=f"umap_{integration_method}",
    neighbors_key=f"neig_{integration_method}",
    color=[
        "line_name_generic",
        "genotype",
        "donor_sex",
        f"leid_{integration_method}_{reso}",
        cell_type,
        "SNAP25",
        "TH",
        "AQP4",
        "COL1A1",
    ],
    layer="tpm",
    sort_order=True,
    cmap="Reds",
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    frameon=False,
    ncols=3,
    legend_fontsize=14,
    legend_fontoutline=2,
    wspace=0.5,
    show=False,
)

for i, j in zip(
    ax,
    [
        "line ID",
        "genotype",
        "donor sex",
        f"leiden {reso}",
        "cell type",
        "SNAP25",
        "TH",
        "AQP4",
        "COL1A1",
    ],
):
    i.set_title(label=j, fontdict={"fontsize": 14, "fontweight": "bold"}, loc="center")
plt.savefig(
    Path(out_dir_path, f"umap_{integration_method}.png"),
    bbox_inches="tight",
    pad_inches=0,
)

ax = sc.pl.embedding(
    adata=adata,
    basis=f"umap_{integration_method}",
    neighbors_key=f"neig_{integration_method}",
    color=[
        cell_type,
        "SNAP25",
        "TH",
        "AQP4",
        "COL1A1",
    ],
    layer="tpm",
    sort_order=True,
    cmap="Reds",
    gene_symbols="gene_symbol",
    vmin=0,
    vmax="p99",
    frameon=False,
    ncols=2,
    legend_fontsize=14,
    legend_fontoutline=2,
    wspace=0.5,
    show=False,
)

for i, j in zip(
    ax,
    [
        "cell type",
        "SNAP25",
        "TH",
        "AQP4",
        "COL1A1",
    ],
):
    i.set_title(label=j, fontdict={"fontsize": 14, "fontweight": "bold"}, loc="center")
plt.savefig(
    Path(out_dir_path, f"umap_{integration_method}_slim.png"),
    bbox_inches="tight",
    pad_inches=0,
)

# %%
