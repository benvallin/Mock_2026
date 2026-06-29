# %% Set up ----

# Import required libraries
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# %% set_pretty_cat() ----


# Function to set obs as category and supply colors for plotting
def set_pretty_cat(adata, obsnm, coldict):
    adata.obs[obsnm] = adata.obs[obsnm].copy().astype("category")
    adata.obs[obsnm] = adata.obs[obsnm].cat.reorder_categories(coldict.keys())
    adata.uns[obsnm + "_colors"] = list(coldict.values())


# %% find_top_expressing_cells() ----


# Function to compute the 'top X% marker gene-expressing cells'
def find_top_expressing_cells(adata, layer: None, gene_symbol: str, top_pct: int):

    if gene_symbol not in list(adata.var["gene_symbol"]):
        return print(f"{gene_symbol} not in adata.var['gene_symbol']")

    top_quantile = 1 - (top_pct / 100)

    if layer is None:
        gene_counts = adata[:, adata.var["gene_symbol"] == gene_symbol].X.copy()
    else:
        if layer not in adata.layers.keys():
            return print(f"{layer} not in adata.layers")
        gene_counts = (
            adata[:, adata.var["gene_symbol"] == gene_symbol].layers[layer].copy()
        )

    if isinstance(gene_counts, sparse.csr_matrix):
        gene_counts = gene_counts.todense()
        gene_counts = [x[0] for x in gene_counts.tolist()]
    else:
        gene_counts = gene_counts[:, 0].copy()

    top_quantile_value = np.quantile(a=gene_counts, q=top_quantile, axis=0)

    top_pct_expressing = np.where(gene_counts >= top_quantile_value, True, False)

    top_pct_expressing = pd.Series(
        data=top_pct_expressing, index=adata.obs_names, dtype="category"
    )

    if len(set(top_pct_expressing)) > 1:
        top_pct_expressing = top_pct_expressing.cat.reorder_categories(
            new_categories=[True, False]
        )
        adata.obs["top" + str(top_pct) + "pct_" + gene_symbol + "_expressing"] = (
            top_pct_expressing
        )
        adata.uns[
            "top" + str(top_pct) + "pct_" + gene_symbol + "_expressing_colors"
        ] = ["#4361ee", "#e7ece5"]
    else:
        adata.obs["top" + str(top_pct) + "pct_" + gene_symbol + "_expressing"] = (
            top_pct_expressing
        )
        adata.uns[
            "top" + str(top_pct) + "pct_" + gene_symbol + "_expressing_colors"
        ] = ["#907ad6"]


# %% gene_metadata_to_dc_net() ----


def gene_metadata_to_dc_net(gene_metadata: pd.DataFrame, collection: str):

    tmp = gene_metadata.set_index("gene_symbol").loc[:, collection].explode().dropna()

    output = (
        pd.json_normalize(tmp)
        .set_index(tmp.index)
        .reset_index(drop=False)
        .loc[:, ["set_name", "gene_symbol"]]
        .rename(columns={"set_name": "source", "gene_symbol": "target"})
        .drop_duplicates()
    )

    return output


# %% run_score_by_cat_var_lmm() ----


def run_scores_by_cat_var_lmm(
    adata,
    scores_df,
    meta_vars,
    formula_rhs,
    cat_var,
    cat_var_ref_lvl,
    cat_var_test_lvl,
    batch_var,
    score_var_nm="score_var_nm",
):

    formula = f"score {formula_rhs}"

    coef_name = (
        f"C({cat_var}, Treatment(reference='{cat_var_ref_lvl}'))[T.{cat_var_test_lvl}]"
    )

    scores_meta_df = pd.merge(
        left=adata.obs.loc[:, meta_vars].copy(),
        right=scores_df,
        left_index=True,
        right_index=True,
    )

    score_names = set(scores_meta_df.columns) - set(meta_vars)

    output = []

    for i in tqdm(score_names):
        df = scores_meta_df.loc[:, [i] + meta_vars].copy().rename(columns={i: "score"})

        try:
            lmm = smf.mixedlm(formula=formula, data=df, groups=batch_var).fit(reml=True)

            output.append(
                {
                    score_var_nm: i,
                    "coef": lmm.params[coef_name],
                    "se": lmm.bse[coef_name],
                    f"mean_{cat_var_ref_lvl}": df.loc[
                        df[cat_var] == cat_var_ref_lvl, "score"
                    ].mean(),
                    f"mean_{cat_var_test_lvl}": df.loc[
                        df[cat_var] == cat_var_test_lvl, "score"
                    ].mean(),
                    "pval": lmm.pvalues[coef_name],
                    "converged": lmm.converged,
                }
            )

        except Exception as e:
            output.append(
                {
                    score_var_nm: i,
                    "coef": np.nan,
                    "se": np.nan,
                    f"mean_{cat_var_ref_lvl}": np.nan,
                    f"mean_{cat_var_test_lvl}": np.nan,
                    "pval": np.nan,
                    "converged": False,
                }
            )

    output = pd.DataFrame(output)

    # Flag non-converged models
    n_failed = (~output["converged"]).sum()
    print(f"{n_failed} / {len(output)} models failed to converge")

    # Apply BH correction only to converged models with valid p-values
    valid = output["converged"] & output["pval"].notna()
    padj = np.full(shape=len(output), fill_value=np.nan)
    padj[valid] = multipletests(pvals=output.loc[valid, "pval"], method="fdr_bh")[1]
    output["padj"] = padj

    # Sort by padj and coefficient
    output = output.sort_values(
        by=["padj", "coef"], ascending=[True, False]
    ).reset_index(drop=True)

    return output


# %% run_score_by_scores_lmm() ----


def run_score_by_scores_lmm(
    adata,
    dep_score_df,
    exp_scores_df,
    batch_var,
    dep_score_var_nm="dep_score_var_nm",
    exp_score_var_nm="exp_score_var_nm",
):

    scores_meta_df = pd.merge(
        left=adata.obs.loc[:, batch_var].copy(),
        right=dep_score_df,
        left_index=True,
        right_index=True,
    ).merge(
        right=exp_scores_df,
        left_index=True,
        right_index=True,
    )

    output = []

    dep_score_name = dep_score_df.columns[0]

    for i in tqdm(exp_scores_df.columns):

        df = (
            scores_meta_df.loc[:, [dep_score_name, i, batch_var]]
            .copy()
            .rename(columns={dep_score_name: "dep_score", i: "exp_score"})
        )

        try:
            lmm = smf.mixedlm(
                formula="dep_score ~ exp_score", data=df, groups=batch_var
            ).fit(reml=True)

            output.append(
                {
                    dep_score_var_nm: dep_score_name,
                    exp_score_var_nm: i,
                    "coef": lmm.params["exp_score"],
                    "se": lmm.bse["exp_score"],
                    "pval": lmm.pvalues["exp_score"],
                    "converged": lmm.converged,
                }
            )

        except Exception as e:
            output.append(
                {
                    dep_score_var_nm: dep_score_name,
                    exp_score_var_nm: i,
                    "coef": np.nan,
                    "se": np.nan,
                    "pval": np.nan,
                    "converged": False,
                }
            )

    output = pd.DataFrame(output)

    # Flag non-converged models
    n_failed = (~output["converged"]).sum()
    print(f"{n_failed} / {len(output)} models failed to converge")

    # Apply BH correction only to converged models with valid p-values
    valid = output["converged"] & output["pval"].notna()
    padj = np.full(shape=len(output), fill_value=np.nan)
    padj[valid] = multipletests(pvals=output.loc[valid, "pval"], method="fdr_bh")[1]
    output["padj"] = padj

    # Sort by padj and coefficient
    output = output.sort_values(
        by=["padj", "coef"], ascending=[True, False]
    ).reset_index(drop=True)

    return output
