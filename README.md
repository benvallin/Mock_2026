# ***Pathological levels of α-synuclein increase basal mitophagy in Parkinson’s dopamine neurons***

**This repository contains the scripts and metadata necessary to replicate the transcriptomics analyses included in the manuscript.**


---


## Data / code availability

- Raw (FASTQ) and processed (Salmon quant.sf and meta_info.json) scRNA-seq data: available from the Gene Expression Omnibus (GEO) under accession number GSE354120.
- Metadata: included in this GitHub repository ('benvallin/Mock_2026'), within folder 'metadata'.
- Scripts / helper functions: included in this GitHub repository ('benvallin/Mock_2026'), within folders 'scripts' and 'python_helper_functions', respectively.
- Software details: included in this GitHub repository ('benvallin/Mock_2026'), within folder 'software_details'.


---


## Code execution steps

1. feature_metadata/make_feature_metadata.R
2. count_processing/summarise_to_gene_level.R
3. update_cell_metadata/add_count_metrics.R
4. update_cell_metadata/mark_barcode_quality.R
5. cell_clustering/build_adata.py
6. cell_clustering/integration.py
7. update_cell_metadata/add_cell_type_annotation.R
8. cell_clustering/summarise_cell_counts.R
9. dge_gsea/dge_mast.R
10. dge_gsea/gsea_dge_mast.R
11. summarise_gsea_dge_mast_results.R


---


## Requirements

- R package `rnaseqtools` available from [GitHub](https://github.com/benvallin/rnaseqtools) with:

```r
remotes::install_github("benvallin/rnaseqtools")
```

- R package `paramstream` available from [GitHub](https://github.com/benvallin/paramstream) with:

```r
remotes::install_github("benvallin/paramstream")
```

- Python package `paramstream` available from [PyPI](https://pypi.org/project/py-paramstream/) with:

```bash
pip install py-paramstream
```


---