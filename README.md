# ***Pathological levels of α-synuclein increase basal mitophagy in Parkinson’s dopamine neurons***

**This repository contains all the scripts and metadata necessary to replicate the transcriptomics analyses included in the manuscript.**


---


## Data / code availability

- Raw sequencing data (fastq files): available from the Sequence Read Archive under the BioProject accession number XXX. 
- Processed data (Salmon’s quant.sf files): accessible from the Gene Expression Omnibus through GEO Series accession number XXX. 
- Scripts / metadata / helper functions: included in this GitHub repository ('benvallin/Mock_2026'), within folders 'scripts', 'metadata' and 'python_helper_functions', respectively.


---


## Code execution steps

1. make_cell_metadata/make_cell_metadata.R (output provided in 'metadata/cell_metadata.csv')
2. feature_metadata/make_feature_metadata.R
3. count_processing/summarise_to_gene_level.R
4. update_cell_metadata/add_count_metrics.R
5. update_cell_metadata/mark_barcode_quality.R
6. cell_clustering/build_adata.py
7. cell_clustering/integration.py
8. update_cell_metadata/add_cell_type_annotation.R
9. cell_clustering/summarise_cell_counts.R
10. dge_gsea/dge_mast.R
11. dge_gsea/gsea_dge_mast.R
12. summarise_gsea_dge_mast_results.R


---


## Requirements

- R package 'rnaseqtools' available from GitHub with 'devtools::install_github("benvallin/rnaseqtools")'
- R package 'paramstream' available from GitHub with 'devtools::install_github("benvallin/paramstream")'
- Python package 'paramstream' available from PyPI with 'pip install py-paramstream'


---