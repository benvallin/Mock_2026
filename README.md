# ***Pathological levels of α-synuclein increase basal mitophagy in Parkinson’s dopamine neurons***

**This repository contains all the scripts and metadata necessary to replicate the transcriptomics analyses included in the manuscript.**


---


## Data / code availability

- Raw sequencing data (fastq files): available from the Sequence Read Archive under the BioProject accession number XXX. 
- processed data (Salmon’s quant.sf files): accessible from the Gene Expression Omnibus through GEO Series accession number XXX. 
- Scripts / metadata: included in this GitHub repository ('benvallin/Mock_2026'), within folders 'scripts' and 'metadata', respectively.


---


## Code execution steps

- make_cell_metadata/make_cell_metadata.R (output provided in 'metadata/cell_metadata.csv')
- feature_metadata/make_feature_metadata.R
- count_processing/summarise_to_gene_level.R
- update_cell_metadata/add_count_metrics.R
- update_cell_metadata/mark_barcode_quality.R
- cell_clustering/build_adata.py
- cell_clustering/integration.py
- update_cell_metadata/add_cell_type_annotation.R
- cell_clustering/summarise_cell_counts.R
- dge_gsea/dge_mast.R
- dge_gsea/gsea_dge_mast.R
- summarise_gsea_dge_mast_results.R


---