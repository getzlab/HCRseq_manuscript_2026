# HCRseq_manuscript_2026

This repository contains code to reproduce the results in Chen et al. 

Scripts rely on installation of the [HCRseq package](https://github.com/getzlab/HCRseq). Preprocessing of the HCR-seq data was performed in Terra using the workflows in wdl/. Normalization, clustering, and annotation of of the single-cell data was performed using the code in preprocess_scrna.py. Analysis of the nanopore data was performed using mmHCRseqNanopore_v2.py. Downstream analyses and figure generation were performed using the notebooks in figures/. 

The raw and processed data can be downloaded from GEO at accessions GSE330134 (amplicon data) and GSE330132 (scRNA data).

