# GiRAFR (Guide RNA Anomalies and Functionality Revealed)

GiRAFR is a python package developed to perform quality control of single cell CRISPR screens, and to assign gRNAs to cells in a sensitive, mutation-aware manner. 

# Introduction/Abstract

Novel single cell RNA-seq analysis combined with CRISPR screens enables the high- throughput characterization of transcriptional changes caused by genetic perturbations. Dedicated software to annotate CRISPR guide RNA (gRNA) libraries and associate them with single cell transcriptomes are however lacking. Here, we generated a CRISPR droplet sequencing dataset. We demonstrate that the current default tool fails to detect mutant gRNAs. We therefore developed GiRAFR, a pysam-based software tool to characterize intact and mutant gRNAs. We show that mutant gRNAs are dysfunctional, and failure to detect and annotate them leads to an inflated estimate of the number of untransformed cells as well as an underestimated multiplet frequency. These findings are mirrored in publicly available datasets, where we find that up to 34 % of cells are transduced with a mutant gRNA. Applying GiRAFR hence stands to improve the annotation and quality of single cell CRISPR screens.

# Installation

``pip install girafr``

# Documentation
* [Single cell CRISPR screen gRNA library mutation profiling](https://girafr.readthedocs.io/en/latest/gRNA_mutation.html)
* [Detect CRISPR-Cas9 editing effects in Single cell RNA-seq](https://girafr.readthedocs.io/en/latest/editing_effect.html)
* [Output files format](https://girafr.readthedocs.io/en/latest/output.html)

# Citation

GiRAFR improves gRNA detection and annotation in single cell CRISPR screens  
Qian Yu, Paulien Van Minsel, Eva Galle, Bernard Thienpont  
bioRxiv 2022.10.24.513352; doi: https://doi.org/10.1101/2022.10.24.513352



# License and copyright


