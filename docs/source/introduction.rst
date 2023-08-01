Introduction
============

GiRAFR is a python package developed to perform quality control of single cell CRISPR screens, and to assign gRNAs to cells in a sensitive, mutation-aware manner.

Citation
--------

GiRAFR improves gRNA detection and annotation in single cell CRISPR screens  
Qian Yu, Paulien Van Minsel, Eva Galle, Bernard Thienpont  
bioRxiv 2022.10.24.513352; doi: https://doi.org/10.1101/2022.10.24.513352  


**Abstract**
Novel single cell RNA-seq analysis combined with CRISPR screens enables the high- throughput characterization of transcriptional changes caused by genetic perturbations. Dedicated software to annotate CRISPR guide RNA (gRNA) libraries and associate them with single cell transcriptomes are however lacking. Here, we generated a CRISPR droplet sequencing dataset. We demonstrate that the current default tool fails to detect mutant gRNAs. We therefore developed GiRAFR, a pysam-based software tool to characterize intact and mutant gRNAs. We show that mutant gRNAs are dysfunctional, and failure to detect and annotate them leads to an inflated estimate of the number of untransformed cells as well as an underestimated multiplet frequency. These findings are mirrored in publicly available datasets, where we find that up to 34 % of cells are transduced with a mutant gRNA. Applying GiRAFR hence stands to improve the annotation and quality of single cell CRISPR screens.

Package workflow
----------------
.. image:: scheme.png
   :height: 400px
   :width: 800px
   :align: left


Quickstart
----------

GiRAFR has two subcommands: ``gRNA_mutation`` and ``editing_effect``. Input file: ``ConfigFile`` is required for both subcommands. 

To create configuration file ``ConfigFile``, see :ref:`create configuration file to detect gRNA mutation <configuration_gRNA>` and :ref:`create configuration file to detect editing effect <configuration_editing>`

* Profile gRNA mutation in Single cell CRIPSR screen
        
.. code-block:: bash

        girafr gRNA_mutation -c abosolute_path/ConfigFile 

GiRAFR takes alignment results (bam file) as its main input. It is best compatible to Cell Ranger alignment output. Therefore, it can be used as downstream analysis after Cell Ranger alignments to provide cell annotations before the transcriptome analysis. GiRAFR is also compatible to Drop-seq toolbox and other alignment results by simply modifying its configuration (tag information).

* Detect CRISPR-Cas9 editting effect in single cell RNA-seq:

.. code-block:: bash
        
        girafr editing_effect -c absolute_path/ConfigFile

If no -c/--config is used, girafr will try ``ConfigFile`` in the current path.


Scalability and Runtime
-----------------------
The runtime of GiRAFR is determined by the sequencing depth and the number of cells in gRNA library. Below figure shows the correlation between GiRAFR runtime and the product of cells number and reads number using down-sampled inhouse CROP-seq data.  

.. image:: scalability.png
   :height: 200px
   :width: 200px
   :align: left

