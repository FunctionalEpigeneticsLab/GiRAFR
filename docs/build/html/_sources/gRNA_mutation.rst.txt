Single cell CRISPR screen gRNA library mutation profiling
=========================================================

.. _configuration_gRNA:

Configuration
-------------

Prepare **ConfigFile** for girafr gRNA_mutation. 

see example 

``gRNA_bam_file``
        gRNA library alignment bam file, mapped to custom reference with gRNA as artificial chromosome (starting ‘_chrom’) 
``filtered_barcode``
        filtered barcode list of expression library mapping results (unzipped). 
        Cellranger corrected cell barcodes (CB tag) is saved as filtered barcodes in the cellranger output. Two alignments can cause cellranger correct same cell barcodes in different ways and cause less cell pass bam filter. 
``min_reads``
        Integer number. By default = 1. UMIs with less than min_reads will be filtered out.
``auto``
        Boolean. By default is True. 
``pool``
        Boolean. By default is True


Requirements
------------
* samtools
* 2bit genome downloaded from UCSC

All coordinates need to be 1-based

gRNA mutation profile
---------------------

.. code-block:: bash

        girafr gRNA_mutation -f path/ConfigFile

Simplified process
------------------

* Step1 (utils.gRNA_bam_filter): 
        gRNA bam file filtration. 
        Script will remove secondary alignments and those which are not aligned to designed gRNA cassette. 
* Step2 (consensus_sequence.generate_consensus_sequence_gRNA):
        Next, script will construct consensus sequence for each UMI. We take the most reads supported sequence as the consensus sequence for the UMI. 
        Arg:
                bam_in: bam file, alignment file of gRNA library after removed secondary alignment and mapped not on gRNA reference
                barcodes: filtered barcode list
        Return:
                consensus.sequence.gRNA.txt: consensus sequence supported by column 7 (n_consensus_reads) > 1
                consensus.seqeunce.gRNA.all_umi.txt: all UMI detected consensus sequence
                consensus.bam: consensus sequence in bam file format
                Non- consensus.bam: alignment not the same as consensus sequence
* Step3 (variant.call_gRNA_variant): call gRNA mutation from consensus.bam. 
        Then, we compare consensus sequence of each UMI with its reference and annotate where the mutation is by structure annotation. Variances are encoded in a similar way like CIGAR in sam format. 
        Arg:
                consensus_seq_file: consensus.sequence.gRNA.txt
                ref_fasta: oligo_pool_plasmid.fa
                structure_gtf: oligo_pool_plasmid_structure.gtf
        Return:
                consensus.sequence.gRNA.variant.txt
* Step4 (assign_gRNA.assign_gRNA_to_cell): Assign gRNA to cell. 
        In the end, we assign found guides to cells. It is required for a cell to have more than min_umi molecule of gRNA so that the script will assign the gRNA to that cell. This umi threshold can be defined by min_umi as fixed threshold for all guides or it can be automatically calculated by fitting a two model mixed gaussian model when auto is set as true in the configuration file. 
        Args:
                in_file: consensus.sequence.gRNA.variant.txt
                min_umi: minum number of UMI, default is 3
                auto: boolean, whether use autodetection or fixed min_umi, default is false
                pool: boolean, whether calculate min umi thresholds together with variant gRNA of the same guide, default is false
        Return:
                Write cells.gRNA.txt and cells.gRNA.single.txt
                consensus.sequence.matrix
                gRNA.umi.threshold.txt
* Step5 (optional): (assign_gRNA.add_variant_type, profile_MT_pattern.py) add mutation details for downstream analysis.
        Args:
                in_file1: consensus.sequence.gRNA.variant.txt
                in_file2: cells.gRNA.single.txt
                in_file3: cells.gRNA.txt
        Return:
                Write: cells.gRNA.single.MT.txt, MT.txt and all.MT.txt

        profile_MT_pattern.py
        Args:
                in_file: consensus.sequence.gRNA.variant.txt
                min_umi2: ninum number of umi for certain Variant to be profiled
        Return:
                Write consensus.sequence.gRNA.MT.txt
 


Output files
------------
See section

Additional information
----------------------
Build custom reference (optional):
Oligo_pool.csv: two columns: oligo_name and sequence, no header.
prepare.py: generate oligo_pool_plasmid.fa, oligo_pool_plasmid.gtf and oligo_pool_plasmid_structure.gtf 
This part gives instruction to build a custom CellRanger reference with designed cassette as artificial chromosome. utils.write_annotation function generates oligo_pool_plasmid.fa and oligo_pool_plasmid.gtf which will be used to generate cellranger reference (see build note as example), and oligo_pool_plasmid_structure.gtf which will be used to profile where the mutations are on the cassatte.

CIGAR-like string:
        Digit numbers represents exact matches, and nucleotides followed are mutated bases. 0 represents no nucleotide. 
        Digit numbers followed by insertions (I), deletions (D) and soft clippings (S) show the number of nucleotides of those events. Hard clippings (H) are not included. The major difference between this string and CIGAR-string is it replaces matches (M) into mismatches and encode detailed mutated nucleotides [ATGC] into the string. 
Mutation structure annotation:
        Annotations begin with oligo structures such as gRNA which are consistent with user input oligo_pool_plasmid_structure.gtf. Then each mutation annotation follows oligo structure with semicolon as separator. Comma separates individual mutation event. Digit numbers represents the distance to the beginning of the structure. Nucleotides followed are mutated bases. 0 represents no nucleotide. Digit numbers in bracket followed insertions (I), deletions (D) and soft clippings (S) represent the number of nucleotides of those events. 


