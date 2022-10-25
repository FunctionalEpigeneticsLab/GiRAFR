Output file format
==================

.. _output:

gRNA_mutation
-------------

consensus.bam:
Bam file contains consensus sequence alignment result for each UMI. 

consensus.sequence.gRNA.variant.txt

.. list-table::  
   :header-rows: 1
   
   * - column number
     - Value

   * - Column 1
     - Cell barcode
   * - Column 2
     - UMI

   * - Column 3
     - Number of UMI in the cell
   * - Column 4
     - Number of reads for this UMI (default min = 1)
   * - Column 5
     - single/multiple. Single refers to UMI with only one type of sequence. Multiple refers to UMI with more than one type of sequence.
   * - Column 6
     -  Consensus sequence. Sequence with most reads supported. 
   * - Column 7
     -  Number of reads that support the consensus sequence. 
   * - Column 8
     -  True/False. True: wild type sequence. False: different from reference. 
   * - Column 9
     -  Name of the guide. 'multiple gene': when UMI mapped to multiple guides with equal number of reads. The progam cannot assign gene to UMI confidently. In this case, UMI will be skipped later.
   * - Column 10
     -  Guide variant type with detailed mismatch information
   * - Column 11
     -  Mutation structure annotation 


cells.gRNA.txt

.. list-table:: 
   :header-rows: 1
   
   * - column number
     - Value

   * - Column 1
     -  Cell barcode
   * - Column 2
     -  Number of guides in the cell
   * - Column 3
     -  Guide name. 
   * - Column 4
     -  WT: Intact gRNA
   * - Column 5
     -  Number of UMI support wild type guide
   * - Column 6
     -  Guide variant type if exist. 
   * - Column 7
     -  Number of UMI support mutant guide. 

cells.gRNA.single.txt: cells with only one gRNA

.. list-table::
   :header-rows: 1

   * - column number
     - Value

   * - Column 1
     -  Cell barcode
   * - Column 2
     -  1 as only one guide in the cell
   * - Column 3
     -  Guide name. Guide name ends with WT represents wild type. 
   * - Column 4
     -  Number of UMI support mutant guide.  

cells.gRNA.single.MT.txt

.. list-table::
   :header-rows: 1

   * - column number
     - Value

   * - Column 1 to Column 4
     -  same as cells.gRNA.single.txt
   * - Column 5
     -  Mutation details. None when wildtype or soft clipping only.
   * - Column 6
     -  Mutant sequence. None when wildtype or soft clipping only.

all.MT.txt: mutations in single cell in a long format, each line contains single mutation event.

.. list-table::
   :header-rows: 1

   * - column number
     - Value

   * - Column 1
     -  target gene
   * - Column 2
     -  structure (defined in structure.gtf)
   * - Column 3
     -  length of the structure
   * - Column 4
     -  Position on the structure
   * - Column 5
     -  mutation (:ref:`CIGAR similar format <CIGAR_string>`)
   * - Column 6
     -  variant name with number suffix. 

editing_effect
--------------
For each gRNA: 
``{region_name}.bam``: alignment within detection window
``{region_name}.bam.featureCounts.sorted.bam``: FeatureCount assigned reads to gene.
``{region_name}.consensus.sequence.txt``: Consensus sequence within detection window for each UMI and cell. 

.. list-table::
   :header-rows: 1

   * - column number
     - Value

   * - Column 1
     -  cell barcode
   * - Column 2
     -  UMI
   * - Column 3
     -  Number of UMI for cell
   * - Column 4
     -  Number of sequences for that UMI
   * - Column 5
     -  reference sequence in the detection window
   * - Column 6
     -  union sequence for that UMI in the detection window
   * - Column 7
     -  union number of reads
   * - Column 8
     -  Boolean True for wildtype sequence, false for mutation exist in sequence
   * - Column 9
     -  mapped gene
   * - Column 10
     -  mapped positions on the genome
   * - Column 11
     -  is_umi_seq_different 
   * - Column 12
     -  reference gRNA id (region_name)

{region_name}.editing_effect.cutsite.txt

.. list-table::
   :header-rows: 1

   * - column number
     - Value

   * - Column 1 to Column 12
     - same as {region_name}.consensus.sequence.txt
   * - Column 13
     -  mutation, ‘No mutation’ given when wildtype sequence. "Deletion at cutsite”, "Insertion/Deletion somewhere else”, "Mismatch” 
   * - Column 14
     -  gRNA assigned to this cell 


