# GiRAFR (Guide RNA Anomalies and Functionality Revealed)


# Requirements:
* umi_tools 
* samtools
* Pysam
* 2bit genome available from UCSC

All coordinate should be 1-based

# gRNA library mutation profiling 
## Input (ConfigFile)
gRNA_bam_file: gRNA library alignment bam file, mapped to custom reference with gRNA as artificial chromosome (starting ‘_chrom’)  
filtered_barcode: filtered barcode list of expression library mapping results (unzipped).  
min_reads=1 (by default)    
auto=True (by default)   
pool=True (by default)   

# Output 
consensus.bam: Bam file contains consensus sequence alignment result for each UMI.   
consensus.sequence.gRNA.variant.txt:   
Column 1: Cell barcode.  
Column 2: UMI.  
Column 3: Number of UMI in the cell.   
Column 4: Number of reads for this UMI (default min = 1).  
Column 5: single/multiple. Single refers to UMI with only one type of sequence. Multiple refers to UMI with more than one type of sequence. ‘Multiple’ will be skipped in the later steps. (take up to xx%) 256794 multiple 153598 singles.   
Column 6: Consensus sequence. Sequence with most reads supported.   
Column 7: Number of reads that support the consensus sequence.   
Column 8: True/False. True: wild type sequence. False: different from reference.   
Column 9: Name of the guide. ‘multiple gene’: when UMI mapped to multiple guides. Those UMI will be skipped later. (15057 UMI out of 410393 mapped to multiple guide)   
Column 10: Guide variant type with detailed mismatch information.  
Column 11: Mutation structure annotation  

cells.gRNA.txt   
Column 1: Cell barcode. 
Column 2: NO. guide in the cell. 
Column 3: Guide name.    
Column 4: WT: Wild type. 
Column 5: Number of UMI support wild type guide
Column 6: Guide variant type if exist. 
Column 7: Number of UMI support mutant guide. 
cells.gRNA.single.txt: Cells with only one guide
Column 1: Cell barcode
Column 2: 1 as only one guide in the cell
Column 3: Guide name. Guide name ends with WT represents wild type. 
Column 4: Number of UMI support mutant guide.  

cells.gRNA.single.MT.txt   
Column 1 – Column 4: same as cells.gRNA.single.txt  
Column 5: Mutation details. None when wildtype or soft clipping only.  
Column 6: Mutant sequence. None when wildtype or soft clipping only.  

MT.txt: mutations in single cell in a long format, each line contains single mutation event.   
Column 1: target gene  
Column 2: structure (defined in structure.gtf)  
Column 3: length of the structure  
Column 4: Position on the structure    
Column 5: mutation (CIGAR similar format)  
Column 6: variant name with number suffix.   

all.MT.txt: mutations of all in-cell guides same format as MT.txt. An important difference is multiple cells with the same variant guide will be recorded once in this file, but multiple times in MT.txt.   
consensus.sequence.gRNA.MT.txt: same format as all.MT.txt but for all mutations in consensus.sequence.gRNA.variant.txt including guides that are not detected in cells. Variants with at least 2 UMIs will be recorded. Each variant only appears once in this file.  




