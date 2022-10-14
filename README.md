# GiRAFR (Guide RNA Aberrations and Functionality Revealed)

All coordinate should be 1-based

# gRNA library mutation profiling 

# Requirements:
* umi_tools
* samtools
* Pysam
* 2bit genome available from UCSC

## Input (ConfigFile)
*Example configfile can be found with the scripts*
gRNA_bam_file: gRNA library alignment bam file, mapped to custom reference with gRNA as artificial chromosome (starting '_chrom')  
filtered_barcode: filtered barcode list of expression library mapping results (unzipped).  
min_reads=1 (by default)    
auto=True (by default)   
pool=True (by default)   

# Output 
**consensus.bam**  
Bam file contains consensus sequence alignment result for each UMI.  
 
**consensus.sequence.gRNA.variant.txt**  
Column 1: Cell barcode.  
Column 2: UMI.  
Column 3: Number of UMI in the cell.   
Column 4: Number of reads for this UMI (default min = 1).  
Column 5: single/multiple. Single refers to UMI with only one type of sequence. Multiple refers to UMI with more than one type of sequence. 
Column 6: Consensus sequence. Sequence with the most reads supported.   
Column 7: Number of reads that support the consensus sequence.   
Column 8: True/False. True: wild type sequence. False: different from reference.   
Column 9: Name of the guide. 'multiple gene': when UMI mapped to multiple guides with equal number of reads. The progam cannot assign gene to UMI confidently. In this case, UMI will be skipped later.    
Column 10: Guide variant type with detailed mismatch information.  
Column 11: Mutation structure annotation.  

**cells.gRNA.txt**   
Column 1: Cell barcode. 
Column 2: Number of guides in the cell. 
Column 3: Guide name.    
Column 4: WT: Intact gRNA
Column 5: Number of UMI support wild type guide
Column 6: Guide variant type if exist. 
Column 7: Number of UMI support mutant guide. 

**cells.gRNA.single.txt**  
Cells with only one guide
Column 1: Cell barcode
Column 2: 1 as only one guide in the cell
Column 3: Guide name. Guide name ends with WT represents wild type. 
Column 4: Number of UMI support mutant guide.  

**cells.gRNA.single.MT.txt**   
Column 1 - Column 4: same as cells.gRNA.single.txt  
Column 5: Mutation details. None when wildtype or soft clipping only.  
Column 6: Mutant sequence. None when wildtype or soft clipping only.  

**MT.txt**  
mutations in single cell in a long format, each line contains single mutation event.   
Column 1: target gene  
Column 2: structure (defined in structure.gtf)  
Column 3: length of the structure  
Column 4: Position on the structure    
Column 5: mutation (CIGAR similar format)  
Column 6: variant name with number suffix.   

**all.MT.txt**  
mutations of all in-cell guides same format as MT.txt. An important difference is multiple cells with the same variant guide will be recorded once in this file, but multiple times in MT.txt.   
**consensus.sequence.gRNA.MT.txt**  
same format as all.MT.txt but for all mutations in consensus.sequence.gRNA.variant.txt including guides that are not detected in cells. Variants with at least 2 UMIs will be recorded. Each variant only appears once in this file.  

# Simplified process
### Step1 (utils.gRNA_bam_filter)  
gRNA bam file filtration. 
Script will remove secondary alignments and those which are not aligned to designed gRNA cassette.   
### Step2 (consensus_sequence.generate_consensus_sequence_gRNA)  
Next, script will construct consensus sequence for each UMI. We take the most reads supported sequence as the consensus sequence for the UMI.  
Arg:  
		bam_in: bam file, alignment file of gRNA library after removed secondary alignment and mapped not on gRNA reference
		barcodes: filtered barcode list
Return:  
		consensus.sequence.gRNA.txt: consensus sequence supported by column 7 (n_consensus_reads) > 1
		consensus.seqeunce.gRNA.all_umi.txt: all UMI detected consensus sequence
		consensus.bam: consensus sequence in bam file format
		Non- consensus.bam: alignment not the same as consensus sequence
### Step3 (variant.call_gRNA_variant)  
call gRNA mutation from consensus.bam. 
Then, we compare consensus sequence of each UMI with its reference and annotate where the mutation is by structure annotation. Variances are encoded in a similar way like CIGAR in sam format. 
	Arg:
		consensus_seq_file: consensus.sequence.gRNA.txt
		ref_fasta: oligo_pool_plasmid.fa
		structure_gtf: oligo_pool_plasmid_structure.gtf
	Return:
		consensus.sequence.gRNA.variant.txt: 
### Step4 (assign_gRNA.assign_gRNA_to_cell)  
Assign gRNA to cell. 
In the end, we assign found guides to cells. It is required for a cell to have more than min_umi molecule of gRNA so that the script will assign the gRNA to that cell. This umi threshold can be defined by min_umi as fixed threshold for all guides or it can be automatically calculated by fitting a two model mixed gaussian model when auto is set as true in the configuration file. 
Args:
in_file: consensus.sequence.gRNA.variant.txt
min_umi: minum number of UMI, default is 3
auto: boolean, whether use autodetection or fixed min_umi, default is false
pool: boolean, whether calculate min umi thresholds together with variant gRNA of the same guide, default is false
Returns:
Write cells.gRNA.txt and cells.gRNA.single.txt
consensus.sequence.matrix
gRNA.umi.threshold.txt

### Step5 (optional) (assign_gRNA.add_variant_type, profile_MT_pattern.py)  
add mutation details for downstream analysis.
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
Returns:
Write consensus.sequence.gRNA.MT.txt

## Build custom reference (optional)
### Input
Oligo_pool.csv: two columns: oligo_name and sequence, no header. 
This part gives instruction to build a custom CellRanger reference with designed cassette as artificial chromosome. utils.write_annotation function generates oligo_pool_plasmid.fa and oligo_pool_plasmid.gtf which will be used to generate cellranger reference (see build note as example), and oligo_pool_plasmid_structure.gtf which will be used to profile where the mutations are on the cassatte.

## Additional Information
CIGAR-like string:
Digit numbers represents exact matches, and nucleotides followed are mutated bases. 0 represents no nucleotide. 
Digit numbers followed by insertions (I), deletions (D) and soft clippings (S) show the number of nucleotides of those events. Hard clippings (H) are not included. The major difference between this string and CIGAR-string is it replaces matches (M) into mismatches and encode detailed mutated nucleotides [ATGC] into the string. 

Mutation structure annotation:
Annotations begin with oligo structures such as gRNA which are consistent with user input oligo_pool_plasmid_structure.gtf. Then each mutation annotation follows oligo structure with semicolon as separator. Comma separates individual mutation event. Digit numbers represents the distance to the beginning of the structure. Nucleotides followed are mutated bases. 0 represents no nucleotide. Digit numbers in bracket followed insertions (I), deletions (D) and soft clippings (S) represent the number of nucleotides of those events. 

Build note:
annotation.csv: gRNA id must end with sgRNA and unique name, protospacer
prepare.py: generate oligo_pool_plasmid.fa, oligo_pool_plasmid.gtf and oligo_pool_plasmid_structure.gtf


# Detect editing effect
## Input
gRNA_region_coordinates_ori.txt
Column 1: chromosome (same as in bam file)
Column 2: start (1-based)
Column 3: end (1-based)
Column 4: Name (unique) eg. Region2, Region3 ...
Column 5: +/- strand sgRNA targets to which strand
Column 6: designed protospacer sequence. By default, PAM is also given following protospacer sequence. 
Column 7: target gene (gene symbol)
By default, cut site is defined as 3nt to 4nt upstream (5') of the Protospacer Adjacent Motif (PAM). 

Detection window size: 31bp by default. Detection window is symmetric around cutsite. 
Barcode list. Filtered cell barcode list. 
cell sgRNA assignment. file that contain cell - sgRNA annotation
For cellRanger feature barcoding analysis, cell_file is protospacer_calls_per_cell.csv or cells.gRNA.txt generated by previous step or file with 
Column 1: cell barcode
Column 3: sgRNA assigned (new row if multiple sgRNA/cell)

2bit genome download: UCSC
Human: hg38
Mouse: mm10

## Output

For each gRNA:
	{region_name}.bam: alignment within detection window
	{region_name}.bam.featureCounts.sorted.bam: FeatureCount assigned reads to gene.
	{region_name}.consensus.sequence.txt: Consensus sequence within detection window for each UMI and cell. 
Column 1: cell barcode
Column 2: UMI
Column 3: Number of UMI for cell
Column 4: Number of sequences for that UMI
Column 5: reference sequence in the detection window
Column 6: union sequence for that UMI in the detection window
Column 7: union number of reads
Column 8: Boolean True for wildtype sequence, false for mutation exist in sequence
Column 9: mapped gene
Column 10: mapped positions on the genome
Column 11: is_umi_seq_different 
Column 12: reference gRNA id (region_name)


region_name}.editing_effect.cutsite.txt: 
all.editing_effect.cutsite.txt: concatenate all {region_name}.editing_effect.cutsite.txt.
Column 13: mutation, ‘No mutation’ given when wildtype sequence. "Deletion at cutsite”, "Insertion/Deletion somewhere else”, "Mismatch” 
Column 14: gRNA assigned to this cell 

# Process
The script first generates coordinates based on input window size, and reference sequence and later they are used by Samtools to split alignment file (bam file) into bam files which contain reads that aligned to window region. FeatureCount from UMI-tools is employed to filter out reads that not align properly to this target gene. Consensus sequences are generated by union reads within detection window with the same UMI and cell barcode. Reference sequence for each region is compared with the consensus sequence to identify mutations. Deletions that pass cutsite are defined as editing effect. 


