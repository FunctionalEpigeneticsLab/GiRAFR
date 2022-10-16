#!/usr/bin/bash

results_dir=$1

cd ${results_dir}


cell_nr=`less barcodes.tsv|wc -l`
gRNA_reads=`samtools view -c gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam`
consensus_bam=`samtools view -c consensus.bam`

#num1=`cut -f 5 consensus.sequence.gRNA.txt|sort|uniq -c|tail -n 2`
#num1=`cut -f 9 consensus.sequence.gRNA.txt|grep -v 'gRNA'|grep -v '_'|sort|uniq -c|tail -n 2`

umi_filtered=`less -S consensus.sequence.gRNA.variant.txt|wc -l`

umi_MT=`less -S consensus.sequence.gRNA.variant.txt|cut -f 8|grep 'False'|wc -l` 
umi_WT=`less -S consensus.sequence.gRNA.variant.txt|cut -f 8|grep 'True'|wc -l`

single_cell=`less cells.gRNA.single.txt|wc -l`
cell_withguide=`cut -f 1 cells.gRNA.txt|sort|uniq|wc -l`

echo -e "cell_nr,gRNA_reads,consensus_bam_reads,umi_filtered,umi_MT,umi_WT,single_cell_nr,cell_withguide_nr"
echo -e "$cell_nr,$gRNA_reads,$consensus_bam,$umi_filtered,$umi_MT,$umi_WT,$single_cell,$cell_withguide"
