#!/usr/bin/bash

processed_cells=`head -n 1 consensus.sequence.gRNA.txt | cut -f 2`
echo "Total processed cells: "$processed_cells > stats.txt

#consensus_cells=`cut -f 1 consensus.sequence.gRNA.txt | sort| uniq|wc -l`
#consensus_cells="$(($consensus_cells - 1))"
#echo "Cells consensus sequence with enough reads: "$consensus_cells >> stats.txt

matrix_cells=`head -n 1 consensus.count.matrix|wc -w`
matrix_cells="$(($matrix_cells - 1))"
echo "Cells consensus count matrix: "$matrix_cells >> stats.txt

less -S cells.gRNA.txt |awk '{ x[$1]=$0 }END{for (i in x ) print x[i]}'|awk '{if($2==1) num+=1;else if($2==2) num2+=1;else num3+=1}END{print "Non-empty cells: "num+num2+num3"\nSingle protospacer cells: "num"\n2 protospacers cells: "num2"\nMore than 2 protospacers: "num3}' >> stats.txt


single_WT=`less -S cells.gRNA.single.txt |grep 'WT' |wc -l`
single_MT=`less -S cells.gRNA.single.txt |grep 'variant' |wc -l`
#echo "Single protospacer cells with wildtype: "$single_WT >> stats.txt
echo "Single protospacer cells with mutant gRNA: "$single_MT >> stats.txt

#echo "" >> stats.txt

n_reads=`samtools view -c gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam`
echo "Reads after filter: "$n_reads >> stats.txt
n_UMI=`wc -l < consensus.sequence.gRNA.all_umi.txt`
echo "All UMI: "$n_UMI >> stats.txt
n_UMI=`wc -l < consensus.sequence.gRNA.variant.txt`
echo "UMI with enough reads: "$n_UMI >> stats.txt
less -S consensus.sequence.gRNA.variant.txt|awk '{if($8 == "False") num+=1; else if($8 == "True") num2+=1}END{print "UMI wildtype: "num2"\nUMI with mutations: "num}' >> stats.txt


echo "calculation finsished"
