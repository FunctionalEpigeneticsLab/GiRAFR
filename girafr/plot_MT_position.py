#!/usr/bin/env python

import re

structure_gtf = open('oligo_pool_plasmid_structure.gtf','r')

structure = {}
for s in structure_gtf:
	s = s.strip()
	chrom = s.split()[0]
	gene_name = re.search('\"(.+)\"', s).group(1)
	if gene_name in structure.keys():
		if s.split()[2] in structure[gene_name].keys():
			print("Duplicated structure annotatiton")
	if gene_name not in structure.keys():
		structure[gene_name] = {s.split()[2]: [chrom, int(s.split()[3]), int(s.split()[4])]}
	else:
		structure[gene_name][s.split()[2]] = [chrom, int(s.split()[3]), int(s.split()[4])]
structure_gtf.close()




MT = open('cells.gRNA.single.MT.txt','r')
#AAAGGATAGGACTGGT	1	MED23_gRNA2_gene_WT	17	None
out = open('MT.txt','w')

for line in MT:
	line = line.strip()
	gRNA = line.split()[2]
	mutation = line.split()[4]
	if mutation == 'None':
		continue
	if bool(re.search('_WT',gRNA)):
		continue
	gRNA_gene = re.split('_variant', gRNA)[0]
	mutations = mutation.split(',')
	for m in mutations:
		s =  m.split(':')[0]
		relative_pos = int(re.findall('(\d+)(.+)',m.split(':')[1])[0][0]) + 1
		mutation = re.findall('(\d+)(.+)',m.split(':')[1])[0][1]	
		length = int(structure[gRNA_gene][s][2]) - int(structure[gRNA_gene][s][1]) + 1
		out_line = gRNA_gene + '\t' + s + '\t'+ str(length) + '\t'+ str(relative_pos) + '\t' + mutation + '\t' + gRNA
		out.write(out_line + '\n')
		
	
	


