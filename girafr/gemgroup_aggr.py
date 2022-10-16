#!/usr/bin/env python

import os
from . import utils
import pysam
import subprocess
from . import consensus_sequence
from datetime import datetime
from . import variant
from . import assign_gRNA


print('Combine cells from different gem groups')

#(samtools, twoBitToFa, featureCounts) = utils.get_tools_config()
#(gRNA_bam_file, barcode, output_dir, n_consensus_reads_min, min_umi, auto, pool, ref_fasta, structure_gtf, is_10x) = utils.get_gRNA_mutation_config()

gemgroup_file = 'gemgroup_B.csv'
out = open('consensus.sequence.gRNA.variant.txt','w')

gRNA_all = {}
with open(gemgroup_file, 'r') as fi:
	num=0
	for l in fi:
		ls = l.strip().split(',')
		group = ls[0]
		group_folder = ls[1]
		num = num + 1
		group_file = group_folder + '/consensus.sequence.gRNA.variant.txt'

		print(num)
		with open(group_file, 'r') as fi2:
			for l2 in fi2:
				ls2 = l2.strip().split()
				cb = ls2[0]
				gRNA = ls2[8]
				gRNA_variant = ls2[9] # ALKBH1_gRNA1_gene_variant_1:21M1D21M2D48M
				gRNA_mutation = ls2[10] # TSO:26D(1)
				if num == 1: # initialise by input the first consensus.sequence.gRNA.variant.txt, generate gRNA_all[gRNA][gRNA_mutation] = gRNA_type
					ls2[0] = cb + '-' + str(num)
					out.write('\t'.join(ls2))
					out.write('\n')
					if gRNA == 'multiple':
						continue
					if gRNA_variant == 'WildType':
						continue
					else: 
						gRNA_type = ':'.join(gRNA_variant.split(':')[0:-1]) # ALKBH1_gRNA1_gene_variant_1
					if gRNA not in gRNA_all.keys():
						gRNA_all[gRNA] = {gRNA_mutation:gRNA_type}
					else:
						if gRNA_mutation in gRNA_all[gRNA].keys():
							#if gRNA_type != gRNA_all[gRNA][gRNA_mutation]: # dispendency caused by soft cliping
							continue
						else:
							gRNA_all[gRNA][gRNA_mutation] = gRNA_type
				else: # add gem group
					ls2[0] = cb + '-' + str(num)
					gRNA_type = ':'.join(gRNA_variant.split(':')[0:-1])
					gRNA_type_2 = gRNA_variant.split(':')[-1]
					#gRNA_mutation = gRNA_variant.split(':')[-1]
					gRNA_mutation = ls2[10]
					if gRNA_variant == 'WildType':
						gRNA_variant
						out.write('\t'.join(ls2))
						out.write('\n')
						continue
					if gRNA not in gRNA_all.keys():
						gRNA_all[gRNA] = {gRNA_type: gRNA_mutation}
						out.write('\t'.join(ls2))
						out.write('\n')
					else:
						if gRNA_mutation in gRNA_all[gRNA].keys():
							ls2[9] = gRNA_all[gRNA][gRNA_mutation]+':'+ gRNA_type_2
							out.write('\t'.join(ls2))
							out.write('\n')
						else:
							new_num = 0
							gRNA_types = []
							gRNA_types = gRNA_all[gRNA].values()
							while True:
								new_num = new_num +1
								gRNA_type = gRNA + '_variant_' + str(new_num)
								if(gRNA_type not in gRNA_types):
									break
							gRNA_all[gRNA][gRNA_mutation] = gRNA_type # new variant
							ls2[9] =  gRNA_type + ':' + gRNA_type_2
							#gRNA_mutation
							out.write('\t'.join(ls2))
							out.write('\n')

exit()

####### Assign gRNAs to cells ###########
#assign_gRNA.assign_gRNA_to_cell(in_file = output_dir + 'consensus.sequence.gRNA.variant.txt', min_umi = min_umi, output_dir = output_dir, auto = auto, pool = pool)
#assign_gRNA.add_variant_type(in_file1 = output_dir + 'consensus.sequence.gRNA.variant.txt', in_file2 = output_dir + 'cells.gRNA.single.txt', in_file3 = output_dir + 'cells.gRNA.txt', structure_gtf = structure_gtf, output_dir = output_dir)


