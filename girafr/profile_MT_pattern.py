#!/usr/bin/env python

from . import utils
import re


def profile_MT_pattern(in_file, min_umi2, output_dir):
	"""
	Args:
		in_file: consensus.sequence.gRNA.variant.txt
		min_umi2: ninum number of umi for certain Variant to be profiled
	Returns:
		Write consensus.sequence.gRNA.MT.txt
	"""
	out = open(output_dir + 'consensus.sequence.gRNA.MT.txt', 'w')
	variants = {} # variants[gRNA] 
	with open(in_file, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			cb = ls[9]
			gRNA = ls[8]

			if gRNA == 'multiple':
				continue
			gRNA_type = ls[9] # ALKBH1_gRNA1_gene_variant_1:21M1D21M2D48M
			if gRNA_type != 'WildType':
				if len(ls) < 11: continue # TODO,for variant like 90M1S, soft  cliping here is skipped
				#gRNA = gRNA_type.split(':')[0] # can have problem if gene name has ':'
				gRNA = ':'.join(gRNA_type.split(':')[0:-1]) # xx_gene_Variant1
				gRNA_MT = ls[10]
			else:
				continue # skip wild type UMIs
			if gRNA not in variants:
				variants[gRNA] = [gRNA_MT, 1]
			else:
				if variants[gRNA][0] != gRNA_MT:
					print('gRNA', gRNA)
					print('variants[gRNA][0]:',variants[gRNA][0])
					print('records do not match')
					exit()
				else:
					variants[gRNA][1] = variants[gRNA][1] + 1 # add #UMI
			
	structure = utils.read_annotation_structure(structure_gtf)
	print(len(variants.keys()))
	for key in variants:
		gRNA = key
		if variants[key][1] >= min_umi2: # at least min_umi2 for this variants across the entire library
			mutation = variants[gRNA][0] # gRNA:6I(1),Rest:8D(1)
			if bool(re.search('_WT', gRNA)):
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
		else:
			continue

(samtools, twoBitToFa, featureCounts) = utils.get_tools_config()
(gRNA_bam_file, barcode, output_dir, n_consensus_reads_min, min_umi, auto, pool, ref_fasta, structure_gtf, is_10x) = utils.get_gRNA_mutation_config()
profile_MT_pattern(in_file = output_dir + 'consensus.sequence.gRNA.variant.txt', min_umi2 = 2, output_dir = output_dir)
