#!/usr/bin/python

import utils
import re

def assign_gRNA_to_cell(in_file, min_umi, output_dir):
	"""
	
	Args:
		in_file: consensus.sequence.gRNA.variant.txt
		min_umi: minum number of UMI 
	Returns:
		Write cells.gRNA.txt and cells.gRNA.single.txt

	"""
	out = open(output_dir + 'cells.gRNA.txt','w')
	out2 = open(output_dir + 'cells.gRNA.single.txt','w')
	cells = {}
	with open(in_file, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			cb = ls[0]
			gRNA = ls[8]

			if gRNA == 'multiple':
				continue #TODO

			gRNA_type = ls[9] # ALKBH1_gRNA1_gene_variant_1:21M1D21M2D48M
			if cb not in cells.keys():
				if gRNA_type == 'WildType':
					cells[cb] = {gRNA: {'WT': 1}}

				else:
					#gRNA_type = gRNA_type.split(':')[0] # can have problem if gene name has ':'
					gRNA_type = ':'.join(gRNA_type.split(':')[0:-1]) 
					cells[cb] = {gRNA: {gRNA_type: 1}}
			else:
				if gRNA not in cells[cb].keys():
					if gRNA_type == 'WildType':
						cells[cb][gRNA] = {'WT': 1}
					else:
						#gRNA_type = gRNA_type.split(':')[0]
						gRNA_type = ':'.join(gRNA_type.split(':')[0:-1])
						cells[cb][gRNA] = {gRNA_type: 1}
				else:
					if gRNA_type == 'WildType':
						if 'WT' not in cells[cb][gRNA].keys():
							cells[cb][gRNA]['WT'] = 1
						else:
							cells[cb][gRNA]['WT'] += 1
					else:
						#gRNA_type = gRNA_type.split(':')[0]
						gRNA_type = ':'.join(gRNA_type.split(':')[0:-1])
						if gRNA_type not in cells[cb][gRNA].keys():
							cells[cb][gRNA][gRNA_type] = 1
						else:
							cells[cb][gRNA][gRNA_type] += 1
	fi.close()
	for cb in cells.keys():
		n_gRNA = 0
		for gRNA in cells[cb].keys():
			WT_umi = 0
			has_gRNA = False
			out_line = ''
			if 'WT' not in cells[cb][gRNA].keys():
				WT_umi = 0
			else:
				if cells[cb][gRNA]['WT'] >= min_umi:
					WT_umi = cells[cb][gRNA]['WT']
					n_gRNA += 1
					has_gRNA = True
				else:
					WT_umi = 0
			for gRNA_type in cells[cb][gRNA].keys():
				if gRNA_type != 'WT':
					if cells[cb][gRNA][gRNA_type] >= min_umi:
						n_gRNA += 1
						has_gRNA = True
						out_line += '\t'+ gRNA_type + '\t'+ str(cells[cb][gRNA][gRNA_type])+ '\t'
			if has_gRNA > 0 :
				out_line = cb + '\t' + str(n_gRNA) + '\t' + gRNA + '\t' + 'WT' + '\t' + str(WT_umi) + out_line
				out.write(out_line + '\n')
		if n_gRNA == 1:
			for gRNA in cells[cb].keys():
				for gRNA_type in cells[cb][gRNA].keys():
					if cells[cb][gRNA][gRNA_type] >= min_umi:
						if gRNA_type == 'WT':
							out2.write(cb + '\t' + str(n_gRNA) + '\t' + gRNA + '_' + gRNA_type + '\t'+ str(cells[cb][gRNA][gRNA_type]) + '\n')
						else:
							out2.write(cb + '\t' + str(n_gRNA) + '\t' + gRNA_type + '\t'+ str(cells[cb][gRNA][gRNA_type]) + '\n')


def add_variant_type(in_file1, in_file2, in_file3, structure_gtf, output_dir):
	"""
	Arg:
		in_file1: consensus.sequence.gRNA.variant.txt
		in_file2: cells.gRNA.single.txt
		in_file3: cells.gRNA.txt

	Return:
		Write cells.gRNA.single.MT.txt
		MT.txt
		all.MT.txt: Number of variants can be less than cells.gRNA.txt as softclipping ones are not counted here.
	"""
	variants = {}
	with open(in_file1, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			gRNA = ls[8]
			seq = ls[5]
			if gRNA == 'multiple': # TODO
				continue
			gRNA_type = ls[9]
			if gRNA_type != 'WildType':
				if len(ls) < 11: continue # TODO,for variant like 90M1S, soft  cliping here is skipped
				#gRNA = gRNA_type.split(':')[0] # can have problem if gene name has ':'
				gRNA = ':'.join(gRNA_type.split(':')[0:-1])
				gRNA_MT = ls[10]
			else:
				gRNA = gRNA + '_WT' # ALKBH1_gRNA1_gene_WT
				gRNA_MT = None
				seq = None # As different wiltype reads are not exact match, some might have a bit more sequence to the upstream or downstream
			if gRNA not in variants:
				variants[gRNA] = [gRNA_MT, seq]
			else:
				if variants[gRNA][0] != gRNA_MT or variants[gRNA][1] != seq:
					print('gRNA', gRNA)
					print('variants[gRNA][0]:',variants[gRNA][0])
					print('gRNA_MT', gRNA_MT)
					print('variants[gRNA][1]', variants[gRNA][1])
					print('seq', seq)
					print('records do not match')
					exit()
				else:
					continue

	fi.close()
	structure = utils.read_annotation_structure(structure_gtf)

	out1 = open(output_dir + 'cells.gRNA.single.MT.txt','w')
	out2 = open(output_dir + 'MT.txt', 'w')	
	out3 = open(output_dir + 'all.MT.txt', 'w')	
	with open(in_file2, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			gRNA = ls[2]
			if gRNA not in variants: #TODO, CHD1L_gRNA2_gene_variant_20  64M27S
				out_line = l.strip() + '\t' + str(None) + '\t' + str(None)
				out1.write(out_line + '\n')
				continue
			else:
				out_line = l.strip() + '\t' + str(variants[gRNA][0]) + '\t' +str(variants[gRNA][1])
				out1.write(out_line + '\n')
			
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
					out2.write(out_line + '\n')

	with open(in_file3, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			if len(ls) == 5: # only WT
				continue
			gRNA = ls[5]
			
			if gRNA not in variants: #TODO, CHD1L_gRNA2_gene_variant_20  64M27S
				continue
			else:
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
					out3.write(out_line + '\n')
				
				variants.pop(gRNA) # remove the gRNA variant which has been written
