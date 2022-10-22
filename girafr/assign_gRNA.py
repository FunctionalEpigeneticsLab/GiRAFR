#!/usr/bin/env python

from . import utils
import numpy as np
import re


def generate_matrix(in_file, output_dir):
	"""
	Generate matrix (row: features, column: cells) from consensus sequence txt

	Args:
		in_file: consensus.sequence.gRNA.variant.txt

	Returns:
		Write consensus.sequence.matrix

	"""
	if in_file.split('/')[-1] != 'consensus.sequence.gRNA.variant.txt':
		print('Input error: Input needs to be consensus.sequence.gRNA.variant.txt for assign gRNA to cell')
		exit()

	out = open(output_dir + 'consensus.count.matrix', 'w')

	(cells, cb_all, gRNA_all) = utils.read_consensus_sequence_txt(in_file) # cells[cb][gRNA] = {gRNA_type: UMI count}, cb_all: all cell barcodes exists in consensus.sequence.gRNA.variant.txt, gRNA_all: all gRNA types including mutant types

	#cells[list(cells.keys())[0]] = {'HAT1_gRNA3_gene': {'WT': 1}, 'HDAC2_gRNA3_gene': {'WT': 1}, 'CTRL00545_gene': {'CTRL00545_gene_variant_1': 1}, 'RPL18_gRNA1_gene': {'WT': 1}, 'SMARCD1_gRNA4_gene': {'WT': 1}}

	gRNA_all = sorted(gRNA_all)	
	header = 'CellBarcode\t' + ('\t').join(cb_all) + '\n'
	out.write(header)
	for gRNA in gRNA_all:
		out_line = gRNA
		for cb in cb_all:
			if gRNA in cells[cb].keys(): # Retrieve Wild type UMI counts
				if 'WT' in cells[cb][gRNA]:
					out_line += '\t' + str(cells[cb][gRNA]['WT'])
				else: # no WT of this gRNA in this cell
					out_line += '\t' + str(0)
			else: # gRNA not WT
				is_found = False # is_found
				for gRNA_main in cells[cb].keys(): # walk all gRNA types
					if gRNA in cells[cb][gRNA_main].keys() and is_found == False:
						assert gRNA != 'WT'
						is_found = True
						out_line += '\t' + str(cells[cb][gRNA_main][gRNA])
					elif gRNA in cells[cb][gRNA_main].keys() and is_found == True:
						print('Multiple records of', gRNA, 'in cells dictionary')
						exit()
				if is_found == False:
					out_line += '\t' + str(0)
		out.write(out_line + '\n')

def detect_umi_threshold(umi_counts):
	"""
		Fit mix gaussian model to determine umi threshold
	
	Args:
		umi_counts:  gRNA(including variant) UMI counts
	
	Return:
		min_umi: >= min_umi to have this gRNA in cell
	"""
	log_umi = np.log(1. + umi_counts)
	log_umi = np.reshape(log_umi, (len(log_umi), 1))

	gmm = utils.fit_gmm(log_umi) # fit gmm model (2 components)
	if not gmm.converged_:
		print("Warning: model does not converge")
		exit()	
	if np.unique(gmm.predict(log_umi)).size == 1: # only 1 cluster fit
		print("1 model fit gmm") #TODO
		eixt()
	else:
		min_umi = min(umi_counts[gmm.predict(log_umi) == np.argmax(gmm.means_)]) # min_umi: the minum umi of distribution with higher mean
	
	min_umi = int(min_umi)
	return min_umi


def detect_gRNA_umi_thresolds(output_dir, in_file = 'consensus.count.matrix'):
	"""
	
	Args:
		output_dir
		in_file: consensus.count.matrix
		
	Returns:
		gRNA_min_umi = {gRNA: min_umi}
		Write gRNA.umi.threshold.txt

	"""
	gRNA_min_umi = {}
	out = open(output_dir + 'gRNA.umi.threshold.txt', 'w')
	with open(in_file, 'r') as fi:
		for l in fi:
			l = l.strip()
			if l[0:11] == 'CellBarcode': 
				continue # skip cell barcode line
			else:
				ls = l.split()
				gRNA = ls[0].strip()
				
				gRNA_umi_counts = np.array([])
				for i in ls[1:]:
					gRNA_umi_counts = np.append(int(i), gRNA_umi_counts)
				
				min_umi = detect_umi_threshold(gRNA_umi_counts) # min_umi threshold >= min_umi
				if gRNA not in gRNA_min_umi.keys():
					gRNA_min_umi[gRNA] = min_umi
					out_line = gRNA + '\t' + str(min_umi)
					out.write(out_line + '\n')
				else:
					print("ERR: Duplicated gRNA: ",gRNA, "\n Check consensus.sequence.matrix")
					exit()
	
	return gRNA_min_umi

	
def detect_gRNA_umi_thresolds_pool(output_dir, in_file = 'consensus.count.matrix'):
	"""
	 calculate min umi thresholds together with variant gRNA of the same guide	
	Args:
		output_dir
		in_file: consensus.count.matrix
		
	Returns:
		gRNA_min_umi = {gRNA: min_umi}
		Write gRNA.umi.threshold.txt

	"""
	gRNA_min_umi = {}
	out = open(output_dir + 'gRNA.umi.threshold.txt', 'w')
	with open(in_file, 'r') as fi:
		
		for l in fi:
			l = l.strip()
			if l[0:11] == 'CellBarcode': 
				ls = l.split()
				cell_number = len(ls) - 1
				gRNA_target = '#' # initial gRNA_target as false gene name 
				gRNA_umi_counts = np.zeros(cell_number) # initial gRNA umi counts
				continue # skip cell barcode line
			else:
				ls = l.split()
				gRNA = ls[0].strip() # eg. ATAD2_gRNA1_gene_variant_10 when is_MT is True, ATAD2_gRNA1_gene when is_MT is False
				is_MT = 'variant' in gRNA # True: MT, False: WT

				if  gRNA_target not in gRNA: # current gRNA variant is not for the same gRNA target gene (includes no wild type for this guide exists)
					if not (gRNA_umi_counts == np.zeros(cell_number)).all(): # calculate previous gRNA target
						#print('calculate previous gRNA targte:', gRNA_target)
						#print(gRNA_target, end = '\t')
						#for x in gRNA_umi_counts:
						#	print(x, end = '\t')
						#print('\n')
						min_umi = detect_umi_threshold(gRNA_umi_counts) # min_umi threshold >= min_umi
						if gRNA_target not in gRNA_min_umi.keys():
							gRNA_min_umi[gRNA_target] = min_umi
							out_line = gRNA_target + '\t' + str(min_umi)
							out.write(out_line + '\n')
						else:
							print("ERR: Duplicated gRNA: ",gRNA_target, "\n Check consensus.sequence.matrix")
							exit()
					
					gRNA_target = gRNA if not is_MT else gRNA[:re.search('variant',gRNA).span()[0]-1] # update new gRNA_target gene name

					gRNA_umi_counts = np.zeros(cell_number) # initial gRNA umi counts

					for idx, val in enumerate(ls[1:]):
						if int(val) > gRNA_umi_counts[idx]:
							gRNA_umi_counts[idx] = int(val)

				else: # variant, gRNA_target in gRNA
					assert len(gRNA_umi_counts) != 0
					for idx, val in enumerate(ls[1:]):
						if int(val) > gRNA_umi_counts[idx]:
							gRNA_umi_counts[idx] = int(val)
		
		print('calculate previous gRNA targte:', gRNA_target)
		min_umi = detect_umi_threshold(gRNA_umi_counts) # min_umi threshold >= min_umi
		if gRNA_target not in gRNA_min_umi.keys():
			gRNA_min_umi[gRNA_target] = min_umi
			out_line = gRNA_target + '\t' + str(min_umi)
			out.write(out_line + '\n')
		else:
			print("ERR: Duplicated gRNA: ",gRNA_target, "\n Check consensus.sequence.matrix")
			exit()
	
	return gRNA_min_umi


def assign_gRNA_to_cell(in_file, output_dir, min_umi = 3, auto = False, pool = False):
	"""
	
	Args:
		in_file: consensus.sequence.gRNA.variant.txt
		min_umi: minum number of UMI, default is 3 
		auto: boolean, whether use autodetection or fixed min_umi
		pool: boolean, whether calculate min umi thresholds together with variant gRNA of the same guide
	Returns:
		Write cells.gRNA.txt and cells.gRNA.single.txt

	"""

	generate_matrix(in_file, output_dir) # generate consensus.sequence.matrix 

	out = open(output_dir + 'cells.gRNA.txt','w')
	out2 = open(output_dir + 'cells.gRNA.single.txt','w')
	
	if auto: # if use min_umi auto detection
		if pool: # if calculate min umi thresholds together with variant gRNA of the same guide
			gRNA_min_umi = detect_gRNA_umi_thresolds_pool(output_dir)
			print('gRNA umi thresold auto and pool detection is enabled')
		else:
			gRNA_min_umi = detect_gRNA_umi_thresolds(output_dir)
			print('gRNA umi thresold auto detection is enabled')

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
			# add WT gRNA umi count
			if 'WT' not in cells[cb][gRNA].keys():
				WT_umi = 0
			else:
				try:
					min_umi = gRNA_min_umi[gRNA] if auto else min_umi # wildtype gRNA min_umi
				except KeyError:
					print(gRNA, 'not found in gRNA_min_umi dictionary\n')
					exit()
				if cells[cb][gRNA]['WT'] >= min_umi:
					WT_umi = cells[cb][gRNA]['WT']
					n_gRNA += 1
					has_gRNA = True
				else:
					WT_umi = 0
			# add Variant gRNA umi count
			for gRNA_type in cells[cb][gRNA].keys():
				if gRNA_type != 'WT':
					try:
						if auto:
							# variant gRNA min_umi
							if pool:
								gRNA_target = gRNA_type[:re.search('variant', gRNA_type).span()[0]-1]
								assert gRNA_type not in gRNA_min_umi.keys()
								min_umi = gRNA_min_umi[gRNA_target] # variant gRNA min_umi saved together with WT
							else:
								min_umi = gRNA_min_umi[gRNA_type] 

					except KeyError:
						print(gRNA_type, 'not found in gRNA_min_umi dictionary\n')
						exit()
					if cells[cb][gRNA][gRNA_type] >= min_umi:
						n_gRNA += 1
						has_gRNA = True
						out_line += '\t'+ gRNA_type + '\t'+ str(cells[cb][gRNA][gRNA_type])+ '\t'
			if has_gRNA > 0 :
				out_line = cb + '\t' + str(n_gRNA) + '\t' + gRNA + '\t' + 'WT' + '\t' + str(WT_umi) + out_line
				out.write(out_line + '\n')
		
		if n_gRNA == 1: # for cell only contains single gRNA
			for gRNA in cells[cb].keys():
				for gRNA_type in cells[cb][gRNA].keys():
					if gRNA_type == 'WT': # wild type
						min_umi = gRNA_min_umi[gRNA] if auto else min_umi # wildtype gRNA min_umi
						if cells[cb][gRNA][gRNA_type] >= min_umi:
							out2.write(cb + '\t' + str(n_gRNA) + '\t' + gRNA + '_' + gRNA_type + '\t'+ str(cells[cb][gRNA][gRNA_type]) + '\n')
					else: # variant type
						if auto:
							if pool:
								gRNA_target = gRNA_type[:re.search('variant', gRNA_type).span()[0]-1]
								assert gRNA_type not in gRNA_min_umi.keys()
								min_umi = gRNA_min_umi[gRNA_target]
							else:
								min_umi = gRNA_min_umi[gRNA_type]

						if cells[cb][gRNA][gRNA_type] >= min_umi:
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
				if variants[gRNA][0] != gRNA_MT:
					print('gRNA', gRNA)
					print('gRNA_MT', gRNA_MT)
					print('variants[gRNA][0]:', variants[gRNA][0])
					print(variants[gRNA])
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
			for index in range(5, len(ls), 2): # fix bug 30 Sep
				gRNA = ls[index]
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
