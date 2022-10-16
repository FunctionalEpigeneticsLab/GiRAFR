#!/usr/bin/env python

import gzip
from . import utils
import sys
import ast



def detect_editing_effect(consensus_in_file, cut_site, cell_file):

	"""Detect edditing effect
	Arg:
		consensus_in: consensus.sequence.txt generted from consensus_sequence.generate_consensuse_sequence
		cell_file: file that contain cell - sgRNA annotation
		For cellRanger feature barcoding analysis, cell_file is protospacer_calls_per_cell.csv
		or cells.gRNA.txt generated by previous step or file with
		Column 1: cell barcode
		Column 3: sgRNA asisgned (new row if multiple sgRNA/cell)
			
	"""
	cut_site = int(cut_site) - 1 # 1-based convert to 0-based
	KO = utils.read_cell_assignment(cell_file)
	consensus_in = open(consensus_in_file, 'r')
	out_file = consensus_in_file.replace('consensus.sequence.txt','editing_effect.cutsite.txt')
	out = open(out_file, 'w')
	
	for l in consensus_in:
		l = l.strip()
		if l.startswith('#'):
			continue  #skip header
		ls = l.split("\t")
		cb = ls[0]
		umi = ls[1]
		number_of_UMI_for_cell = ls[2]
		number_of_UMI_sequence_for_umi = ls[3]
		union_seq = ls[5]
		union_nreads = ls[6]
		is_WT = ls[7]	
		gene = ls[8]
		union_positions = ast.literal_eval(ls[9])

		if union_seq == 'None': # remove records that do not cover this region
			continue

		if cb in KO.keys():
			KO_gene = str(KO[cb])
		else:
			KO_gene = 'empty'

		if is_WT == 'True':
			result = l + "\t" + "No mutation\t" + KO_gene + "\n"
			out.write(result)
			continue
		#print(l)
		#print(union_positions)
		#print(cut_site)
		
		# new version
		left_pos = union_positions[0]
		right_pos = union_positions[-1]
		i = 0
		while left_pos == None:
			i = i + 1
			left_pos = union_positions[i]
		i = -1
		while right_pos == None:
			i = i - 1
			right_pos = union_positions[i]

		#print(union_positions)
		#print(left_pos, right_pos)
		
		if left_pos < cut_site and right_pos > cut_site: # over cut site #TODO union_positions[0] <= cut_site and union_positions[-1] >= cut_site????
			if cut_site not in union_positions or cut_site+1 not in union_positions or cut_site-1 not in union_positions:
				result = l + "\t" + "Deletion at cutsite\t" + KO_gene + "\n"
				out.write(result) # must mean insertion or deletion
			elif right_pos - left_pos +1 != len(union_positions):
				result = l + "\t" + "Insertion/Deletion somewhere else\t"  + KO_gene + "\n"
				out.write(result)
			elif right_pos - left_pos +1 == len(union_positions):
				result = l + "\t" + "Mismatch\t"  + KO_gene + "\n"
				out.write(result)
			else:
				print("unseen type #TODO")
				exit()
		else: # sequence does not cover cut site
			result = l + "\t" + "No cover\t"  + KO_gene + "\n"
			out.write(result)


		# old version:
		#if union_positions[0] < cut_site and union_positions[-1] > cut_site: # over cut site #TODO union_positions[0] <= cut_site and union_positions[-1] >= cut_site????
		#	if cut_site not in union_positions or cut_site+1 not in union_positions or cut_site-1 not in union_positions:
		#		result = l + "\t" + "Deletion at cutsite\t" + KO_gene + "\n"
		#		out.write(result) # must mean insertion or deletion?
		#	elif union_positions[-1] - union_positions[0] +1 != len(union_positions):
		#		result = l + "\t" + "Insertion/Deletion somewhere else\t"  + KO_gene + "\n"
		#		out.write(result)
		#	elif union_positions[-1] - union_positions[0] +1 == len(union_positions):
		#		result = l + "\t" + "Mismatch\t"  + KO_gene + "\n"
		#		out.write(result)
		#	else:
		#		print("unseen type #TODO")
		#		exit()
		#else: # sequence does not cover cut site
		#	result = l + "\t" + "No cover\t"  + KO_gene + "\n" 
		#	out.write(result)

