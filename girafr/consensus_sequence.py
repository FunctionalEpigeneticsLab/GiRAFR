#!/usr/bin/env python

import gzip
from . import utils
import sys
import pysam

def generate_consensus_sequence_gRNA(bam_in_file, barcode_in, output_dir, n_consensus_reads_min = 1, is_10x = True):
	"""
	########### gRNA consensus sequence ########
	Generate consensus sequence for gRNA library

	Arg:
		bam_in: bam file, alignment file of gRNA library after removed secondary alignment and mapped not on gRNA reference
		barcodes: filtered barcode list
	return:
		Write consensus.sequence.gRNA.txt: consensus sequence supported by column 7 (n_consensus_reads) > 1
		      consensus.bam: consensus sequence in bam file format
		Wirte consensus.seqeunce.gRNA.all_umi.txt: all UMI detected consensus sequence 
		      Non-consensus.bam: alignment not the same as consensus sequence

	"""
	cb_tag = 'CB' if is_10x else 'XC'
	umi_tag = 'UB' if is_10x else 'XM'
	gene_tag = 'GN'if is_10x else 'gn'

	barcodes = utils.read_barcode(barcode_in) # filtered barcode in list
	print('Total filtered barcodes: ', len(barcodes))
	bam_in = utils.create_bam_infile(bam_in_file)
	out_file = output_dir + 'consensus.sequence.gRNA.txt'
	out_file2 = output_dir + 'consensus.sequence.gRNA.all_umi.txt'
	out = open(out_file,'w')
	out2 = open(out_file2, 'w')

	if is_10x: # cell barcodes discrepancy between gene expression library and crispr capture library, especially for feature barcoding technology
		lookup_barcodes = utils.lookup_barcodes()

	dict_2d={}
	# Read all records
	for r in bam_in:
		read_id = r.query_name
		seq = r.query_sequence
		cb = utils.get_read_tag(r, cb_tag) # for now, assume its after corrected ### TODO
		umi = utils.get_read_tag(r, umi_tag) # for now, assume its after corrected ### TOD
		gene = utils.get_read_tag(r, gene_tag)
		cigar = r.get_cigar_stats()
		is_WT = utils.define_WT(r, is_10x)
		if cb == None: # empty
			continue
		if umi == None: # no UMI
			continue
		if is_10x & bool(cb):
			cb = cb.replace('-1','')
		if cb not in barcodes:
			if is_10x and cb in lookup_barcodes.keys():
				cb1 = lookup_barcodes[cb] 
				if cb1 in barcodes: # really not  filtered cell barcodes
					cb = cb1 # corrected cell barcodes =>new CB
				else:
					continue	
			else:
				continue # skip cells which are not in filtered barcode list
		if cb in dict_2d.keys():
			if umi in dict_2d[cb].keys():
				dict_2d[cb][umi].append((seq, gene, cigar, is_WT))
			else:
				dict_2d[cb][umi]=[(seq, gene, cigar, is_WT)]
		else:
			dict_2d[cb] = {umi: [(seq, gene, cigar, is_WT)]}
	
	print(str(len(dict_2d.keys())) + "\tcells in filtered barcodes list")
	header = "#\t"+str(len(dict_2d.keys()))+"\tCells processed\n"
	out.write(header)
	
	consensus_2d = {}
	for cb in dict_2d.keys():
		#print(len(dict_2d[cb].keys())) # number of UMI in the cell
		for umi in dict_2d[cb].keys():
			#dict_umi_seq = utils.collapse_umi(dict_2d[cb][umi])

			gene = utils.collapse_mapped_gene(dict_2d[cb][umi])
			gene_mapped = utils.which_mapped_gene(gene) # When multiple aligned genes, take the one with most reads. 
			dict_umi_seq = utils.collapse_umi_v2(dict_2d[cb][umi], gene_mapped) # count reads that mapped tp gene_mapped
			
			try:
				max_seq = max(dict_umi_seq, key = dict_umi_seq.get) # most frequent sequence
			except TypeError:
				print(TypeError)
				print(umi,cb)
				print(dict_2d[cb][umi])
				continue
			
			umi_class = utils.classify_umi(dict_umi_seq)

			n_consensus_reads = dict_umi_seq[max_seq][0] # number of reads that support consensus sequence
			result = str(cb) + "\t" + str(umi) + "\t" + str(len(dict_2d[cb].keys())) + "\t" + str(len(dict_2d[cb][umi])) + "\t" + str(umi_class) + "\t"+ str(max_seq) + "\t" + str(n_consensus_reads) + "\t" + str(dict_umi_seq[max_seq][1]) + "\t" + str(gene_mapped) + "\n"

			if n_consensus_reads > n_consensus_reads_min:
				out.write(result) # write consensus.seqeunce.gRNA.txt
			out2.write(result) # write consensus.seqeunce.gRNA.all_umi.txt

			consensus_seq = str(max_seq)
			is_WT = str(dict_umi_seq[max_seq][1])

			if cb in consensus_2d.keys():
				if umi in consensus_2d[cb].keys():
					print('duplicated umi in same cell:',cb,umi)
					exit()
				else:
					consensus_2d[cb][umi] = [consensus_seq, is_WT] 
			else:
				consensus_2d[cb] = {umi: [consensus_seq, is_WT]}
	out.close()
	out2.close()

	bam_in = utils.create_bam_infile(bam_in_file)
	bam_out = pysam.AlignmentFile(output_dir + 'consensus.bam','wb', template = bam_in)
	bam_out2 = pysam.AlignmentFile(output_dir + 'Non-consensus.bam','wb', template = bam_in)
	flag = False # mean this value has been recorded
	for r in bam_in:
		seq = r.query_sequence	
		cb = utils.get_read_tag(r, cb_tag) 
		umi = utils.get_read_tag(r, umi_tag) 
		if is_10x & bool(cb):
			cb = cb.replace('-1','')
		if cb not in consensus_2d.keys(): 
			if is_10x and cb in lookup_barcodes.keys():
				cb1 = lookup_barcodes[cb] # corrected cell barcodes =>new CB 
				if cb1 in consensus_2d.keys(): # really not  filtered cell barcodes
					cb = cb1
					r.set_tag('UB', cb)
				else:
					continue
			else:
				continue # cell barcode not in metadata, meaning its not in filtered barcode list
		if umi == None: # no UMI
			continue

		consensus_seq = consensus_2d[cb][umi][0]
		
		is_WT = consensus_2d[cb][umi][1]
		r.tags += [('WT',is_WT)]
		
		if bool(consensus_seq) & (seq == consensus_seq): # if this query sequence is same as consensus sequence
			bam_out.write(r) # write alignment to consensus.bam
			consensus_2d[cb][umi][0] = flag # this umi in this cell has been recorded with consensus sequence in output
		elif bool(consensus_seq) == False:
			continue
		elif seq != consensus_seq: # not the same as consensus sequence
			bam_out2.write(r) # write alignment to Non-consensus.bam
		else:
			print('What situation is this?')
			exit()

	

def generate_consensus_sequence(bam_in_file, barcode_in, region_chr, region_start, region_end, region_ref, ref_gene, ref_gRNA_id, is_10x = True):

	"""Generate consensus sequence
	
	Arg:
		bam_in: bam file split by samtools by detection window, already processed by samtools view region + umi_tools eg. Region629.bam.featureCounts.sorted.bam
		barcodes: filtered barcode list 
		region_chr: region chromosome eg. chr15
		region_start: 1-based window start eg. 72208851
		region_end: 1-based window end eg. 72208881
		region_ref: reference sequence in detection window
		ref_gene: target gene, gene symbol eg. PKM
		ref_gRNA_id: eg. PKM_gRNA2
		is_10x: whether bam file generated by CellRanger or Dropseq toolbox (True by default)
	return:
		Write consensus.sequence.txt
		
	"""
	cb_tag = 'CB' if is_10x else 'XC'
	umi_tag = 'UB' if is_10x else 'XM'
	gene_tag = 'XT' # assigned by umi_tools

	region_start = int(region_start) - 1
	region_end = int(region_end) - 1
	region_ref = region_ref.upper()
	
	barcodes = utils.read_barcode(barcode_in) # filtered barcode in list
	bam_in = utils.create_bam_infile(bam_in_file) # already processed by samtools view region + umi_tools eg. Region629.bam.featureCounts.sorted.bam 

	out_file = bam_in_file.replace('bam.featureCounts.sorted.bam','consensus.sequence.txt')
	out = open(out_file,'w')

	dict_2d = {}
	# Read all records
	for r in bam_in:
		# r is AlignedSegment
		#print(r.query_name) # Read ID
		#print(r.get_reference_sequence()) # reference sequence
		#print(r.query_sequence) # query sequence
		#print(r.get_tag('XC')) # cell barcode
		#print(r.get_tag('XM')) # umi

		read_id = r.query_name
		seq = r.query_sequence
		cb = utils.get_read_tag(r, cb_tag) # for now, assume its after corrected ### TODO
		umi = utils.get_read_tag(r, umi_tag) # for now, assume its after corrected ### TODO
		gene = utils.get_read_tag(r, gene_tag)
		cigar = r.get_cigar_stats()
		is_rev = r.is_reverse # if reads mapped to reverse strand
		# for reads that mapped to reverse strand, we treat them normally since records have been reverse complemented in bam file
	
		mapped_positions = r.get_reference_positions(full_length=True) # [71893693,...] 1-based
		number_overlapped = r.get_overlap(region_start, region_end + 1)  # return number of aligned bases of read overlapping the interval start and end on the reference sequence.

		seq_overlapped = ''
		positions_overlapped = []
		is_WT = True
		left_overlapped_pos = region_end
		left_overlapped_i = None
		right_overlapped_pos = 0
		right_overlapped_i = None
		for i,pos in enumerate(mapped_positions):
			if pos == None: # None values will be included for any soft-clipped or unaligned positions within the read.
			# What to do if D/I exist within 31bp region
			# when clipped => None as pos => these reads included or not? TODO
			# when there is a D => reads skip N Deletion bases on ref => pos will be within ref => results should be fine
			# when there is a I => reads with additional bases inserted => pos will be none  => TODO
				continue
			#if pos >= region_start and pos < region_end:
			#	seq_overlapped += seq[i]
			#	positions_overlapped.append(pos)
			#	if seq[i] != region_ref[pos - region_start]:
			#		is_WT = False
			if pos >= region_start and pos <= region_end:
				if pos <= left_overlapped_pos:
					left_overlapped_pos = pos
					left_overlapped_i = i
				if pos >= right_overlapped_pos:
					right_overlapped_pos = pos
					right_overlapped_i = i
		
		if left_overlapped_i!= None and right_overlapped_i!=None:
			seq_overlapped = seq[left_overlapped_i:right_overlapped_i+1] # include right_overlapped_o
			positions_overlapped = mapped_positions[left_overlapped_i:right_overlapped_i+1] # include right_overlapped position
			is_WT = True if seq_overlapped == region_ref[left_overlapped_pos - region_start: right_overlapped_pos - region_start + 1] else False

		# if insertion exist (N insertions), then get 31 + N bp sequence instead of 31bp from reference, but overlapped seq is still 31bp maximum. 
		number_overlapped = number_overlapped + positions_overlapped.count(None)
		assert len(seq_overlapped) == number_overlapped
		assert len(positions_overlapped) == len(seq_overlapped)
		if number_overlapped == 0:
			seq_overlapped = None
			is_WT = False
		if cb == None: # empty
			continue
		if umi == None: # no UMI
			continue
		if is_10x & bool(cb):
			cb = cb.replace('-1','')
		if cb not in barcodes:
			continue # skip cells which are not in barcode whitelist

		if gene != ref_gene:
			continue # skip alignemnt that are not mapped on targeted gene
	
		if cb in dict_2d.keys():
			if umi in dict_2d[cb].keys():
				dict_2d[cb][umi].append((seq_overlapped, gene, cigar, is_WT, positions_overlapped))
			else:
				dict_2d[cb][umi]=[(seq_overlapped, gene, cigar, is_WT, positions_overlapped)]
		else:
			dict_2d[cb] = {umi: [(seq_overlapped, gene, cigar, is_WT, positions_overlapped)]}	
	#print(len(dict_2d.keys()), "Cells processed") # number of cell barcde
	header = "#\t"+str(len(dict_2d.keys()))+"\tCells processed\n"+"#cell_barcode\tUMI\tnumber_of_UMI_for_cell\tnumber_of_UMI_sequence_for_umi\tref_seq\tunion_seq\tunion_nreads\tis_WT\tgene\tunion_positions\tis_umi_seq_different\tref_gRNA_id\n"
	out.write(header)

	for cb in dict_2d.keys():
		#print(cb)
		#print(len(dict_2d[cb].keys())) # number of UMI in the cell 
		for umi in dict_2d[cb].keys():
			dict_umi_seq = utils.collapse_umi(dict_2d[cb][umi]) # dict_umi_seq[seq] = [num_reads, is_WT, positions_overlapped]
			#print(dict_umi_seq)
			umi_seq_different = True # whether seq from same umi mapped to same coordinates
			union_seq = ''
			union_positions = [0]
			union_nreads = 0
			is_WT = True
			if False in sum(list(dict_umi_seq.values()),[]):
				is_WT = False
			if len(dict_umi_seq.keys()) == 1: # only one sequence for one UMI
				union_seq = list(dict_umi_seq.keys())[0] if list(dict_umi_seq.keys())[0]!= None else ''
				union_positions = list(dict_umi_seq.values())[0][2]
				union_nreads = list(dict_umi_seq.values())[0][0]
			else:
				all_seq_umi = list(dict_umi_seq.keys()) # whether sequences are the same for one same UMI
				all_seq_umi = [i for i in all_seq_umi if i != None]
				all_seq_umi.sort(key =len, reverse = True) # order by the length of the sequence
				union_seq = all_seq_umi[0] # longest sequence as initial
				union_positions = dict_umi_seq[union_seq][2]
				union_nreads = dict_umi_seq[union_seq][0]
				#print('=========================')
				#print(dict_umi_seq)
				for seq in all_seq_umi[1:]: # for each sequence of one same UMI, union them
					#print("current seq:", seq)
					positions_overlapped = dict_umi_seq[seq][2]
					union_nreads = union_nreads + dict_umi_seq[seq][0]
					for i,pos in enumerate(positions_overlapped):
						#print(i,pos)
						#print("union_seq:",union_seq)
						#print("union_positons:",union_positions)
						if pos == None: # when insertion / cliping happens 
							if i != 0 or i != len(positions_overlapped): # not at the most left or most right
								pos_left = positions_overlapped[i-1]
								if pos_left in union_positions:
									i2 = union_positions.index(pos_left)
									#print('seq:', seq)
									#print('positions_overlapped', positions_overlapped)
									#print('union_seq',union_seq)
									#print('union_positions', union_positions)
									if i2+1 == len(union_seq): # pos_left at the most right of the union_seq
										union_seq = union_seq + seq[i] # add to the right of the union sequence	
										union_positions.append(pos)
										continue
									elif union_seq[i2+1] == None: # if in union_seq this position is also None
										continue
									else: # 
										union_seq = union_seq[:i2] + seq[i] + union_seq[i2:]
										union_positions.insert(i2,pos)
										#print("insertion")
										#print(union_seq)
										#print(union_positions)
										continue
								else:
									sys.exit("Unkown ERROR")
							else:
								print("Insertion or clipling happens at the most left or right => skip")
								continue
						union_left = union_positions[0]
						union_right = union_positions[-1]

						if union_left == None: # begin with insertion/cliping
							#print("Insertion or clipling happens at the most left => skip") # updated 10 July
							#continue
							union_left = union_positions[1]
						
						tmp_i = -1
						while(union_right == None): # end with insertion/cliping
							# caused by previous union (previous pos = None)
							tmp_i = tmp_i - 1
							union_right = union_positions[tmp_i]
						
						if pos < union_left: # left to the sequence
							union_seq = seq[i] + union_seq
							union_positions = [pos] + union_positions # to the most left of union_positions
						elif pos > union_right: # right to the sequence
							union_seq = union_seq + seq[i]
							union_positions.append(pos)
						elif union_left <= pos and pos <= union_right: # test if overlapped part is the same
							if pos in union_positions: # position already in the union_positions
								if seq[i] != union_seq[union_positions.index(pos)]: # NOT same as union_seq # union_positions.index(pos)-1
									#sys.exit('Different UMI sequence') # different sequences for the same UMI
									#print("Different UMI sequence")
									if dict_umi_seq[seq][0] > union_nreads - dict_umi_seq[seq][0]: # if more reads support this sequence
										union_nreads = dict_umi_seq[seq][0]
										union_seq = union_seq[:union_positions.index(pos)] + seq[i] + union_seq[union_positions.index(pos)+1:]# change this position into seq[i]
										union_positions = union_positions # union_positions remain the same
									else:
										break
									umi_seq_different = False
							else: # position not in union_positions
								for i2,pos2 in enumerate(union_positions): # look for where to put in the list, positions should be ordered
									if pos2 == None: # when insertion / cliping happens
										continue
									if pos2 > pos: 
										union_seq = union_seq[:i2] + seq[i] + union_seq[i2:]
										union_positions.insert(i2,pos)
										break
						else:
							print('Unknown situation with ',pos,'union_positions: ',union_positions)
							exit()
						union_positions_tmp1 = [i for i in union_positions if i != None] # cannot order list which contain Nonetype
						union_positions_tmp2 = union_positions_tmp1
						union_positions_tmp1.sort()
						assert union_positions_tmp2 == union_positions_tmp1
						#print("*********"+union_seq)
			union_positions = [i for i in union_positions if i != 0] # remove the first added zero, otherwise the length will not pass assertion
			assert len(union_seq) == len(union_positions)
			if union_seq == '':
				union_seq = None
			gene = utils.collapse_mapped_gene(dict_2d[cb][umi])
			assert len(gene.keys()) == 1
			gene_mapped = str(list(gene.keys())[0])
			result = str(cb) + "\t" + str(umi) + "\t" + str(len(dict_2d[cb].keys())) + "\t" + str(len(dict_2d[cb][umi]))  + "\t" + region_ref + "\t"+ str(union_seq) + "\t" + str(union_nreads) + "\t" + str(is_WT) + "\t" + str(gene_mapped) + "\t" + str(union_positions) + "\t" + str(umi_seq_different)+ "\t" + str(ref_gRNA_id) + "\n"
		
			out.write(result)
