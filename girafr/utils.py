#!/usr/bin/env python

import pysam
import os
import subprocess
import re
from Bio import SeqIO
import pandas as pd
import configparser
from sklearn.mixture import GaussianMixture
import numpy as np

def create_bam_infile(file_name):
	bam_file = pysam.AlignmentFile(file_name, 'rb')
	if bam_file.check_index():
		return bam_file
	else:
		return 'Need to create index file first'


def get_read_tag(read, tag):
    try:
        r = read.get_tag(tag)
        if r == '':
            r = None
        return r
    except KeyError: # tag not found return None
        return None


def collapse_umi(reads):
	"""
	Arg:
		dict_2d[cb][umi]=[(seq, gene, cigar, is_WT, positions_overlapped)] for detect editing effect
		dict_2d[cb][umi]=[(seq, gene, cigar, is_WT)] for gRNA mutation detection	
	Return:
		dictionary: out[num, is_WT] for gRNA mutation detection
		dictionary: out[num, is_WT, positions_overlapped] for detect editing effect
	"""
	out ={}
	if(len(reads) == 1):
		seq = reads[0][0] # first record's sequence]
		if len(reads[0]) == 5: # dict_2d[cb][umi]=[(seq, gene, cigar, is_WT, positions_overlapped)]
			out[seq] = [1, reads[0][3], reads[0][4]]
		else: # dict_2d[cb][umi]=[(seq, gene, cigar, is_WT)]
			out[seq] = [1, reads[0][3]]
		return out
	else:
		for read in reads:
			seq = read[0]
			is_WT = read[3]
			if seq in out.keys():
				try:
					assert is_WT == out[seq][1]
				except AssertionError:
					print("AssertionError\n")
					print(out)
					print(read)
					exit()
				out[seq][0] += 1
			else:
				if len(read) == 5: # dict_2d[cb][umi]=[(seq, gene, cigar, is_WT, positions_overlapped)]
					positions_overlapped = read[4]
					out[seq] = [1, is_WT, positions_overlapped]
				else: # dict_2d[cb][umi]=[(seq, gene, cigar, is_WT)]
					out[seq] = [1, is_WT]

		return out


def collapse_umi_v2(reads, gene):
	"""
	Only count reads that aligned to specified gene
	if specified gene is None, take all reads. 
	if specified gene is 'multiple genes' (no consident mapped genes), take all reads 
	Called by generate_consensus_sequence_gRNA
	Arg:
		dict_2d[cb][umi]=[(seq, gene, cigar, is_WT)] for gRNA mutation detection	
		gene: specifed mapped gene, string, output of utils.which_mapped_gene
	Return:
		dictionary: out[num, is_WT] for gRNA mutation detection
	"""
	out ={}
	if(len(reads) == 1):
		seq = reads[0][0] # first record's sequence]
		out[seq] = [1, reads[0][3]]
		return out
	else:
		for read in reads:
			seq = read[0]
			mapped_gene = str(read[1]) # convert None to 'None'
			if mapped_gene != gene and gene != 'multiple genes': # read mapped to different gene, not count
				continue
			is_WT = read[3]
			if seq in out.keys():
				try:
					assert is_WT == out[seq][1]
				except AssertionError:
					print("AssertionError\n")
					print(out)
					print(read)
					exit()
				out[seq][0] += 1
			else:
				out[seq] = [1, is_WT]

		return out

def classify_umi(umi_seq):
	if(len(umi_seq.keys()) == 1):
		return 'single'
	else:
		return 'multiple'

def collapse_mapped_gene(reads):
	#dict_2d[cb][umi]=[(seq, gene, cigar)]
	dict={}
	if(len(reads) == 1):
		gene = reads[0][1]
		dict[gene] = 1
		return dict
	else:
		for read in reads:
			gene = read[1]
			if gene in dict.keys():
				dict[gene]+=1
			else:
				dict[gene]=1
		return dict

def which_mapped_gene(gene_dict):
	"""
	Decide which mapped gene when UMI in CB with multiple aligned gene results.
	The maximum reads supported results as mapped gene
	When tie, label as multiple genes as the gene cannot be confidently assigned
	Arg:
		gene_dict: gene name as key, value: the number of reads mapped to this gene. return of utils.collapse_mapped_gene()
	Return:
		gene_mapped
	"""
	if len(gene_dict.keys()) == 1:
		gene_mapped = str(list(gene_dict.keys())[0])
		return gene_mapped
	else:
		#gene_mapped = 'multiple genes' # old version
		gene_mapped = max(gene_dict, key = gene_dict.get) # maximum reads supported 
		num_reads = gene_dict[gene_mapped]
		gene_dict[gene_mapped] = 0 
		if num_reads in gene_dict.values(): # tie 
			return 'multiple genes'
		else:
			return str(gene_mapped)	

#def filter_mapping_record():	


def define_WT(record, is_10x):
	r = record
	is_WT = bool
	q_seq = r.query_sequence
	cigar = r.cigartuples # 43M1D48M => [(0, 43), (2, 1), (0, 48)]
	if len(cigar) == 1:
		cigar_WT = True # wildtype
	elif set([x[0] for x in cigar]) == set([0,4,5]) or set([x[0] for x in cigar]) == set([0,4]) or set([x[0] for x in cigar]) == set([0,5]): # consider soft cliping and hard cliping
		cigar_WT = True # cigar_WT is True is not necessary mean no mistmatch, because mismatches are marked as M in bam file. But mismatches will be counted in nM
	else:
		cigar_WT = False # mutant

	if (is_10x):
		nM = r.get_tag('nM') # 10x bam file does not contain MD, therefore, cannot use get_reference_sequence()
		try:
			if (nM == 0) and (cigar_WT): # if no mismtach and cigar_WT = True
				is_WT = True
			else:
				is_WT = False
		except TypeError:
			print('TypeError\n')
			print(r)
	else:
		r_seq = r.get_reference_sequence()
		if (q_seq == r_seq) and (cigar_WT):
			is_WT = True # perfect match
		else:
			is_WT = False
	return is_WT



def compare_two_sequences(seq1, seq2):
	# mismatch = utils.compare_two_sequences(ref[pos_r: (pos_r+n)], seq[pos_q: (pos_q+n)])
	# only compare two sequences which have same length
	# seq1= 'ATGCatcc' # ref
	# seq2= 'AAGCaTgc' # query
	# return 1A4C
	assert len(seq1) == len(seq2)
	seq1 = seq1.upper() # reference sequence normally
	seq2 = seq2.upper()
	pos = 0
	out = ''
	for i in range(0,len(seq1)):
		if seq1[i] ==  seq2[i]:
			pos += 1
		elif seq1[i] == 'N': # When N in reference sequence, report as match no matter for any seq2[i]. 
			pos += 1
		else:
			#if pos == 0:
			#	out += str(seq2[i])
			#else:
			#	out += str(pos)+str(seq2[i])
			out += str(pos) + str(seq2[i])
			pos = 0 # initalise counting
	out += str(pos) + 'M' # maybe should use a new symple to represent match #TODO
	return out


def which_structure(structure_dict, gene, pos):
	if gene not in structure_dict.keys(): 
		print('query is not in structure_dict dictionary')
		exit(1)
	out = ''
	for s in structure_dict[gene].keys():
		if s == 'gene': continue
		start = int(structure_dict[gene][s][1])
		end = int(structure_dict[gene][s][2])
		if (start <= pos) & (pos <= end):
			if out != '': 
				print('Overlapped structure')
				exit(1)
			else:
				out = s
	return out

def gRNA_bam_filter(input_file, samtools, output_dir):
	"""Filter gRNA library bam file
	Arg:
		input_file: gRNA library bam file mapped to custom reference with oligo pool as artificial chromosome
		samtools: path to samtools
	Return:
		Filtered gRNA lirbary bam file removed secondary alignment and mapped results not on gRNA reference
		
	"""
	subprocess.call('mkdir -p ' + output_dir + '/_tmp', shell = True)
	tmp_dir = output_dir + '/_tmp'
	subprocess.call('%s sort %s -o %s/gRNA.sorted.bam' % (samtools, input_file, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.bam' % (samtools, tmp_dir), shell = True)
	subprocess.call('%s view -H %s/gRNA.sorted.bam > %s/header' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -b -F 4 %s/gRNA.sorted.bam > %s/gRNA.sorted.mapped.bam' % (samtools, tmp_dir, tmp_dir), shell = True)

	subprocess.call('%s view -b -F 256 %s/gRNA.sorted.mapped.bam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call("%s view %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam|grep 'chrom' > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam" % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call("cat %s/header %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam > %s/tmp.sam" % (tmp_dir, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -S -b %s/tmp.sam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, tmp_dir, output_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, output_dir), shell = True)
	subprocess.call('rm -r ' + tmp_dir, shell = True)

def gRNA_bam_filter_v2(input_file, ref_fasta, samtools, output_dir):

	subprocess.call('mkdir -p ' + output_dir + '/_tmp', shell = True)
	tmp_dir = output_dir + '/_tmp'
	subprocess.call('%s sort %s -o %s/gRNA.sorted.bam' % (samtools, input_file, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.bam' % (samtools, tmp_dir), shell = True)
	subprocess.call('%s view -H %s/gRNA.sorted.bam > %s/header' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -b -F 4 %s/gRNA.sorted.bam > %s/gRNA.sorted.mapped.bam' % (samtools, tmp_dir, tmp_dir), shell = True)

	subprocess.call('%s view -b -F 256 %s/gRNA.sorted.mapped.bam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir), shell = True)
	fasta_sequences = SeqIO.parse(open(ref_fasta),'fasta')
	for fasta in fasta_sequences:
		oligo = fasta.id
		subprocess.call("%s view %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam %s >> %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam" % (samtools, tmp_dir, oligo, tmp_dir), shell = True)
	subprocess.call("cat %s/header %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam > %s/tmp.sam" % (tmp_dir, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -S -b %s/tmp.sam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, tmp_dir, output_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, output_dir), shell = True)
	subprocess.call('rm -r ' + tmp_dir, shell = True)


def gRNA_bam_filter_v3(input_file, ref_fasta, samtools, output_dir):

	subprocess.call('mkdir -p ' + output_dir + '/_tmp', shell = True)
	tmp_dir = output_dir + '/_tmp'
	subprocess.call('%s sort %s -o %s/gRNA.sorted.bam' % (samtools, input_file, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.bam' % (samtools, tmp_dir), shell = True)
	#subprocess.call('%s view -H %s/gRNA.sorted.bam > %s/header' % (samtools, tmp_dir, tmp_dir), shell = True)
	#subprocess.call('%s view -b -F 4 %s/gRNA.sorted.bam > %s/gRNA.sorted.mapped.bam' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -b -F 4 -F 256 %s/gRNA.sorted.bam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir, tmp_dir), shell = True)

	#subprocess.call('%s view -b -F 256 %s/gRNA.sorted.mapped.bam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir), shell = True)
	fasta_sequences = SeqIO.parse(open(ref_fasta),'fasta')
	cmd = samtools + ' merge -f ' + output_dir + '/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam'
	for fasta in fasta_sequences:
		oligo = fasta.id
		subprocess.call("%s view %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam %s -b >> %s/%s.bam" % (samtools, tmp_dir, oligo, tmp_dir, oligo), shell = True)
		cmd = cmd + ' ' + tmp_dir + '/' + oligo + '.bam'

	subprocess.call(cmd, shell = True)
	#subprocess.call("cat %s/header %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam > %s/tmp.sam" % (tmp_dir, tmp_dir, tmp_dir), shell = True)
	#subprocess.call('%s view -S -b %s/tmp.sam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, tmp_dir, output_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, output_dir), shell = True)
	subprocess.call('rm -r ' + tmp_dir, shell = True)

def detect_PAM(input_file):
	"""Detect if gRNA_region_coordinates_ori.txt contain protospacer sequence with PAM
	input_file: gRNA_region_coordinates_ori.txt
	
	return: boolean is_PAM
		is_PAM = True if all protosapcer sequence ended with PAM (5'-NGG-3' for Cas9) else is_PAM = False
	
	is_PAM is used later for cutsite determination
	"""
	is_PAM = bool 
	with open(input_file, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			seq = ls[5]
			if seq[-2:] == 'GG': # last 2 nt == GG
				is_PAM = True
			elif is_PAM == True: # if last protospacer end with GG, send a warning to user. 
				print("WARNING: not all protospacer sequence given end with PAM or without PAM")
				print("Please make sure all sgRNAs given in gRNA_region_coordinates_ori.txt end with PAM or all without PAM")
				print("Continue as all ithout PAM")
				is_PAM = False	
				return is_PAM
			else: # if all sgRNAs are without PAM
				is_PAM = False

	return is_PAM
	

def detect_editing_effect_input(input_file, out_dir, window_size = 31):
	"""prepare input file for detecing editing effect

        Arg:
                input_file: gRNA_region_coordinates_ori.txt

	chr14   77674163        77674186        Region2 -       ACGGCCATGTTTATGCACAGTGG ALKBH1
	chr14   77675744        77675767        Region3 -       AGGATTTCCGAGCTGAAGCAGGG ALKBH1
	chr14   77674113        77674136        Region4 +       GGACTGCGTGGTTCAAGAGGCGG ALKBH1
	chr14   77675696        77675719        Region5 +       GTCTACGTGGATTCCCAGTGTGG ALKBH1

        Returns:
		output_file: gRNA_region_coordinates_31bp.bed, gRNA_region_coordinates_31bp.txt
                Generate coordinate file by request detection window size (31bp by default)
	"""

	fo_name1 = out_dir + '/' + 'gRNA_region_coordinates_'+ str(window_size) + 'bp.bed'
	fo1 = open(fo_name1, 'w')
	fo_name2 = out_dir + '/' +'gRNA_region_coordinates_' + str(window_size) + 'bp.txt'
	fo2 = open(fo_name2, 'w')
	fo2.write('#chrom\tstart\tend\tregion_name\tstrand\tgRNA_seq\ttarget_gene\twindow_start\twindow_end\tcut_site\n')
	
	is_PAM = detect_PAM(input_file) # if PAM already there, return True. 
	if is_PAM:
		print('PAM at the end of protospacer sequence in input file: gRNA_region_coordinates_ori.txt')
	else:
		print('no PAM at the end of protospacer sequence in input file: gRNA_region_coordinates_ori.txt')

	window_size = int(window_size)
	with open(input_file, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			chrom = ls[0]
			start = int(ls[1])
			end = int(ls[2])
			name = ls[3]
			strand = ls[4]
			seq = ls[5]
			target_gene = ls[6]

			if strand == '+':
				if is_PAM: # if protospacer sequence includes PAM as the end
					cut_site =  end - 6
				else:
					cut_site = end - 3
			elif strand == '-':
				if is_PAM:
					cut_site = start + 6 - 1
				else:
					cut_site = start + 3 - 1
			else:
				print('invalid strand')
				exit()
			start2 = cut_site - round((window_size - 1)/2) # For 31bp window, minus 10
			end2 = cut_site + round((window_size - 1)/2)
			pf1 = '%s\t%s\t%s\t%s\n' % (chrom, start2-1, end2, name)
			fo1.write(pf1) # bed file: start2 - 1
			
			pf2 = '%s\t%s\t%s\t%s\n' % (l.strip(), start2, end2, cut_site)
			fo2.write(pf2)
	print("Generate files:\n" + fo_name1)
	print(fo_name2)
	return(fo_name1, fo_name2)


def read_barcode(barcode_file):

	"""Process barcode file
	Arg:
		barcode_file: filtered barcode list
	Return:
		barcodes: filtered barcode in list
	"""
	barcodes = []
	with open(barcode_file,'r') as fi:
		for barcode in fi:
			barcodes.append(barcode.strip().replace('-1',''))
	fi.close()
	return barcodes

def read_cell_assignment(cell_file):
	
	"""Process cell gRNA assignmet file
	Arg:
		cell_file: file that contain cell - sgRNA annotation
		For cellRanger feature barcoding analysis, cell_file is protospacer_calls_per_cell.csv
		or cells.gRNA.txt generated by previous step 
		or file with 
		Column 1: cell barcode
		Column 3: sgRNA asisgned (new row if multiple sgRNA/cell)
		
	Return:
		KO: dictionary KO[cb] = [sgRNA]
	"""
	KO = {}
	fi = open(cell_file, 'r')
	for l in fi:
		ls = l.strip().split(',')
		if len(ls) == 1:
			ls = l.strip().split()
		assert len(ls) > 1
		cb = ls[0].replace('-1','')
		gRNA = ls[2]
		if cb in KO.keys():
			KO[cb].append(gRNA)
		else:
			KO[cb]=[gRNA]
	fi.close()
	
	return KO

def write_annotation(oligo_pool, output_fasta, output_gtf, output_gtf2):
	
	"""Generate FASTA file and GTF to build cellranger reference. Modified from https://github.com/epigen/crop-seq
	Datlinger, P., Rendeiro, A. F., Schmidl, C., Krausgruber, T., Traxler, P., Klughammer, J., Schuster, L. C., Kuchler, A., Alpar, D., & Bock, C. (2017). Pooled CRISPR screening with single-cell transcriptome readout. Nature Methods, 14(3), 297301. https://doi.org/10.1038/nmeth.4177

	Arg:
		annotation.csv: oligo_name, sequence no header
			ALKBH1_gRNA1,TGGAAAGGACGAAACACCGGATTTCCGAGCTGAAGCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC
			ALKBH1_gRNA2,TGGAAAGGACGAAACACCGACTGCGTGGTTCAAGAGGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC

	Return:
		oligo_pool_plasmid.fa
		oligo_pool_plasmid.gtf
		oligo_pool_plasmid_structure.gtf

	"""
	df = pd.read_csv(oligo_pool, header=None, names=['oligo_name','sequence'])
	fasta_entries = list()
	gtf_entries = list()
	gtf2_entries = list()

	fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF"
	gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""
	gtf2_template = """{chrom}_chrom\thavana\t{type}\t{start}\t{length}\t.\t+\t.\tgene_name "{id}_gene";
"""
	U6_seq = "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTG" 
	Rest_seq = "TAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGACACTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGT"

	for i in df.index:
		oligo_name = df.iloc[i]["oligo_name"]
		guide_sequence = df.iloc[i]["sequence"]

		# add fasta entry
		sequence = U6_seq + guide_sequence + Rest_seq
		U6_seq_length = len(U6_seq) + len("TGGAAAGGACGAAACACC")
		gRNA_start = len(U6_seq) + len("TGGAAAGGACGAAACACC") + 1
		gRNA_end = len(U6_seq+guide_sequence) - len("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC")
		Rest_start = len(U6_seq+guide_sequence) - len("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC") + 1
		Rest_end = len(sequence)

		header = fasta_header_template.format(chrom=oligo_name, length=len(sequence))
		fasta_entries.append(header)
		fasta_entries.append(sequence)
		# add gtf entry
		gtf_entries.append(gtf_template.format(chrom=oligo_name, id=oligo_name, length=len(sequence)))
		gtf2_entries.append(gtf2_template.format(chrom=oligo_name, type='gene',id=oligo_name, start=1, length=len(sequence)))
		gtf2_entries.append(gtf2_template.format(chrom=oligo_name, type='U6',id=oligo_name, start=1, length=U6_seq_length))
		gtf2_entries.append(gtf2_template.format(chrom=oligo_name, type='gRNA',id=oligo_name, start=gRNA_start, length=gRNA_end))
		gtf2_entries.append(gtf2_template.format(chrom=oligo_name, type='Rest',id=oligo_name, start=Rest_start, length=Rest_end))

	# write to file
	with open(output_fasta, "w") as fasta_handle:
		fasta_handle.writelines("\n".join(fasta_entries))
	with open(output_gtf, "w") as gtf_handle:
		gtf_handle.writelines(gtf_entries)
	with open(output_gtf2, "w") as gtf2_handle:
		gtf2_handle.writelines(gtf2_entries)
	

def read_annotation_structure(input_file = 'oligo_pool_plasmid_structure.gtf'):
	"""Read in oligo pool structure annotation generated from utils.write_annotation
	Arg:
		input_file: oligo_pool_plasmid_structure.gtf

	Return:
		structure = {} Dictionary: structure[gene_name] = {s.split()[2]: [chrom, int(s.split()[3]), int(s.split()[4])]}
	"""
	structure = {}
	with open(input_file, 'r') as structure_gtf:
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
	return structure

def read_gRNA_reference(ref_fasta = 'oligo_pool_plasmid.fa', structure_gtf = 'oligo_pool_plasmid_structure.gtf'):

	"""Oligo pool reference, since cellranger alignmment result does not contain MD:Z, we have to extract reference sequence from fasta and gtf file of oligo pool. 
	Change to STAR? #TODO
	Arg:
		ref_fasta: oligo_pool_plasmid.fa
		structure_gtf: oligo_pool_plasmid_structure.gtf
	Return:
		ref_dict={} Dictionary ref_dict[gRNA] = reference_sequence

	"""
	# dropseq will be easier, by using MD:Z: TODO
	genome = {}
	fasta_sequences = SeqIO.parse(open(ref_fasta),'fasta')
	for fasta in fasta_sequences:
		genome[fasta.id] = str(fasta.seq)
	ref_dict={}
	structure = read_annotation_structure(structure_gtf) # structure = {} Dictionary: structure[gene_name] = {s.split()[2]: [chrom, int(s.split()[3]), int(s.split()[4])]}
	for gene in structure.keys():
		chrom = structure[gene]['gene'][0]
		start = structure[gene]['gene'][1] # 1-based to 0-based
		end = structure[gene]['gene'][2] # 1-based
		sequence = genome[chrom][(start-1):end ]
		ref_dict[gene] = sequence
	del genome
	return ref_dict


def get_gRNA_mutation_config(config_file):
	parser = configparser.ConfigParser()
	#ConfigFile = os.getcwd() + '/ConfigFile'
	ConfigFile = config_file
	parser.read(ConfigFile)
	gRNA_bam_file = parser.get('config_gRNA_mutation', 'gRNA_bam_file')
	barcode = parser.get('config_gRNA_mutation', 'filtered_barcode')
	output_dir = parser.get('config_gRNA_mutation', 'output_dir')
	output_dir = output_dir + '/'
	try:
		n_consensus_reads_min = int(parser.get('config_gRNA_mutation', 'min_reads'))
	except configparser.NoOptionError:
		n_consensus_reads_min  = 1	

	try:
		min_umi = int(parser.get('config_gRNA_mutation', 'min_umi'))
	except configparser.NoOptionError:
		min_umi = 3 # set default min_umi = 3
	try:
		auto = parser.get('config_gRNA_mutation', 'auto')
		auto = auto.lower() == 'true'
	except configparser.NoOptionError:
		auto = False # set default auto as False

	try:
		pool = parser.get('config_gRNA_mutation', 'pool')
		pool = pool.lower() == 'true'

	except configparser.NoOptionError:
		pool = False # set default pool as False
	
	ref_fasta = parser.get('config_annotation', 'ref_fasta')
	structure_gtf = parser.get('config_annotation', 'structure_gtf')
	try:
		is_10x = parser.get('config_annotation', 'is_10x')
		is_10x = is_10x.lower() == 'true'
	except configparser.NoOptionError:
		is_10x = True # set default is_10x as True
	return(gRNA_bam_file, barcode, output_dir, n_consensus_reads_min, min_umi, auto, pool, ref_fasta, structure_gtf, is_10x)

def get_tools_config(config_file):
	parser = configparser.ConfigParser()
	#ConfigFile = os.getcwd() + '/ConfigFile'
	ConfigFile = config_file
	parser.read(ConfigFile)
	samtools = parser.get('config_tools', 'samtools')
	twoBitToFa = parser.get('config_tools', 'twoBitToFa')
	featureCounts = parser.get('config_tools', 'featureCounts')
	return(samtools, twoBitToFa, featureCounts)

def get_detect_editing_effect_config(config_file):
	parser = configparser.ConfigParser()
	#ConfigFile = os.getcwd() + '/ConfigFile'
	ConfigFile = config_file
	parser.read(ConfigFile)
	input_file = parser.get('config_detect_editing_effect', 'gRNA_region')
	bam_file = parser.get('config_detect_editing_effect', 'expression_bam_file')
	barcode = parser.get('config_detect_editing_effect', 'filtered_barcode')
	cell_file = parser.get('config_detect_editing_effect', 'cell_file')
	window_size = parser.get('config_detect_editing_effect', 'window_size')
	output_dir = parser.get('config_detect_editing_effect', 'output_dir')
	
	genome = parser.get('config_annotation', 'genome')
	genome_gtf = parser.get('config_annotation', 'genome_gtf')

	return(input_file, bam_file, barcode, cell_file, genome, genome_gtf, output_dir, window_size)


def read_consensus_sequence_txt(in_file):
	"""
	Process consensus sequence txt file and return dictionary cells: [cb][gRNA]: {gRNA_type: UMI count}
	Args:
		consensus.sequence.gRNA.variant.txt
	Returns:
		Dictionary cells, cells[list(cells.keys())[0]] = {'HAT1_gRNA3_gene': {'WT': 1}, 'HDAC2_gRNA3_gene': {'WT': 1}, 'CTRL00545_gene': {'CTRL00545_gene_variant_1': 1}, 'RPL18_gRNA1_gene': {'WT': 1}, 'SMARCD1_gRNA4_gene': {'WT': 1}}
		cb_all, all cell barcodes that consensus.sequence.gRNA.variant.txt has
		gRNA_all all types of gRNA including mutant type
	"""
	cells = {}
	gRNA_all = []
	with open(in_file, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			cb = ls[0]
			gRNA = ls[8]
			if gRNA == 'multiple':
				continue #TODO

			gRNA_type = ls[9] # ALKBH1_gRNA1_gene_variant_1:21M1D21M2D48M
			
			if gRNA_type == 'WildType':
				if gRNA not in gRNA_all:
					gRNA_all.append(gRNA)
			else: # gRNA_type is mutant, treat it as different gRNA
				gRNA_type = ':'.join(gRNA_type.split(':')[0:-1]) #gRNA_type = gRNA_type.split(':')[0] can have problem if gene name has ':'
				if gRNA_type not in gRNA_all:
					gRNA_all.append(gRNA_type)

			if cb not in cells.keys():
				if gRNA_type == 'WildType':
					cells[cb] = {gRNA: {'WT': 1}}
                                
				else:
					cells[cb] = {gRNA: {gRNA_type: 1}}
			else:
				if gRNA not in cells[cb].keys():
					if gRNA_type == 'WildType':
						cells[cb][gRNA] = {'WT': 1}
					else:
						cells[cb][gRNA] = {gRNA_type: 1}
				else:
					if gRNA_type == 'WildType':
						if 'WT' not in cells[cb][gRNA].keys():
							cells[cb][gRNA]['WT'] = 1
						else:
							cells[cb][gRNA]['WT'] += 1
					else:
						if gRNA_type not in cells[cb][gRNA].keys():
							cells[cb][gRNA][gRNA_type] = 1
						else:
							cells[cb][gRNA][gRNA_type] += 1
	fi.close()

	cb_all = list(cells.keys())
	return(cells, cb_all, gRNA_all)


def fit_gmm(log_umi):
	"""
		Fit gaussian mixture model using scikit learn gaussianmixture
		fit 2 components 1-D univariate gaussian mixture model
	Args:
		log_umi: np.arrary

	Return:
		gmm: gmm object 
	"""
	gmm = GaussianMixture(n_components=2, covariance_type='tied').fit(log_umi)
	
	return gmm

def lookup_barcodes():
	"""
	Lookup function to search cell barcodes in cell ranger translation lookup table
	Especailly for cell ranger feature barcoding technology
	When crispr capture library alignment does not correct cell barcode, then lookup in 3M-february-2018.txt.gz
	Args:
		Default lookup_table is 3M-february-2018.txt.gz from cellranger
	Return:
		barcodes: dictionary with GEX capture sequence barcode variant1 as value, and faeture barcoding capture sequence barcode variant2 as key
		barcodes[cb2] = cb1
	"""
	infile = os.path.dirname(__file__) + '/3M-february-2018.txt'
	barcodes = {}
	with open(infile, 'r') as fi:
		for l in fi:
			ls = l.strip().split()
			cb1 = ls[0]
			cb2 = ls[1]
			if cb2 not in barcodes.keys():
				barcodes[cb2] = cb1
			else:
				print("duplicated cell barcode in lookup table")
				print("Check",lookup_table)
				exit(1)
	return barcodes



