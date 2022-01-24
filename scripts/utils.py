#!/usr/bin/python

import pysam
import os
import subprocess
import re
from Bio import SeqIO
import pandas as pd
import configparser

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
    except KeyError:
        return None


def collapse_umi(reads):
	#dict_2d[cb][umi]=[(seq, gene, cigar, is_WT, positions_overlapped)]
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
#def filter_mapping_record():	


def define_WT(record, is_10x):
	r = record
	is_WT = bool
	q_seq = r.query_sequence
	cigar = r.cigartuples # 43M1D48M => [(0, 43), (2, 1), (0, 48)]
	if (is_10x):
		nM = r.get_tag('nM') # 10x bam file does not contain MD, therefore, cannot use get_reference_sequence()
		try:
			if (nM == 0)  & (len(cigar) == 1):
				is_WT = True
			else:
				is_WT = False
		except TypeError:
			print('TypeError\n')
			print(r)
	else:
		r_seq = r.get_reference_sequence()
		if (q_seq == r_seq) & (len(cigar) == 1):
			is_WT = True # perfect match
		else:
			is_WT = False
	return is_WT



def compare_two_sequences(seq1, seq2):
	# only compare two sequences which have same length
	# seq1= 'ATGCatcc'
	# seq2= 'AAGCaTgc'
	# return 1A4C
	assert len(seq1) == len(seq2)
	seq1 = seq1.upper()
	seq2 = seq2.upper()
	pos = 0
	out = ''
	for i in range(0,len(seq1)):
		if seq1[i] ==  seq2[i]:
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
	subprocess.call('%s sort %s > %s/gRNA.sorted.bam' % (samtools, input_file, tmp_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.bam' % (samtools, tmp_dir), shell = True)
	subprocess.call('%s view -H %s/gRNA.sorted.bam > %s/header' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -b -F 4 %s/gRNA.sorted.bam > %s/gRNA.sorted.mapped.bam' % (samtools, tmp_dir, tmp_dir), shell = True)

	subprocess.call('%s view -b -F 256 %s/gRNA.sorted.mapped.bam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam' % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call("%s view %s/gRNA.sorted.mapped.removedSecondaryAlignment.bam|grep 'chrom' > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam" % (samtools, tmp_dir, tmp_dir), shell = True)
	subprocess.call("cat %s/header %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.sam > %s/tmp.sam" % (tmp_dir, tmp_dir, tmp_dir), shell = True)
	subprocess.call('%s view -S -b %s/tmp.sam > %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, tmp_dir, output_dir), shell = True)
	subprocess.call('%s index %s/gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam' % (samtools, output_dir), shell = True)
	subprocess.call('rm -r ' + tmp_dir, shell = True)
	



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
				cut_site =  end - 6
			elif strand == '-':
				cut_site = start + 6 - 1
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


def get_gRNA_mutation_config():
	parser = configparser.ConfigParser()
	ConfigFile = os.getcwd() + '/ConfigFile'
	parser.read(ConfigFile)
	gRNA_bam_file = parser.get('config_gRNA_mutation', 'gRNA_bam_file')
	barcode = parser.get('config_gRNA_mutation', 'filtered_barcode')
	output_dir = parser.get('config_gRNA_mutation', 'output_dir')
	min_umi = int(parser.get('config_gRNA_mutation', 'min_umi'))
	
	ref_fasta = parser.get('config_annotation', 'ref_fasta')
	structure_gtf = parser.get('config_annotation', 'structure_gtf')
	return(gRNA_bam_file, barcode, output_dir, min_umi, ref_fasta, structure_gtf)

def get_tools_config():
	parser = configparser.ConfigParser()
	ConfigFile = os.getcwd() + '/ConfigFile'
	parser.read(ConfigFile)
	samtools = parser.get('config_tools', 'samtools')
	twoBitToFa = parser.get('config_tools', 'twoBitToFa')
	featureCounts = parser.get('config_tools', 'featureCounts')
	return(samtools, twoBitToFa, featureCounts)

def get_detect_editing_effect_config():
	parser = configparser.ConfigParser()
	ConfigFile = os.getcwd() + '/ConfigFile'
	parser.read(ConfigFile)
	input_file = parser.get('config_annotation', 'gRNA_region')
	bam_file = parser.get('config_detect_editing_effect', 'expression_bam_file')
	barcode = parser.get('config_detect_editing_effect', 'filtered_barcode')
	cell_file = parser.get('config_detect_editing_effect', 'cell_file')
	window_size = parser.get('config_detect_editing_effect', 'window_size')
	output_dir = parser.get('config_detect_editing_effect', 'output_dir')
	
	genome = parser.get('config_annotation', 'genome')
	genome_gtf = parser.get('config_annotation', 'genome_gtf')

	return(input_file, bam_file, barcode, cell_file, genome, genome_gtf, output_dir, window_size)
