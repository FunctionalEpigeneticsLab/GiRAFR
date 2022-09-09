#!/usr/bin/python

import utils
import os
import subprocess
from datetime import datetime
import consensus_sequence
import detect_cut_site
from Bio import SeqIO

time_start = datetime.now()

print('Detecing eiting effect start: ' + str(time_start) + '\n')



#input_file = 'gRNA_region_coordinates_ori.txt'
#output_dir = '/staging/leuven/stg_00064/projects/cropseq/267genes/develop/detect_editing_effect/cleanup/61bp_detection_window'
#twoBitToFa = '/staging/leuven/stg_00064/tools/ucsc/twoBitToFa'
#genome = '/staging/leuven/stg_00064/projects/cropseq/267genes/develop/detect_editing_effect/hg38.2bit' 
#samtools = '/data/leuven/334/vsc33470/miniconda3/bin/samtools'
#bam_file = 'possorted_genome_bam.bam'
#featureCounts = '/data/leuven/334/vsc33470/miniconda3/envs/umi_tools/bin/featureCounts'
#gtf = '/staging/leuven/stg_00064/tools/reference/cellranger/Homo_sapiens.GRCh38.93.gtf'
#barcode = '../apply_deep2/barcodes.tsv'
#cell_file = '/staging/leuven/stg_00064/master_thesis_projects/crop_seq/deep/deep-2/deep/outs/crispr_analysis/protospacer_calls_per_cell.csv'
#window_size = 61

(input_file, bam_file, barcode, cell_file, genome, gtf, output_dir, window_size) = utils.get_detect_editing_effect_config()
(samtools, twoBitToFa, featureCounts) = utils.get_tools_config()


print('Input gRNA region coorinates file: ' + input_file)
print('Detection window size: ' + str(window_size) + 'bp')
print('Output directory: ' + str(output_dir))
subprocess.call('mkdir -p ' + output_dir, shell = True)

# generate gRNA_region_coordinates_31bp.bed

print(input_file, output_dir, window_size)
(bed_file, txt_file) = utils.detect_editing_effect_input(input_file, output_dir, window_size)


cmd = '%s %s -bed=%s stdout > %s/gRNA_region_sequence.fa' % (twoBitToFa, genome, bed_file, output_dir)
print(cmd)
subprocess.call(cmd, shell = True)

reference = {}
fasta_sequences = SeqIO.parse(open(output_dir + '/gRNA_region_sequence.fa'), 'fasta')
for fasta in fasta_sequences:
	reference[fasta.id] = str(fasta.seq)

fi1 = open(txt_file,'r')
fi1.readline() # skip the header
fo_name = txt_file.replace('.txt', '_sequence.txt')
fo = open(fo_name, 'w')
for l1 in fi1:
	l1 = l1.strip()
	region_name = l1.split()[3]
	region_seq = reference[region_name]
	fo.write(l1 + '\t' + region_seq + '\n')

fo.close()
fi1.close()

# step 1 generate alignment file for each region
subprocess.call('mkdir -p ' + output_dir + '/region_bams', shell = True)
#subprocess.call(samtools + ' index '+ bam_file, shell = True)

fi_name = fo_name
with open(fi_name, 'r') as fi:
	for l in fi:
		if l[0] == '#': continue # skip the header
		ls = l.strip().split()
		chrom = ls[0]
		start = ls[7] # window_start
		end = ls[8] # window_end
		region_name = ls[3] # region name (unique)
		region_ref = ls[10] # region sequence reference
		target_gene = ls[6] # target gene
		cut_site = ls[9] # cut site
		
		print("Processing\t"+ region_name + "\n")

		# check if chr match with bam file
		check_bam = utils.create_bam_infile(bam_file)
		if check_bam.header.get('SQ')[0]['SN'][0:3] == 'chr' and chrom[0:3] == 'chr':
			pass
		elif check_bam.header.get('SQ')[0]['SN'][0:3] == 'chr' and chrom[0:3] != 'chr':
			chrom = 'chr' + chrom
		elif check_bam.header.get('SQ')[0]['SN'][0:3] != 'chr' and chrom[0:3] == 'chr':
			chrom = chrom.replace('chr','')
		else:
			print('bam file chromesome and gRNA_region coordinate cannot match')
			exit()

		cmd1 = '%s view -b %s %s:%s-%s > %s/region_bams/%s.bam' % (samtools, bam_file, chrom, start, end, output_dir, region_name)
		subprocess.call(cmd1, shell = True)
		# step 2 use featureCounts to annotate gene
		cmd2 = '%s -a %s -o %s/region_bams/%s.gene_assigned -R BAM %s/region_bams/%s.bam -T 4 -g gene_name >/dev/null 2>&1' % (featureCounts, gtf, output_dir, region_name, output_dir, region_name) 
		subprocess.call(cmd2, shell = True)
		
		sub_bam_file = output_dir + '/region_bams/' + region_name + '.bam.featureCounts.bam'
		sub_bam_sorted_file = output_dir + '/region_bams/' + region_name + '.bam.featureCounts.sorted.bam'
		subprocess.call(samtools + ' sort ' + sub_bam_file + ' > ' + sub_bam_sorted_file, shell = True)
		subprocess.call(samtools + ' index '+ sub_bam_sorted_file, shell = True)

		consensus_sequence.generate_consensus_sequence(sub_bam_sorted_file, barcode, chrom, start, end, region_ref, target_gene, region_name)	
		
		consensus_file = output_dir + '/region_bams/' + region_name + '.consensus.sequence.txt'
		detect_cut_site.detect_editing_effect(consensus_file, cut_site, cell_file)
		subprocess.call('rm ' + sub_bam_file, shell = True)

fi.close()	
subprocess.call('cat %s/region_bams/*.editing_effect.cutsite.txt > %s/all.editing_effect.cutsite.txt' % (output_dir, output_dir), shell = True)

print('Detecing eiting effect finished! Cost time:  ' + str(datetime.now() - time_start) + '\n')
