#!/usr/bin/python

import utils
import os
import subprocess
from datetime import datetime
import consensus_sequence
import detect_cut_site
from Bio import SeqIO


sub_bam_sorted_file='/lustre1/project/stg_00064/projects/cropseq/267genes/develop/paper/exp2_deep/51bp_detection_window/region_bams/Region65.bam.featureCounts.sorted.bam'
barcode='/staging/leuven/stg_00064/projects/cropseq/267genes/develop/custom/apply_to_exp2_deep/barcodes.tsv'
chrom='chr19'
start=48615376
end=48615426
region_ref='TTTTTGTAGCCTCGGCTGGCCCGTCGGCCTCTGGCACGCTCGAACTTCCGG'
target_gene='RPL18'
region_name='Region65'
#consensus_sequence.generate_consensus_sequence(sub_bam_sorted_file, barcode, chrom, start, end, region_ref, target_gene, region_name)

consensus_file='/lustre1/project/stg_00064/projects/cropseq/267genes/develop/paper/exp2_deep/51bp_detection_window/region_bams/Region65.consensus.sequence.txt'
cut_site=48615401
cell_file='/staging/leuven/stg_00064/projects/cropseq/267genes/develop/paper/exp2_deep/auto_pool/cells.gRNA.txt'


detect_cut_site.detect_editing_effect(consensus_file, cut_site, cell_file)
