#!/usr/bin/python

import os
import utils
import pysam
import subprocess
import consensus_sequence
from datetime import datetime
import variant
import assign_gRNA

time_start = datetime.now()
print('gRNA mutation profiling:', str(time_start),'\n')


(samtools, twoBitToFa, featureCounts) = utils.get_tools_config()
(gRNA_bam_file, barcode, output_dir, n_consensus_reads_min, min_umi, auto, pool, ref_fasta, structure_gtf, is_10x) = utils.get_gRNA_mutation_config()

######## Filtering of mapped reads ########
#utils.gRNA_bam_filter(gRNA_bam_file, samtools, output_dir) # time consuming
#print('Prepare bam file. Cost time: ' + str(datetime.now() - time_start) + '\n' )

#time_start = datetime.now()
#utils.gRNA_bam_filter_v2(gRNA_bam_file, ref_fasta, samtools, output_dir) # time consuming
#print('Prepare bam file. Cost time: ' + str(datetime.now() - time_start) + '\n' )


time_start = datetime.now()
utils.gRNA_bam_filter_v3(gRNA_bam_file, ref_fasta, samtools, output_dir) # time consuming
print('Prepare bam file. Cost time: ' + str(datetime.now() - time_start) + '\n' )
