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

assign_gRNA.assign_gRNA_to_cell(in_file = output_dir + 'consensus.sequence.gRNA.variant.txt', min_umi = min_umi, output_dir = output_dir, auto = auto, pool = pool)
