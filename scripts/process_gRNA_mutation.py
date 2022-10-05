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
print('Prepare bam file. Cost time: ' + str(datetime.now() - time_start) + '\n' )

######## gRNA consensus sequence ##########
print('Generating consensus sequence for gRNA library\n')
#bam_in_file = output_dir + 'gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam'
#consensus_sequence.generate_consensus_sequence_gRNA(bam_in_file, barcode, output_dir, n_consensus_reads_min, is_10x = is_10x)

#subprocess.call('%s index %s/consensus.bam' % (samtools, output_dir), shell = True)

####### Identification of mutations in the gRNA consensus ########
#variant.call_gRNA_variant(output_dir, output_dir + 'consensus.sequence.gRNA.txt', ref_fasta, structure_gtf, is_10x = is_10x)

####### Assign gRNAs to cells ###########
assign_gRNA.assign_gRNA_to_cell(in_file = output_dir + 'consensus.sequence.gRNA.variant.txt', min_umi = min_umi, output_dir = output_dir, auto = auto, pool = pool)

assign_gRNA.add_variant_type(in_file1 = output_dir + 'consensus.sequence.gRNA.variant.txt', in_file2 = output_dir + 'cells.gRNA.single.txt', in_file3 = output_dir + 'cells.gRNA.txt', structure_gtf = structure_gtf, output_dir = output_dir)

print('gRNA mutation finished! Cost time: ' + str(datetime.now() - time_start) + '\n' )



