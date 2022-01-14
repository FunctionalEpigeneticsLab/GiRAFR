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
(gRNA_bam_file, barcode, output_dir, min_umi, ref_fasta, structure_gtf) = utils.get_gRNA_mutation_config()

#utils.gRNA_bam_filter(gRNA_bam_file, samtools, output_dir) # time consuming
print('Prepare bam file. Cost time: ' + str(datetime.now() - time_start) + '\n' )

bam_in_file = output_dir + 'gRNA.sorted.mapped.removedSecondaryAlignment.onlyMappedToGrnaChrom.bam'
#consensus_sequence.generate_consensus_sequence_gRNA(bam_in_file, barcode, output_dir)

#subprocess.call('%s index %s/consensus.bam' % (samtools, output_dir), shell = True)

variant.call_gRNA_variant(output_dir, output_dir + 'consensus.sequence.gRNA.txt', ref_fasta, structure_gtf)

assign_gRNA.assign_gRNA_to_cell(in_file = output_dir + 'consensus.sequence.gRNA.variant.txt', min_umi = min_umi, output_dir = output_dir)

assign_gRNA.add_variant_type(in_file1 = output_dir + 'consensus.sequence.gRNA.variant.txt', in_file2 = output_dir + 'cells.gRNA.single.txt', structure_gtf = structure_gtf, output_dir = output_dir)

print('gRNA mutation finished! Cost time: ' + str(datetime.now() - time_start) + '\n' )



