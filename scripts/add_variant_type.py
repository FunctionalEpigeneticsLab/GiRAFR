#!/usr/bin/env python

metadata = open('/staging/leuven/stg_00064/projects/cropseq/267genes/develop/custom/apply_to_deep2/metadata.variant.txt','r')

outdir
in_file = open('consensus.sequence.gRNA.variant.txt', 'r')

variants ={}

for line in metadata:
	line = line.strip()
	lines = line.split()
	gRNA = lines[8]
	seq = lines[5]
	if gRNA == 'multiple': # TODO
		continue
	gRNA_type = lines[9]
	if gRNA_type != 'WildType':
		if len(lines) < 11: continue # TODO,for variant like 90M1S, soft  cliping here is skipped
		gRNA = gRNA_type.split(':')[0]
		gRNA_MT = lines[10]
	else:
		gRNA = gRNA + '_WT' # ALKBH1_gRNA1_gene_WT
		gRNA_MT = None
		seq = None # As different wiltype reads are not exact match, some might have a bit more sequence to the upstream or downstream
	if gRNA not in variants:
		variants[gRNA] = [gRNA_MT, seq]
	else:
		if variants[gRNA][0] != gRNA_MT or variants[gRNA][1] != seq:
			print('records do not match')
			exit()
		else:
			continue
			


cells_gRNA = open('/staging/leuven/stg_00064/projects/cropseq/267genes/develop/custom/apply_to_deep2/cells.gRNA.single.txt','r')
out = open('/staging/leuven/stg_00064/projects/cropseq/267genes/develop/custom/apply_to_deep2/cells.gRNA.single.MT.txt','w')

for line in cells_gRNA:
	line = line.strip()
	lines = line.split()
	gRNA = lines[2]
	if gRNA not in variants: #TODO, CHD1L_gRNA2_gene_variant_20  64M27S  
		out_line = line + '\t' + str(None) + '\t' + str(None)
	else:
		out_line = line + '\t' + str(variants[gRNA][0]) + '\t' +str(variants[gRNA][1])
	out.write(out_line + '\n')


