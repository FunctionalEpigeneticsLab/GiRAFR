#!/usr/bin/env python

with open('cells.gRNA.txt','r') as fi:
	num_cells = 0 
	num_MT = 0
	cb_ini = ''
	num_gRNA = 0
	for l in fi:
		ls = l.strip().split()
		cb = ls[0]
		if(cb != cb_ini):
			num_cells = num_cells + num_gRNA
			cb_ini = cb
			num_gRNA = int(ls[1])
		else:
			num_gRNA = int(ls[1])
		num_MT = num_MT + ((len(ls)-5)/2)
	num_cells = num_cells + num_gRNA
	print(num_cells)
	print(num_MT)
