#!/usr/bin/env python

import sys


def usage():
	print >> sys.stderr, '''usage: %s <genemap> <column file> 

where a <column file> is any file which has a <chr>[TAB]<index>[TAB]<etc...> layout''' % sys.argv[0]
	exit(-1)


if len(sys.argv) < 3:
	usage()


## Read in Map file into a dictionary of dictionaries
######################################################################
def maplibrary(mapfile, intergenic=True):
	try:
		m = open(mapfile,'r')
	except IOError:
		print >> sys.stderr, "Could not open map:", mapfile
		exit(-1);

	## How it works:
	#
	#	Retrieve chr ----> Array of [50 000 000, 100 000 000 .... , 250 000 000 ]  (50 elements)
	#							 index =  position%(50 000 000),
	#
						# ---> Further array of [ 100 000, 200 000, ..., 


	chromemap = {}

	for line in m:
		tokens = line.split('\t');

		chrom = tokens[0].strip()
		pos1 = int(tokens[1].strip()); 
		pos2 = int(tokens[2].strip())
		gene = tokens[3].split('_')[0].strip()

		if chrom not in chromemap:
			genemap = {}
			chromemap[chrom] = genemap;

		else:
			# Chrom exists, check genemap
			if gene not in genemap:
				chromemap[chrom][gene] = [pos1, pos2]
			else:
				one,two = chromemap[chrom][gene]
				if pos1 < one:
					one = pos1
				if pos2 > two:
					two = pos2
				chromemap[chrom][gene] = [one, two]


	m.close()
	return chromemap

## Read in VCF file and append appropriate gene names
#########################################################################
def appendtoVCF(colfile, chromemap):

	try:
		c = open(colfile,'r')
	except IOError:
		print >> sys.stderr, "Could not open column file:", colfile
		exit(-1);


	count = 0

	line = c.readline()
	while (line[0]=='#'):
		count +=1
		line = c.readline()


	for line in c:
		count +=1

		print >> sys.stderr, "\rLine: ", count,

		try:
			tokens = line.split('\t')
		except IndexError:
			print >> sys.stderr, line
			continue

		chrom = tokens[0].strip()
		pos = tokens[1].strip()

		for gene in chromemap[chrom].keys():
			pos1, pos2 = chromemap[chrom][gene]

			if (int(pos1) <= int(pos) <= int(pos2)):
				print line.splitlines()[0], '\t', gene
				break


	c.close()


### MAIN ###

mapfile = sys.argv[1];
colfile = sys.argv[2];

chromemap = maplibrary(mapfile)
appendtoVCF(colfile, chromemap)

