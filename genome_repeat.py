#!/usr/bin/python

"""Calculated given genome repeat content and repeat distribution(optional).

The input is given a genome sequence in fasta format. If the genome contains IR,
it must be removed.

The output is repeat content in the format of percentage of the whole genome.
Version 0.1

Repeat content calculation: using mask calculate the percentage of repeat lengths
on the whole genome.

Repeat number/distribution: Each pair of repeat is considerred as an unique group,
and the number/length distribution of these repeats were calculated. For example, if
we have three identical segments A, B, C, then the number of repeat groups are 3: AB,
AC, BC.

Prerequisite: local blast commands
"""

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/08/13 13:10:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

import argparse
import matplotlib.pyplot as plt
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


def calc_repeat(infile,dbname,num):
	species = infile.replace(".fasta","")

	# running blast with transcriptome data against reference database, report the result in the xml
	# the blast results will be in xml format since it's the most portable blast results for different version of blast.
	blastncline = NcbiblastnCommandline(query=infile, db=dbname, num_threads=num, word_size=16, outfmt=5, out=species+'_blastn.out')
	stdout,stderr =  blastncline()

	# parsing the blast results
	result_handle = open(species+'_blastn.out')
	blast_records = NCBIXML.parse(result_handle)

	repeat_groups = {}
	repeat_length = []
	genome_mask = [0]
	genome_length = 0

	# calculate the repeat content using mask
	# since here we have only one query and one hit, so there is only one blast_record and one alignment
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				# skip the first record which is the genome blast against itself, 
				if not genome_length:
					genome_length = blast_record.query_length
					genome_mask = [0] * genome_length
					continue
				else:
					# mask the genome based on query alignment
					genome_mask[hsp.query_start-1:hsp.query_end-1]  = [1] * (hsp.query_end - hsp.query_start + 1)

					# identify the record by query start/end and hit start/end coordinates
					lhs = str(hsp.query_start) + "-" + str(hsp.query_end)
					if hsp.sbjct_start < hsp.sbjct_end:
						rhs = str(hsp.sbjct_start) + "-" + str(hsp.sbjct_end)
					else:
						rhs = str(hsp.sbjct_end) + "-" + str(hsp.sbjct_start)
						
					if lhs < rhs:
						idx = rhs+','+lhs
					else:
						idx = lhs+','+rhs
					
					if idx in repeat_groups:
						continue

					repeat_groups[idx] = hsp.align_length
					repeat_length.append(hsp.query_end - hsp.query_start)

	repeat_content = sum(genome_mask)*1.0/genome_length
	return repeat_content,repeat_groups,repeat_length
	


if __name__ == "__main__":
	# generate the arg parser
	parser = argparse.ArgumentParser(
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Evaluate transcriptome quality.')

	# add arguments
	parser.add_argument('-i', metavar='genomefile', nargs=1, required=True, help='genome file')
	parser.add_argument('-db', metavar='reference', nargs=1, required=True, help='reference database name')
	parser.add_argument('-n', metavar='number of threads', nargs=1, default=12, help='number of threads you want to use, recommended to be the number of cores you can use ')
		
	args = parser.parse_args()	# take sys.args as default

	# dictionary that stores all the arguments
	arglist = vars(args)

	infile = arglist['i'][0]
	db = arglist['db'][0]
	species = infile.replace(".fasta","")


	result = calc_repeat(infile,db,arglist['n'])

	print "repeat content: ", result[0]
	print "number of repeat pairs: ", len(result[1])

	plt.hist(result[2],bins=10)
	plt.savefig(species+"_repeats.png")


