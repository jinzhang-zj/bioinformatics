#!/usr/bin/python

"""Calculated given genome repeat content/number/distribution.

The input is given a genome sequence in fasta format. If the genome contains IR,
it must be removed.

The output is repeat content in the format of percentage of the whole genome.
Version 0.2

removing the genome mask step, which take O(N), N is the total number of nucleotides find all repeats
now the calculation of repeat content take O(n log n), n is the number repeats

some definitions:
	repeat_content:		percentage of nucleotides in repeat regions (i.e. 0.2)
	repeat_segs:		disjoint repeat regions, shown with index in forward strand
	repeat_pairs:		unique repeat pairs between forward and reverse strand, it has four indices
					- forward start
					- forward end
					- reverse start
					- reverse end
Prerequisite: local blast commands
"""

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/10/23 13:10:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

import argparse
import matplotlib.pyplot as plt
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

class repeat:
	def __init__ (self, i, j):
		self.start = int(i);
		self.end   = int(j);
		self.leng  = int(j) - int(i) + 1;
	
	def __lt__ (self, other):
		if self.start == other.start:
			return self.end < other.end
		else:
			return self.start < other.start
	
	def __str__ (self):
		return str(self.start) + "-" + str(self.end)

def merge_interval (reps):
	reps = sorted(reps);

	s = e = 0
	first = True
	result = []
	
	for rep in reps:
		if first:
			first = False
			s = rep.start
			e = rep.end
			continue
		else:
			if e < rep.start:
				result.append(repeat(s, e) )
				s = rep.start
				e = rep.end
			else:
				e = max(e, rep.end)

	result.append(repeat(s, e) )
	return result


def calc_repeat (infile,dbname,num):
	species = infile.replace(".fasta","")

	# running blast with transcriptome data against reference database, report the result in the xml
	# the blast results will be in xml format since it's the most portable blast results for different version of blast.
	blastncline = NcbiblastnCommandline(query=infile, db=dbname, num_threads=num, word_size=16, outfmt=5, out=species+'_blastn.out')
	stdout,stderr =  blastncline()

	# parsing the blast results
	result_handle = open(species+'_blastn.out')
	blast_records = NCBIXML.parse(result_handle)

	repeat_segs     = []	# number and content of disjoint repeat region
	repeat_pairs    = 0	# number of unique repeat pairs, i.e. "1-10,80-90" (forward strand 1-10, and reverse strand 80-90)
	genome_length   = 0	# total genome length
	repeat_content  = 0	# percentage of nucleotides of the repeat region
	temp		= {}	# storing temporatory repeat pair information

	# calculate the repeat content using mask
	# since here we have only one query and one hit, so there is only one blast_record and one alignment
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				# skip the first record which is the genome blast against itself, 
				if not genome_length:
					genome_length = blast_record.query_length
					continue
				else:
					
					# identify the record by query start/end and hit start/end coordinates
					repeat_segs.append(repeat(hsp.query_start, hsp.query_end) )
			
					lhs = str(hsp.query_start) + "-" + str(hsp.query_end)

					if hsp.sbjct_start < hsp.sbjct_end:
						rhs = str(hsp.sbjct_start) + "-" + str(hsp.sbjct_end)
					else:
						rhs = str(hsp.sbjct_end) + "-" + str(hsp.sbjct_start)
						
					if lhs < rhs:
						idx = rhs+','+lhs
					else:
						idx = lhs+','+rhs
					
					if idx not in temp:
						temp[idx] = 1
						repeat_pairs += 1
					else:
						pass

	repeat_segs = merge_interval(repeat_segs)
	
	repeat_content = sum(i.leng for i in repeat_segs) * 1.0 / genome_length
	return repeat_content, repeat_pairs, repeat_segs
	


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
	print "number of repeat pairs: ", result[1]
	print "number of disjoint repeat segments: ", len(result[2])
	print "disjoint repeat segments: "

	#plt.hist(result[2],bins=10)
	#plt.savefig(species+"_repeats.png")


