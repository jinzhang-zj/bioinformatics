#!/usr/bin/python

"""A script to automatically evaluate the transcriptome assembly quality.

The assembly quality was evaluated with metrics of contiguity
and completeness.  Both of the metrics were defined in the paper:

Martin JA, Wang Z:Next-generation transcriptome assembly.Nat Rev Genet.2011,12:671-682.

The input is the transcriptome assembly file and the reference
database name (protein database).

Version 0.1

The output contains the blast results, and the contiguity/completeness.

Currently version is only for on the translated assembly files.

It is recommended to use a well annotated complete transcriptome data of 
a species such as Arabidopsis as reference since contiguity/completeness
evaluate the quality of the transcriptome assembly in the sense of a
"complete" transcriptome data.

The default cutoff for the blast is 1e-10, and the cutoff for the contiguity
and completeness is 0.8.

Prerequisite: local blast commands

"""

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/03/21 23:19:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"


import sys
import argparse
import os
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

# the blast results will be in xml format since it's the most portable blast
# results for different version of blast.

# generate the arg parser
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Evaluate transcriptome quality.')

if __name__ == "__main__":
	# add arguments
	parser.add_argument('-i', metavar='inputfile', nargs=1, required=True, help='input assembly file')
	parser.add_argument('-db', metavar='refdatabase', nargs=1, required=True, help='reference database path/name')
	parser.add_argument('-p', metavar='projectname', nargs=1, required=True, help='project name used to name the output files')
	parser.add_argument('-e', metavar='evalue', nargs=1, default=1e-10, help='blast evalue')
		
	args = parser.parse_args()	# take sys.args as default

	# dictionary that stores all the arguments
	arglist = vars(args)

# running blast with transcriptome data against reference database, report the result in the xml
blastpcline = NcbiblastpCommandline(query=arglist['i'][0], db=arglist['db'][0], max_target_seqs=1, evalue=arglist['e'], outfmt=5, out=arglist['p'][0]+'.blastp.out')
stdout,stderr =  blastpcline()


# parsing the blast results
result_handle = open(arglist['p'][0]+'.blastp.out')
blast_records = NCBIXML.parse(result_handle)

cont=defaultdict(int)
dbnum = 0;



for blast_record in blast_records:
	# skip the query if not hit was found
	

	if not blast_record.alignments:
		continue
	
	if not dbnum:
		dbnum = blast_record.database_sequences
		print dbnum

	cur_query = blast_record.query
	
	# each alignment corresponding to blast reuslts of one hit for the query
	# here we are only interested in the top hit, or the first alignment results

	alignment = blast_record.alignments[0]
	cur_hit = alignment.title[:20]
	cur_hit_length = alignment.length

	# each hsp corresponds to blast results of one hsp of one hit for one query
	# here we are only interested in the top hsp
	hsp = alignment.hsps[0]
	align_length = hsp.sbjct_end - hsp.sbjct_start + 1
	tmp_cont = align_length * 1.0 / cur_hit_length
	
	if tmp_cont > cont[cur_hit]:
		cont[cur_hit] = tmp_cont

# contiguity output
for cut in np.arange(0,1,0.1):
	count = sum([1 for i in cont.values() if i > cut]) * 1.0 / dbnum
	print str(cut) + ' ' + str(count)












# calculating contiguity

# calculating completeness
