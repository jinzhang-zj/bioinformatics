#!/usr/bin/python

"""A script to automatically evaluate the transcriptome assembly quality.

The assembly quality was evaluated with metrics of contiguity
and completeness.  Both of the metrics were defined in the paper:

Martin JA, Wang Z:Next-generation transcriptome assembly.Nat Rev Genet.2011,12:671-682.

The input is the transcriptome assembly file and the reference
database name (protein database).

Version 0.1

The output contains the blast results, and the contiguity/completeness.

Currently version is only focusing on the translated assembly files.

It is recommended to use a well annotated complete transcriptome data of 
a species such as Arabidopsis as reference since contiguity/completeness
evaluate the quality of the transcriptome assembly in the sense of a
"complete" transcriptome data.

The default cutoff for the blast is 1e-10, and the cutoff for the contiguity
and completeness is 0.8.
"""

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/03/21 23:19:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"


import sys
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastxCommandline

# the blast results will be in xml format since it's the most portable blast
# results for different version of blast.

# generate the arg parser
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Evaluate transcriptome quality.')

if __name__ == "__main__":
	# add arguments
	parser.add_argument('-i', metavar='inputfile', nargs=1, type=argparse.FileType('r'), required=True, help='input assembly file')
	parser.add_argument('-db', metavar='refdatabase', nargs=1, required=True, help='reference database path/name')
	parser.add_argument('-p', metavar='projectname', nargs=1, required=True, help='project name used to name the output files')

	args = parser.parse_args()	# take sys.args as default

	# dictionary that stores all the arguments
	arglist = vars(args)

# running blast with transcriptome data against reference database

