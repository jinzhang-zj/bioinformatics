#!/opt/apps/python/epd/7.2.2/bin/python

'''Extract the reads mapped to the reference in given sam file

Take the sam files and extract the read pairs from which either
one or two reads get aligned to the reference file. Additional
opitions such as 'disconcordant' can be specified to also extract
the reads that map disconcordantly to the reference file.

More opitions will be added later...

Version 2.0
Modifications: 
Take one more argument which is the prefix of orginial read files.

Instead of direct extract all information from 
single sam file, the script here will try to get read name from
sam file and extract corresponding read sequences from given original
read files.

The algorithm of is based on the observation that reads in the sam
file are in the same order as that of the original read files.

Note: the script assumes the input sam file is the original one without
filtering, in which case the order/list of reads will be the same between
sam file and read file!
'''

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 2.0 $"
__date__ = "$Date: 2012/12/07 16:55:10"
__copyright__ = "Copyright: (c) 2012 Jin Zhang"
__license__ = "Python"


import sys
import argparse
from Bio.Seq import Seq

# generate the arg parser
parser = argparse.ArgumentParser(
	#usage='%(prog)s [options]',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='Extract reads from mapping results')

if __name__ == "__main__":
	# add arguments
	parser.add_argument('-s', metavar='samfile', nargs=1, type=argparse.FileType('r'), required=True, help='input sam file')
	parser.add_argument('-i', metavar='oripre', nargs=1, help='prefix of input read files')
	parser.add_argument('-p', metavar='etxpre', default='sam_read', nargs='?', help='prefix of output reads files')
	parser.add_argument('-disconcordant', action='store_true', help='including the disconcordantly mapped reads')
	parser.add_argument('-v','--version',action='version', version='%(prog)s 2.0')

	args = parser.parse_args()	# take sys.args as default

	# dictionary that stores all the arguments
	arglist = vars(args)

fwd_fq = []
bak_fq = []	# forward and backward fastq sequences

samfile = arglist['s'][0]

# input read
orifilel = open(arglist['i'][0] + '1.fq', 'r')
orifiler = open(arglist['i'][0] + '2.fq', 'r')

# output read
fwd_file = open(arglist['p'] + '1.fq', 'w')
bak_file = open(arglist['p'] + '2.fq', 'w')

# process read the same file
while True:
	line =  samfile.readline()
	
	if len(line) == 0:
		break

	if line[0] == '@':	# skip the sam header
		continue
		
	sline = samfile.readline()	# read in the read pair

	# read pair
	fread = line.split()
	bread = sline.split()

	# corresponding read seqs
	orireadl = ''
	orireadr = ''
	for i in range(4):	# read four lines each time
		orireadl += orifilel.readline()
		orireadr += orifiler.readline()

	# skip the unmapped reads
	if fread[1] == '77' or fread[1] == '141':
		continue
			
	# mapped concordantly
	if fread[1] == '99' or fread[1] == '147' or fread[1] == '83' or fread[1] == '163':
		fwd_file.write(orireadl)	
		bak_file.write(orireadr)	
	# mapped disconcordantly
	elif fread[1] == '65' or fread[1] == '129' or fread[1] == '113' or fread[1] == '177' or fread[1] == '81' or fread[1] == '161':
		if arglist['disconcordant']:			
			fwd_file.write(orireadl)	
			bak_file.write(orireadr)	
	# single read mapping
	elif fread[1] == '73' or fread[1] == '137' or fread[1] == '89' or fread[1] == '153' or fread[1] == '133' or fread[1] == '69':
		fwd_file.write(orireadl)	
		bak_file.write(orireadr)	
	else:
		print("Unexpected number" + str(fread[1]))
		continue
		
fwd_file.close()
bak_file.close()
