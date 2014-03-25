#!/usr/bin/python

"""A script to automatically evaluate the transcriptome assembly quality.

The blast using multithreading and the following steps is calculated 
sequentially.

The assembly quality was evaluated with metrics of contiguity
and completeness.  Both of the metrics were defined in the paper:

Martin JA, Wang Z:Next-generation transcriptome assembly.Nat Rev Genet.2011,12:671-682.

The input is the transcriptome assembly file and the reference
database name (protein database).

Version 0.1

The output contains:
blast results,    project_blast[x|p].out
contiguity,	  project_cont.out
completeness,	  project_comp.out
log,	 	  project.log
summary,	  project.summary

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

import time
import argparse
import numpy as np
from collections import defaultdict
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline

start = time.time()
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

logfile = open(arglist['p'][0] + '.log','w')
contfile = open(arglist['p'][0] + '_cont.out','w')
compfile = open(arglist['p'][0] + '_comp.out','w')
sumfile = open(arglist['p'][0] + '.summary','w')

logfile.write("""This is script is used to evaluate contiguity/completeness with given transcriptome assembly.

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/03/21 23:19:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

""")

logfile.write("running blast...\n")
logfile.write("database: " + arglist['db'][0] + "\n")
# running blast with transcriptome data against reference database, report the result in the xml
blastpcline = NcbiblastpCommandline(query=arglist['i'][0], db=arglist['db'][0], max_target_seqs=1, evalue=arglist['e'], num_threads=16,outfmt=5, out=arglist['p'][0]+'_blastp.out')
stdout,stderr =  blastpcline()
logfile.write("done...\n")
logfile.write("database information: ")

# parsing the blast results
result_handle = open(arglist['p'][0]+'_blastp.out')
blast_records = NCBIXML.parse(result_handle)

cont = defaultdict(int)		# record contiguity for each gene
comp_seg = defaultdict(set)	# record aligned hsp start and end position pairs
ref = {}			# record reference length
cont_list = {}			# record the contig with highest aligned length to each reference gene
comp_list = defaultdict(list)	# record the contigs aligned to each reference


dbnum = 0;

for blast_record in blast_records:
	if not dbnum:
		dbnum = blast_record.database_sequences
		logfile.write(str(dbnum) + " reference sequences\n")
		logfile.write("\n" + "now parsing the blast results...\n")
	
	# skip the query if not hit was found
	if not blast_record.alignments:
		continue
	
	# each alignment corresponding to blast reuslts of one hit for the query
	# we are only interested in the top hit, or the first alignment results

	alignment = blast_record.alignments[0]
	cur_hit = " ".join(alignment.title.split()[1:])	# remove additional index formated by the database
	ref[cur_hit] = alignment.length

	# each hsp corresponds to blast results of one hsp of one hit for one query
	# for contiguity we are only interested in the top hsp
	# calculating contiguity for each contig
	hsp = alignment.hsps[0]
	align_length = hsp.sbjct_end - hsp.sbjct_start + 1
	tmp_cont = align_length * 1.0 / ref[cur_hit]
	if tmp_cont > cont[cur_hit]:
		cont[cur_hit] = tmp_cont
		cont_list[cur_hit] = blast_record.query

	# for contiguity we are interested in all hsps over the evalue cutoff
	# calculating completeness for each contig
	comp_list[cur_hit].append(blast_record.query)
	for hsp in alignment.hsps:
		if hsp.expect > arglist['e']:
			continue
		tmp_tuple = (hsp.sbjct_start,  hsp.sbjct_end) if (hsp.sbjct_start < hsp.sbjct_end) else (hsp.sbjct_end, hsp.sbjct_start)
		comp_seg[cur_hit].add(tmp_tuple)


logfile.write("parsing done.\n\n")
logfile.write("summarizinging contiguity...\n")

# summarize contiguity over all the references
cont_counts = [0]*10
for idx in range(10):
	cont_counts[idx] = sum (1 for v in cont.values() if v*10 > idx) * 1.0/ dbnum

logfile.write("done.\n\n")
logfile.write("summarizinging completeness...\n")

# summarize completeness over all the references
comp = defaultdict(int)

for gene in comp_seg.keys():
	# if only one hsp was found for the gene, completeness = contiguity
	if len(comp_seg[gene]) == 1:
		comp[gene] = cont[gene]
	else:
		marker = np.array([0]*ref[gene])
		for (left, right) in comp_seg[gene]:
			marker[left-1:right-1] = 1
		comp[gene] = sum(marker) * 1.0 / ref[gene]		

comp_counts = [0]*10
for idx in range(10):
	comp_counts[idx] = sum (1 for v in comp.values() if v*10 > idx) * 1.0/ dbnum

logfile.write("done.\n\n")
logfile.write("writing results to the output files...\n")

print cont_counts
print comp_counts

for gene in cont_list.keys():
	contfile.write(gene + "\t" + str(cont[gene]) + "\t" + cont_list[gene] + "\n")

for gene in comp_list.keys():
	compfile.write(gene + "\t" + str(comp[gene]) + "\t")
	for i in range(len(comp_list[gene]) - 1):
		compfile.write(comp_list[gene][i] + ",")
	compfile.write(comp_list[gene][-1] + "\n")

sumfile.write(" ".join(map(str, cont_counts)))
sumfile.write("\n")
sumfile.write(" ".join(map(str, comp_counts)))
sumfile.write("\n")

timecost = time.time() - start
logfile.write("done.\n")
logfile.write("all done. Took " + str(timecost) + " seconds\n")
logfile.write("Goodbye!\n")




logfile.close()
contfile.close()
compfile.close()
sumfile.close()
