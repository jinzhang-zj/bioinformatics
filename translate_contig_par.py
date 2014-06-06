#!/opt/apps/python/epd/7.3.2/bin/python

"""Translate the transcriptome assemblies use Biopython (multithreading).

This a driver that translate the DNA sequences, and fix the common
assembly problem such as fragment, frame shift, the output includes
the normal translated contigs, the translated contigs with fixed 
frame shift or fragmentation, and the putative CDS.

The input is a file containing all the assembled coding sequences, 
and a referece (protein) database you want to use for blast.

Version 1.0

The output file will be in fasta format, the coordinates of hsp regions
and length of the extending regions is shown at the end of  the sequence 
description line.

This final output sequences are categorized in different region.
Those hsp (high-scoring segment pairs) regions are in upper case,
and those linker seqeunces between hsps and extending regions on both
5' and 3' ends are in lower case.

The Ns are added to correct the frame shift between different hsp
regions. Ns are added to the end of linker regions (sequences between hsp
regions). Since the linker region can be due to frame shift/insertion of
introns, it's not sure where to put Ns. Empirically the translated linker
sequences contain less stop codons when Ns are added to its end than to its
head.

The summary file reports the blast hit chosen, reading frame of first hsp,
used hsps, 5' and 3' extension on the final cds.

Differences with script from Dr. Mikhail Matz. (If anyone's interested)
---------------------------------------------------------------------------------
Problems:				my_solution / Matz's solution

Reference hit selection: pick up the reference containing highest number
of conserverd residues aligned to the query / pick up the one that maximize the 
fraction of conserverd residues times non overlapping hsp length

Linker regions:	add Ns to end of linker and keep linker sequences in cds / replace
the whole linker region with Ns and add more Ns to correct frame shift

Linker translation:	translate linker regions and show them in lower case / doesn't
translate linker regions

Multithreading / Singlethreading
---------------------------------------------------------------------------------
"""

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 1.0 $"
__date__ = "$Date: 2014/06/03 14:02:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

import sys
import re
import os
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastxCommandline
from multiprocessing import Pool
import time
#import math

# although xml is not the easiest human readable file format of blast results
# its overall structure remain stable cross different blast versions, and we
# will only take xml format in our script.

if len(sys.argv) < 3:
	sys.exit('Usage: ' + sys.argv[0] + ' [fastafile] [database] [num cpus]\n' +
	'''
The output cotains four files, assume the input fasta file is A.fna
The output files will be:
	A_cds.fna		coding sequences identified by blast
	A_translated.faa	translated protein sequences from cds file using first reading frame
	A_summary.out		summary of blast results, selected hits, frame, and hsp regions
	A_complete.faa		translated sequence of the input in all six reading frame

This script requires Biopython installation.

Jin Zhang
University of Texas at Austin
zj@utexas.edu
Nov 2014
''')

# iterator generation
def inputItr(blast_records, qseqs):
	path_pattern = ' path=\[.*\]'
	for blast_record in blast_records:
		cur_query = re.sub(path_pattern, '', blast_record.query)	#query name
		qseq = qseqs[cur_query]
		yield (blast_record, qseq, cur_query)

# function to process blast results
def processBlast ((blast_record, qseq, cur_query)):
	skip_pattern = re.compile('hypothetical|predicted|unnamed|unknown', re.I)
	ecut = 1e-5;
	hsp_region = []	# used hsps
	max_iden = 0	# maximal identities for the alignment between query and certain hit
	seq_length = len(qseq)
	cds_seq = ''
	pep_seq = ''
	comOut = ''
	cdsOut = ''
	pepOut = ''
	sumOut = ''
	# skip this query if there's not hit
	if not blast_record.alignments:
		return

	# each alignment corresponding to blast results of one hit for one query
	for alignment in blast_record.alignments:
		cur_hit = alignment.title	# hit name
		cur_iden = 0	# total identical residues for current hit
		# skip those reference sequences named hypothetical/predicted
		if skip_pattern.search(cur_hit):
			continue		
		
		# each hsp corresponding to blast results of one hsp of one hit for one query
		for hsp in alignment.hsps:
			if hsp.expect > ecut:
				continue
			# the intuition is to find the hit with maximum number of indentical residues
			cur_iden += hsp.identities

		# find hits with most identical alignment residues
		if cur_iden > max_iden:	
			max_iden = cur_iden
			raw_hit = alignment

	# If there's hit for this contig but none of them pass the alignment score threshold, we shall also skip this contig
	try:
		raw_hit
	except NameError:
		return

	raw_region = sorted([str(hsp.query_start) + '-' + str(hsp.query_end) for hsp in raw_hit.hsps])	# all hsps


	# filter those low quality hsps
	best_hit = [hsp for hsp in raw_hit.hsps if hsp.expect < ecut]	

	if len(best_hit) == 0:
		return

	# translate the sequences in all six reading frame
	comOut += '>' + cur_query + ' frame 1\n'
	comOut += qseq.translate().tostring() + '\n'
	comOut += '>' + cur_query + ' frame 2\n'
	comOut += qseq[1:].translate().tostring() + '\n'
	comOut += '>' + cur_query + ' frame 3\n'
	comOut += qseq[2:].translate().tostring() + '\n'
	comOut += '>' + cur_query + ' frame -1\n'
	comOut += qseq.reverse_complement().translate().tostring() + '\n'
	comOut += '>' + cur_query + ' frame -2\n'
	comOut += qseq.reverse_complement()[1:].translate().tostring() + '\n'
	comOut += '>' + cur_query + ' frame -3\n'
	comOut += qseq.reverse_complement()[2:].translate().tostring() + '\n'

	# find the subject start in the best alignment
	sub_start = 9999
	for hsp in best_hit:
		if hsp.sbjct_start < sub_start:
			sub_start = hsp.sbjct_start

	is_join_hsp = False
	is_frame_shift = False
	is_fragment = False
	is_inverse = False

	# translate based on the found best hit
	if len(best_hit) > 1:
		hsp_list = sorted(best_hit, key=lambda hsp: hsp.query_start)	# sort by query start
		hsp_region = [ str(hsp.query_start) + "-" + str(hsp.query_end) for hsp in hsp_list]	# used hsps
		
		pre = 0;		
		# check all the hsp in this hit is in the same direction
		for hsp in hsp_list:
			if hsp.frame[0] * pre < 0:	# two hsps are in different direction
				print 'Query: ' + cur_query + 'has inverse translated hsps.'
				is_inverse = True		
				break
			else:
				pre = hsp.frame[0]
		
		if is_inverse:
			return

		frame = hsp_list[0].frame[0]

		if frame > 0:
			qs = hsp_list[0].query_start
			qe = hsp_list[0].query_end

			cds_seq = qseq[qs - 1 : qe].tostring()	# doesn't include qe
			pep_seq = Seq(cds_seq).translate().tostring()
			# merge hsps
			for hsp in hsp_list[1:]:
				if hsp.query_end < qe:	
					continue	# ingore this shorter hsp

				shift = (hsp.query_start - qs) % 3
			
				if shift:
					is_frame_shift = True
	
				if hsp.query_start < qe:	# overlap
					is_join_hsp = True
					if not shift: # in frame
						add_region = qseq[qe : hsp.query_end]
						cds_seq += add_region.tostring()
						pep_seq += add_region.translate().tostring()
					else:			# out of frame, insert Ns to make it in the same frame
						insert = 'N' * (3 - shift)
						add_region = Seq(insert) + qseq[qe : hsp.query_end]
						cds_seq += add_region.tostring()
						pep_seq += add_region.translate().tostring()
				else:	# disjoint, regions between hsp will be in lower case
					linker = qseq[qe : hsp.query_start - 1].lower()	# linker region between hsps
					if not shift: # in frame
						add_region = qseq[hsp.query_start - 1 : hsp.query_end]
						cds_seq += linker.tostring() + add_region.tostring()
						pep_seq += linker.translate().tostring().lower() + add_region.translate().tostring()
					else:
						insert = 'N' * (3 - shift)
						# whether to insert the Ns before or after linker region? It's open to question
						add_region = qseq[hsp.query_start - 1: hsp.query_end]
						cds_seq += linker.tostring() + insert + add_region.tostring()
						pep_seq += (linker + Seq(insert)).translate().tostring().lower() + add_region.translate().tostring()
				# update the coding region coverage
				qe = hsp.query_end
		else:	
			hsp_list = sorted(hsp_list, key=lambda hsp: hsp.query_end)
			hsp_list.reverse()
			qs = hsp_list[0].query_start
			qe = hsp_list[0].query_end

			cds_seq = qseq[qs - 1 : qe].reverse_complement().tostring()	# get the reverse complementary sequences
			pep_seq = Seq(cds_seq).translate().tostring()
			
			# merge hsps
			for hsp in hsp_list[1:]:
				if hsp.query_start > qs:
					continue	# ignore shorted hsp			
			
				shift = (qe - hsp.query_end) % 3

				if shift:
					is_frame_shift = True
				
				if hsp.query_end > qs:	# overlap
					is_join_hsp = True
					if not shift:
						add_region = qseq[hsp.query_start - 1 : qs - 1].reverse_complement()
						cds_seq += add_region.tostring()
						pep_seq += add_region.translate().tostring()
					else:
						insert = 'N' * (3 - shift)
						add_region = Seq(insert) + qseq[hsp.query_start - 1: qs - 1].reverse_complement()
						cds_seq += add_region.tostring()
						pep_seq += add_region.translate().tostring()
				else:	# disjoint hsps
					linker = qseq[hsp.query_end : qs - 1].reverse_complement().lower()
					if not shift:	# in frame
						add_region = qseq[hsp.query_start - 1 : hsp.query_end].reverse_complement()
						cds_seq += linker.tostring() + add_region.tostring()
						pep_seq += linker.translate().tostring().lower() + add_region.translate().tostring()
					else:
						insert = 'N' * (3 - shift)
						add_region = qseq[hsp.query_start - 1 : hsp.query_end].reverse_complement()
						cds_seq += linker.tostring() + insert + add_region.tostring()
						pep_seq += (linker.tostring() + Seq(insert)).translate().tostring().lower() + add_region.translate().tostring()
				# update the coding region coverage
				qs = hsp.query_start

	else:	# only single hsp	
		hsp = best_hit[0]
		hsp_region = [str(hsp.query_start) + "-" + str(hsp.query_end)]
		frame = hsp.frame[0]
		qs = hsp.query_start
		qe = hsp.query_end
		cds_seq = (frame < 0) and qseq[qs - 1 : qe].reverse_complement().tostring() or qseq[qs - 1 : qe].tostring()
		pep_seq = Seq(cds_seq).translate().tostring()

	# extend the translated region
	fextend = Seq('')	# intialization
	bextend = Seq('')	# initialization

	is_bad_start = True
	if frame > 0:	# forward reading frame
		# check if current cds start with ATG/CTG/GTG, if so, it's not a bad start
		start_codon = qseq[qs - 1 : qs + 2].tostring().upper()
		if start_codon == 'ATG' or start_codon == 'CTG' or start_codon == 'GTG':
			is_bad_start = False

		# extend to 5' till meet the furthest ATG before meeting any stop codons
		startm = -1	# start region of codon 'ATG'
		start = qs - 4

		# if it's bad start (not ATG) or subject start is not close the head, extend the 5'
		if is_bad_start or sub_start > 5:	# extends
			while start >= 0:
				cur_codon = qseq[start : start + 3].tostring().upper()
				if cur_codon == 'TAA' or cur_codon == 'TAG' or cur_codon == 'TGA':
					break

				if cur_codon == 'ATG':
					startm = start
				start -= 3	# move forward	
		start += 3	# step back 3 bases
	
		# extend to 3' till meet stop codons such as TAA, TAG, TGA
		end = qe
		while end + 2 < seq_length:
			cur_codon = qseq[end : end + 3].tostring().upper()
			if cur_codon == 'TAA' or cur_codon == 'TAG' or cur_codon == 'TGA':
				end += 3
				break
			end += 3

		if start != qs -1:
			fextend = (startm == -1) and qseq[start : qs - 1].lower()  or qseq[startm : qs -1].lower()
			cds_seq = fextend.tostring() + cds_seq			
			pep_seq = fextend.translate().tostring().lower() + pep_seq

		if end != qe:
			bextend = qseq[qe : end].lower()
			cds_seq += bextend.tostring()
			pep_seq += bextend.translate().tostring().lower()

	else:
		# check if current cds starts with ATG/CTG/GTG, if so, it's not a bad start
		start_codon = qseq[qe - 3 : qe].reverse_complement().tostring().upper()
		if start_codon == 'ATG' or start_codon == 'GTG' or start_codon == 'CTG':
			is_bad_start = False

		# extend to 5' till meet the furthest ATG before meeting any stop codons
		startm = -1
		start = qe
		
		# if it's bad start (not ATG) or subject start is not close the head, extend the 5'
		if is_bad_start or sub_start > 5:
			while start + 2 < seq_length:
				cur_codon = qseq[start : start + 3].reverse_complement().tostring().upper()
				if cur_codon == 'TAA' or cur_codon == 'TAG' or cur_codon == 'TGA':
					break

				if cur_codon == 'ATG':
					startm = start + 3
				start += 3

		# extend to 3' till meet the stop codons such as TAA, TAG, TGA
		end = qs - 1
		while end -3 >= 0:
			end -= 3
			cur_codon = qseq[end : end + 3].reverse_complement().tostring().upper()
			if cur_codon == 'TAA' or cur_codon == 'TAG' or cur_codon == 'TGA':
				break
			
		if start != qe:
			fextend = (startm == -1) and qseq[qe : start].reverse_complement().lower() or qseq[qe : startm].reverse_complement().lower()
			cds_seq = fextend.tostring() + cds_seq
			pep_seq = fextend.translate().tostring().lower() + pep_seq
			
		if end != qs - 1:
			bextend = qseq[end : qs - 1].reverse_complement().lower()
			cds_seq += bextend.tostring()
			pep_seq += bextend.translate().tostring().lower()

	cdsOut += '>' + cur_query + " hsps: "
	cdsOut += ' '.join(hsp for hsp in hsp_region)
	cdsOut += ', extends: 5\' '+ str(len(fextend))+ 'bp' 
	cdsOut += ' 3\' ' + str(len(bextend)) + 'bp\n'
	cdsOut += cds_seq + '\n'

	pepOut += '>' + cur_query + ' translated'
	pepOut += ' extends: N '+ str(len(fextend)/3)+ 'aa' 
	pepOut += ' C ' + str(len(bextend)/3) + 'aa\n'
	pepOut += pep_seq+"\n"

	sumOut += '----------------------------------\n';
	sumOut += 'query: ' + cur_query + '\n'
	sumOut += 'best hit: ' + raw_hit.title + '\n'
	sumOut += 'frame: ' + str(frame) + '\n'
	sumOut += 'total hsps: '
	sumOut += ' '.join(hsp for hsp in raw_region)
	sumOut += '\n'

	if is_join_hsp:
		sumOut += 'joined hsp: y\n'
	else:
		sumOut += 'joined hsp: n\n'

	if is_frame_shift:
		sumOut += 'frame shift: y\n'
	else:
		sumOut += 'frame shift: n\n'

	sumOut += 'used hsps: '
	sumOut += ' '.join(hsp for hsp in hsp_region)
	sumOut += '\n'
	sumOut += '5\' extension: ' + str(len(fextend)/3) + ' codons\n'
	sumOut += '3\' extension: ' + str(len(bextend)/3) + ' codons\n'
	sumOut += 'in frame stop codons: ' + str(pep_seq.count('*')) + '\n'


	#rout = reduce(lambda a, b: math.log(a+b), xrange(10**5))
	
	return (comOut, cdsOut, pepOut, sumOut)

if __name__ == "__main__":
	start = time.time()
	qseqs = {}

	# read in the fasta file
	for contig in SeqIO.parse(sys.argv[1], "fasta"):
		path_pattern = ' path=\[.*\]'
		tmp_description = re.sub(path_pattern, '', contig.description)
		qseqs[tmp_description]= contig.seq

	file = os.path.basename(sys.argv[1]).replace('.fna','')				# input query file in fasta format
	file = file.replace('.fasta','')				# input query file in fasta format

	cds_file = open(str(file) + '_cds.fna', 'w')			# final output of cds file
	pep_file = open(str(file) + '_translated.faa', 'w')	# final output of the translated sequence file
	sum_file = open(str(file) + '_summary.out', 'w')		# summary of the translation results
	com_file = open(str(file) + '_complete.faa', 'w')	# translation of all six reading frame

	# running blastx to identify the reading frame
	blastxcline = NcbiblastxCommandline(query=sys.argv[1], db=sys.argv[2], max_target_seqs=10, num_threads=sys.argv[3],evalue=1e-5, outfmt=5, out=file+".blastx.out")
	stdout, stderr = blastxcline()

	blasttime = time.time() - start

	# open/parse blast results
	result_handle = open(file+".blastx.out")
	blast_records = NCBIXML.parse(result_handle)

	# process each blast results
	results = []
	pool = Pool(processes=int(sys.argv[3]))
	print "Using " + sys.argv[3] + " processes"
	myItr = inputItr(blast_records, qseqs)
	r = pool.map_async(processBlast, myItr, callback=results.extend)
	r.wait()

	processtime = time.time() - start - blasttime

	# write the results into corresponding files
	for result in results:
		if result:
			com_file.write(result[0])
			cds_file.write(result[1])
			pep_file.write(result[2])
			sum_file.write(result[3])

	timecost = time.time() - start
	writetime = time.time() - start - blasttime - processtime

	print "time performance:\n"
	print "read data & blas: " + str(blasttime) + " seconds\n";
	print "process blast results: " + str(processtime) + " seconds\n";
	print "write results: " + str(writetime) + " seconds\n";
	print "total cost: " + str(timecost) + " seconds\n"

	cds_file.close()
	pep_file.close()
	sum_file.close()
	com_file.close()
