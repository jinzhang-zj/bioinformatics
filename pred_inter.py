'''Calculating distance of all atoms between given protein chains.

Version 0.1

The script require multiprocessing.
'''

__author__ = "Jin Zhang(zj@utexas.edu)"
__version__ = "$Version: 0.1 $"
__date__ = "$Date: 2014/04/05 11:27:15"
__copyright__ = "Copyright: (c) 2014 Jin Zhang"
__license__ = "Python"

import sys
import argparse
import Bio.PDB
import numpy as np
from multiprocessing import Pool

def calc_res_dist((cutoff, A, B,res_one, res_two)) :
    """Returns the C-alpha distance between two residues"""
    result = []
    for atomA in res_one:
        for atomB in res_two:
            diff_vec = atomA.coord - atomB.coord
            dis = np.sqrt(np.sum(diff_vec * diff_vec))
	    if dis <= cutoff:
            	output = (res_one.id[1], A, res_one.resname, 
                	                        atomA.name, res_two.id[1], B, res_two.resname, atomB.name, dis)
		result.append(output)
    return result

# generate the arg parser
parser = argparse.ArgumentParser(
	#usage='%(prog)s [options]',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	description='get all interaction residues under certain cutoff')
	
# add arguments
parser.add_argument('-f', metavar='pdbfile', nargs=1, required=True,help='pdf structure file')
parser.add_argument('-c1', metavar='chain_one', nargs=1, required=True,help='''
		chain list one for pairwise comparison; 
		if c1 = abc, c2=de, the chain pairs to be compared is:
		ad, ae, bd, be, cd, ce	
		''')
parser.add_argument('-c2', metavar='chain_two', nargs=1, required=True,help='''
		chain list one for pairwise comparison; 
		if c1 = abc, c2=de, the chain pairs to be compared is:
		ad, ae, bd, be, cd, ce	
		''')
parser.add_argument('-o', metavar='output_file', nargs=1, required=True,help='name of output files')
parser.add_argument('-cpu', metavar='num_of_threads', default=1, help='number of threads to use')
parser.add_argument('-cut', metavar='distance_cutoff', default=10, help='distance cutoff between two residues')
parser.add_argument('-v','--version',action='version', version='%(prog)s 2.0')

args = parser.parse_args()	# take sys.args as default

# dictionary that stores all the arguments
arglist = vars(args)

if __name__ == "__main__":
    	pdb_filename = arglist['f'][0]
    	pdb_code = pdb_filename[:4]
	alist = arglist['c1'][0].upper()
    	blist = arglist['c2'][0].upper()
    	num = int(arglist['cpu'])
    	ofile = open( arglist['o'][0], 'w' )
    	cutoff = int(arglist['cut'])
   
	
    	structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
   	model = structure[0]
    	pool = Pool(processes=num)
    	res_pairs = []

    	for A in alist:
		for B in blist:
			chain_one = model[A]
			chain_two = model[B]

			for res_one in chain_one:
				for res_two in chain_two:
					if res_one.resname == ' MG' or res_one.resname == 'HOH' or res_two.resname == ' MG' or res_two.resname == 'HOH':
						continue
					else:
						res_pairs.append((cutoff, A,B,res_one,res_two))

	# calculating the residue distance using multiprocessing		
    	results = []
    	r = pool.map_async(calc_res_dist, res_pairs, callback=results.append)
    	r.wait()

    	for entry in results[0]:
		for atomlevel in entry:
			ofile.write("{:<4d} {:<3s} {:<6s} {:<6s} {:<4d} {:<3s} {:<6s} {:<6s} {:<6.2f}\n".format(atomlevel[0], atomlevel[1], 
							atomlevel[2], atomlevel[3], atomlevel[4], atomlevel[5], atomlevel[6], atomlevel[7], float(atomlevel[8])))
	ofile.close()
