#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import math
from math import sqrt
from ParserPDB import parserPDB
import sys
from Usage import usage
from structureTools import distancePoints

	
	
def computeRMSD (dPDB1, dPDB2, l_atoms):
	"""Compute RMSD between the two PDB dictionnaries on the given atom list
	input : 2 PDB dictionaries
	output : RMSD score"""

	N = 0
	
	distanceSum = 0
	
	for chain in dPDB1 :
		
		if chain != "nchaine" :
			
			for pos in dPDB1[chain] :
				
				if pos != "position" :
			
					for atom in dPDB1[chain][pos] :
				
						if atom in l_atom :
							
							distanceSum += pow(distancePoints((dPDB1[chain][pos][atom]["x"],dPDB1[chain][pos][atom]["y"], dPDB1[chain][pos][atom]["z"]),(dPDB2[chain][pos][atom]["x"], dPDB2[chain][pos][atom]["y"], dPDB2[chain][pos][atom]["z"])), 2)
							N += 1		
							
				
	RMSD = sqrt(distanceSum/N)		
	return RMSD


	
if __name__ == '__main__':

	try:
		infile1 = sys.argv[sys.argv.index("-pdb1")+1]
		infile2 = sys.argv[sys.argv.index("-pdb2")+1]
		atomParam = sys.argv[sys.argv.index("-atom")+1]
		print "pdb files to treat:", infile1, infile2, atomParam
		d1 = parserPDB(fichier1)
		d2 = parserPDB(fichier2)
		l_atomParam = []
		if atomParam == NULL :
			l_atomParam = ["CA", "N", "C", "O"]
			print computeRMSD(d1, d2, l_atomParam)
		elif atomParam == ALL :
			l_atomParam = ["atom"]
			print computeRMSD(d1, d2, l_atomParam)
		else :
			for x in atomParam.split("_"):
				l_atomParam.append(x); 
			print computeRMSD(d1, d2, l_atomParam)
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files\n"
		sys.exit()


