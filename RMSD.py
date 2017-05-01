#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import math
from math import sqrt
import sys
from Usage import usage
from StructureTools import distancePoints
from parsing import parserPDB

	
def computeRMSD (dPDB1, dPDB2, l_atoms, bfactor):
	"""Compute RMSD between the two PDB dictionnaries on the given atom list
	input : 2 PDB dictionaries
	output : RMSD score"""

	N = 0
	
	distanceSum = 0
	
	RMSD = 100000
	
	for chain in dPDB1 :
		
		if chain != "nchaine" :
			
			for pos in dPDB1[chain] :
				
				if pos != "position" :
			
					for atom in dPDB1[chain][pos] :
				
						if atom in l_atoms :
							
							if bfactor :
								
								if dPDB1[chain][pos]["bfactor"] :
									
									distanceSum += pow(distancePoints((dPDB1[chain][pos][atom]["x"],dPDB1[chain][pos][atom]["y"], dPDB1[chain][pos][atom]["z"]),(dPDB2[chain][pos][atom]["x"], dPDB2[chain][pos][atom]["y"], dPDB2[chain][pos][atom]["z"])), 2)
									N += 1		
							else :
								
								distanceSum += pow(distancePoints((dPDB1[chain][pos][atom]["x"],dPDB1[chain][pos][atom]["y"], dPDB1[chain][pos][atom]["z"]),(dPDB2[chain][pos][atom]["x"], dPDB2[chain][pos][atom]["y"], dPDB2[chain][pos][atom]["z"])), 2)
								N += 1	
				
	if N!=0:
		RMSD = sqrt(distanceSum/N)		
	
	return RMSD


'''
if __name__ == '__main__':

	try:
		infile1 = sys.argv[sys.argv.index("-pdb1")+1]
		infile2 = sys.argv[sys.argv.index("-pdb2")+1]
		print "pdb files to treat:", infile1, infile2
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files\n"
		sys.exit()
		
	try :
		d1 = parserPDB(infile1)
		
	except : 
		print "ERROR : Echec du parsing du fichier "+infile1
		
	try :
		d2 = parserPDB(infile2)
	except : 
		print "ERROR : Echec du parsing du fichier "+infile2
	
	
	try :
		atomParam = sys.argv[sys.argv.index("-atom")+1]
		l_atomParam = []
		if atomParam == ALL :
			l_atomParam = ["atom"]
			print computeRMSD(d1, d2, l_atomParam)
		else :
			for x in atomParam.split("_"):
				l_atomParam.append(x); 
			print computeRMSD(d1, d2, l_atomParam)
	except : 
		l_atomParam = ["CA", "N", "C", "O"]
		print computeRMSD(d1, d2, l_atomParam)
		
		
		
'''
