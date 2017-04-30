#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import sys
from math import sqrt
from structureTools import distancePoints, centerMassOfResidue, computeDist_dico



def initBfactor(dPDB) :
	for chain in dPDB["nchaine"]:
		for res in dPDB[chain]["position"]:
			dPDB[chain][res]["bfactor"] = 0
    

def computeInterface(dPDB, threshold, mode) :
    
    initBfactor(dPDB)
    
    
    for chain1 in dPDB :
		if not chain1 == "nchaine":
			
			for chain2 in dPDB :
				if not chain2 == "nchaine" and not chain2 == chain1:
					
					for resi in dPDB[chain1]["position"] :
						
						for resj in dPDB[chain2]["position"] :
							
							#Calculate the distance between all residues
							distance = computeDist_dico(dPDB[chain1][resi], dPDB[chain2][resj], mode = mode)
							
							if not distance > threshold :# means, the two residues belong to the interface --> bfactor = 1
								if dPDB[chain1][resi]["bfactor"] == 0 :
									dPDB[chain1][resi]["bfactor"] = 1
								if dPDB[chain2][resj]["bfactor"] == 0 :
									dPDB[chain2][resj]["bfactor"] = 1
																					
							#Else the two residues doesn't belong to the interface --> bfactor = 0


def compareInterface (dPDB1, dPDB2) :
	
	nbInter = 0
	
	for chain in dPDB1 :
		if not chain == "nchaine":
			
			for res in dPDB1[chain]["position"] :
				
				#If both residues belong to the interface, the number of interface shared residues increases
				if dPDB1[chain][res]["bfactor"] and dPDB2[chain][res]["bfactor"] :
					nbInter += 1
					print dPDB1[chain][res]["bfactor"]
					
	return nbInter




