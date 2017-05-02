#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import sys
from math import sqrt
from StructureTools import distancePoints, centerMassOfResidue, computeDist_dico
from ResiduesProperties import hydrophobicResidue
from parsing import parserPDB



def initBfactor(dPDB) :
	"""
	Initialize a "bfactor" key to "position" dictionary to label those which belong to the interface
    Input: dico dPDB 
    Output: dico dPDB with bfactor for each residue
    """
	for chain in dPDB["nchaine"]:
		for res in dPDB[chain]["position"]:
			dPDB[chain][res]["bfactor"] = 0


def computeInterfaceScore(dPDB1, dPDB2, threshold, mode) :
    """
    compute interface hydrophobic/hydrophilic score proportion by comparing hydrophobic properties of each atom
    which belong to the interface
    Input: dico dPDB1 and dPDB2 corresponding to the receptor and ligand, threshold distance for the interface,
    mode for the mode of calculation of the distance (center or atom)
    Output: the hydrophobic/hydrophilic proportion at the interface between ligand and receptor
    """
    #Initialize bfactor
	initBfactor(dPDB1)
	initBfactor(dPDB2)
    
	nbContactHH = 0
	nbResInterface = 0
	proportion = 0
    
	for chain1 in dPDB1 :
		if not chain1 == "nchaine":
			for res1 in dPDB1[chain1]["position"] :
				for chain2 in dPDB2 :
					if not chain2 == "nchaine" :
						for res2 in dPDB2[chain2]["position"] :
							#Calculate the distance between all residues
							distance = computeDist_dico(dPDB1[chain1][res1], dPDB2[chain2][res2], mode = mode)
							if not distance > threshold : # means, the two residues belong to the interface							
								nbResInterface += 1
								#Compare the interface residues hydrophobicity
								if ((hydrophobicResidue(dPDB1[chain1][res1]['residu'])) and not(hydrophobicResidue(dPDB2[chain2][res2]['residu']))) or (not(hydrophobicResidue(dPDB1[chain1][res1]['residu'])) and hydrophobicResidue(dPDB2[chain2][res2]['residu'])):
									nbContactHH += 1
	#If there is an interface								
	if not nbResInterface == 0:
		proportion = (nbContactHH)/(nbResInterface)		
	
	return proportion


def interfaceHydrophobicity(listBestHits, listBestScores, dico_Rec, thresholdParam, interfaceModeParam):

	listPercentage=[]
	for files in listBestHits:
		dico=parserPDB(files)
		percentage=int(computeInterfaceScore(dico_Rec, dico, thresholdParam, interfaceModeParam))
		listPercentage.append(percentage)
	
	listnewScore=[(1-a)*b for a,b in zip(listPercentage,listBestScores)]
	
	return listnewScore
    

def computeInterface(dPDB, threshold, mode) :
    """
    compute interface of a ligand-receptor complex
    Input: dico dPDB corresponding to the ligand-receptor complex, threshold distance for the interface, mode 
    for the mode of calculation of the distance (center or atom)
    Output: dico dPDB with modified bfactor at the interface and the number of residues which belong to the 
    interface
    """
    #Initialize bfactor
	initBfactor(dPDB)
	#The number of residues with belong to the interface
	nbResInter = 0
	
	for chain1 in dPDB :
		if not chain1 == "nchaine":
			
			for chain2 in dPDB :
				if not chain2 == "nchaine" and not chain2 == chain1:
					
					for resi in dPDB[chain1]["position"] :
						
						for resj in dPDB[chain2]["position"] :
							
							#Calculate the distance between all residues
							distance = computeDist_dico(dPDB[chain1][resi], dPDB[chain2][resj], mode = mode)
							
							if not distance > threshold :# means, the two residues belong to the interface --> bfactor = 1, nbResInter++
								if dPDB[chain1][resi]["bfactor"] == 0 :
									dPDB[chain1][resi]["bfactor"] = 1
								if dPDB[chain2][resj]["bfactor"] == 0 :
									dPDB[chain2][resj]["bfactor"] = 1
								nbResInter += 1
	return nbResInter													





def compareInterface (dPDB1, dPDB2) :
	"""
    compare the residues with belong to the interface between two ligand-receptor complexes
    Input: dico dPDB1 and dPDB2 corresponding to ligand-receptor complexes to compare
    Output: the number of complexes interface shared residues
    """
    
	nbInter = 0

	
	for chain in dPDB1 :
		if not chain == "nchaine":
			
			for res in dPDB1[chain]["position"] :
				
				#If both residues belong to the interface, the number of interface shared residues increases
				if dPDB1[chain][res]["bfactor"] and dPDB2[chain][res]["bfactor"] :
					nbInter += 1
								
	return nbInter


