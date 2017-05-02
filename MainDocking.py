from StructureTools import scorelist, writePDB, extract100Bests, writeInterfacePDB
from parsing import parserPDB, preparePDB, parseDirectory
from Usage import usage
from RMSD import computeRMSD
import sys
from path import path
import re
from Interface import computeInterface, compareInterface, computeInterfaceScore, interfaceHydrophobicity



if __name__ == '__main__':
	
	#################################
	#	ARGUMENTS AND PARAMETERS	#
	#		   RETRIEVAL			#
	#################################
	
	'''
	Receptor pdb file and ligand solutions directory retrieval
	'''
	try:
		infile = sys.argv[sys.argv.index("-pdb1")+1]
		indir = sys.argv[sys.argv.index("-pdb2")+1]
		
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files to analyze\n"
		sys.exit()
	
	'''
	Native receptor, native ligand and native complex pdb files retrieval
	'''
	
	try:
		nativeComplex = sys.argv[sys.argv.index("-pdb3")+1]
		nativeLigand = sys.argv[sys.argv.index("-pdb4")+1]
		nativeReceptor = sys.argv[sys.argv.index("-pdb5")+1]
		
	except :
		usage()
		print "ERROR: please enter the names of the pdb files to compare\n"
		sys.exit()
	
	'''
	Configuration of the computeInterface function
	'''	
	try :
		thresholdParam = float(sys.argv[sys.argv.index("-threshold")+1])
	except :
		thresholdParam = 2
		
		
	try :
		interfaceModeParam = (sys.argv[sys.argv.index("-mode")+1])
	except :
		interfaceModeParam = "atom"
	
	'''
	Configuration of the computeRMSD function
	'''	
	try :
		atomParam = sys.argv[sys.argv.index("-atom")+1]
		l_atomParam = []
		if atomParam == ALL :
			l_atomParam = ["atom"]
		else :
			for x in atomParam.split("_"):
				l_atomParam.append(x); 
	
	except : 
		l_atomParam = ["CA", "N", "C", "O"]

	#########################
	#	INITIALISATION		#
	#	AND PARSING			#
	#########################
	
	'''
	Native receptor, natve ligand and natif complex pdb files parsing
	'''
	try :
		dico_lig_natif = parserPDB(nativeLigand)
		
	except : 
		print "ERROR : File parsing of "+nativeLigand+" failed"
		
		
	try :
		dico_rec_natif = parserPDB(nativeReceptor)
		
	except : 
		print "ERROR : File parsing of "+nativeReceptor+" failed"
	
	try :
		dico_compl_natif = parserPDB(nativeComplex)
		
	except : 
		print "ERROR : File parsing of "+nativeComplex+" failed"
			
	'''
	Preparation of the receptor pdb file
	'''
	try:
		dico_Rec = preparePDB(infile)
	except:
		print "ERROR : Parsing of "+infile+" failed"
	
	
	#########################
	#	   FIRST METHOD     #
	#########################
	
	#DETERMINATION OF THE RIGHT COMPLEX SOLUTION
	
	print "Start of the first method"
	
	'''
	Ligand directory parsing
	Score determination
	'''
	try:
		listFiles, listScores = parseDirectory(indir, dico_Rec)
	except:
		print "ERROR : Parsing of "+indir+" failed"	
	

	'''
	Best score determination
	'''
	try:
		filenum = int(scorelist(listFiles, listScores, "Scoring1.txt"))
	except:
		print "ERROR : Could not find the file corresponding to the best score"
	
	filename="1BRS_A_1BRS_B_allatom_"+str(filenum)+"_DP.pdb"
	pathfile=[indir,'/',filename]
	dico_Lig = parserPDB(''.join(pathfile)) #Dictionnary of the ligand solution
	
	'''
	Writing the file containing the score ranking according to the first method
	'''
	try:
		writePDB(dico_Rec, dico_Lig, prediction = "complexe_predit_score1.pdb")
	except:
		print "ERROR : Could not create list of scoring"
	
	'''
	First method predicted complex pdb file parsing
	'''
	try :
		dico_compl_predit = parserPDB("complexe_predit_score1.pdb")
		
	except :
		print "ERROR : File parsing of complexe_predit_score1.pdb failed"
	
	
	#SOLUTION EVALUATION
	
	'''
	Calculation of the native and predicted complex interfaces
	'''
	try :
		computeInterface(dico_compl_natif, thresholdParam, interfaceModeParam)
		computeInterface(dico_compl_predit, thresholdParam, interfaceModeParam)
		
	except :
		print "ERROR : Interface calculation failed."
		
	'''
	Interface comparison
	'''
	try :
		nbresInter = compareInterface(dico_compl_natif, dico_compl_predit)
		
	except :
		print "ERROR : Interface comparison failed."

	
	print "Number of native contacts: "+str(nbresInter)
		
	'''
	RMSD calculation
	dico_Lig and dico_lig_natif are dictionaries corresponding to the best ligand solution and native ligand (-pdb4) pdb files
	dico_compl_predit and dico_compl_natif are dictionaries corresponding to the predicted complexe and native complex (-pdb3) pdb files
	bfactor indicate if the calculation is made only on the interface
	'''
	try :
		print "First solution :"
		print "\t->ligand RMSD : "+str(computeRMSD(dico_Lig, dico_lig_natif, l_atomParam, bfactor = False))
		print "\t->complex RMSD : "+str(computeRMSD(dico_compl_predit, dico_compl_natif, l_atomParam, bfactor = False))
		print "\t->interface RMSD : "+str(computeRMSD(dico_compl_predit, dico_compl_natif, l_atomParam, bfactor = True))

	except :
		print "ERROR : RMSD calculation failed."
	
	'''
	Editing of the complexe_predit_score1.pdb to add the bfactor
	For easy manipulation with Pymol
	'''
	try:
		writeInterfacePDB(dico_compl_predit,"complexe_predit_score1.pdb")
	except:
		print "ERROR : Could not edit complexe_predit_score1.pdb"
		
	#########################
	#	  SECOND METHODE    #
	#########################
	
	#SOLUTION DETERMINATION
	
	print "Start of the second method"
	
	'''
	Extraction of the 100 best results by thefirst method
	'''
	try:
		listBestHits, listBestScores, listBestFiles = extract100Bests(indir)
	except:
		print "ERROR : Could not extract the 100 best hits"
	
	'''
	Calculation of the second score considering hydrophobicity of the interface residues
	'''
	try:
		listNewScore=interfaceHydrophobicity(listBestHits, listBestScores, dico_Rec, thresholdParam, interfaceModeParam)
	except:
		print "ERROR : Could not compute the new score taking the hydrophobicity in paramaters"

	filenum = int(scorelist(listBestFiles, listNewScore, "Scoring2.txt"))
	
	filename="1BRS_A_1BRS_B_allatom_"+str(filenum)+"_DP.pdb"
	pathfile=[indir,'/',filename]
	dico_Lig = parserPDB(''.join(pathfile)) #Dictionnary of the ligand solution
	
	'''
	Writing the file containing the score ranking according to the second method
	'''
	try:
		writePDB(dico_Rec, dico_Lig, prediction = "complexe_predit_score2.pdb")
	except:
		print "ERROR : Could not create list of scoring"
	
	'''
	Second method predicted complex pdb file parsing
	'''
	try :
		dico_compl_predit = parserPDB("complexe_predit_score2.pdb")
		
	except :
		print "ERROR : File parsing of complexe_predit_score1.pdb failed"

	#EVALUATION DE LA SOLUTION
	
	'''
	Calculation of the native and predicted complex interfaces
	'''
	try :
		computeInterface(dico_compl_natif, thresholdParam, interfaceModeParam)
		computeInterface(dico_compl_predit, thresholdParam, interfaceModeParam)
		
	except :
		print "ERROR : Interface calculation failed."
		
	'''
	Interface comparison
	'''
	try :
		nbresInter = compareInterface(dico_compl_natif, dico_compl_predit)
		
	except :
		print "ERROR : Interface comparison failed."

	print "Number of native contacts: "+str(nbresInter)
		
	'''
	RMSD calculation
	dico_Lig and dico_lig_natif are dictionaries corresponding to the best ligand solution and native ligand (-pdb4) pdb files
	dico_compl_predit and dico_compl_natif are dictionaries corresponding to the second method predicted complexe and native complex (-pdb3) pdb files
	bfactor indicate if the calculation is made only on the interface
	'''
	try :
		print "Second solution :"
		print "\t->ligand RMSD : "+str(computeRMSD(dico_Lig, dico_lig_natif, l_atomParam, bfactor = False))
		print "\t->complex RMSD : "+str(computeRMSD(dico_compl_predit, dico_compl_natif, l_atomParam, bfactor = False))
		print "\t->interface RMSD : "+str(computeRMSD(dico_compl_predit, dico_compl_natif, l_atomParam, bfactor = True))

	except :
		print "ERROR : RMSD calculation failed."
	
		'''
	Editing of the complexe_predit_score1.pdb to add the Bfactor
	For easy manipulation with Pymol
	'''
	try:
		writeInterfacePDB(dico_compl_predit,"complexe_predit_score2.pdb")
	except:
		print "ERROR : Could not edit complexe_predit_score1.pdb"
