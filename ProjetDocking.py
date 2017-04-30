from structureTools import scorelist, writePDB
from parsing import parserPDB, preparePDB, parseDirectory
from Usage import usage
from RMSD import computeRMSD
import sys
from path import path
import re
from Interface import computeInterface, compareInterface



if __name__ == '__main__':

	'''
	recuperation des arguments du fichier recepteur et du dossier solutions ligand
	'''
	try:
		infile = sys.argv[sys.argv.index("-pdb1")+1]
		indir = sys.argv[sys.argv.index("-pdb2")+1]
		
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files to analyze\n"
		sys.exit()

	'''
	Parsing et preparation du fichier du recepteur
	'''
	try:
		dico_Rec = preparePDB(infile)
	except:
		print "ERROR : Parsing of "+infile+" failed"
	
	'''
	Parsing du repertoire contenant les configurations de ligand
	Determination des scores
	'''
	try:
		listFiles, listScores = parseDirectory(indir, dico_Rec)
	except:
		print "ERROR : Parsing of "+indir+" failed"

	'''
	Determination du meilleur score
	'''
	try:
		filenum = int(scorelist(listFiles, listScores))
	except:
		print "ERROR : Could not find the file corresponding to the best score"
	
	filename="1BRS_A_1BRS_B_allatom_"+str(filenum)+"_DP.pdb"
	pathfile=[indir,'/',filename]
	dico_Lig = parserPDB(''.join(pathfile))
	
	'''
	Ecriture du fichier contenant le classement des scores
	'''
	try:
		writePDB(dico_Rec, dico_Lig, prediction = "complexe_predit_score1.pdb")
	except:
		print "ERROR : Could not create list of scoring"
		
	'''
	Recuperation des fichiers pdb du ligand, recepteur et complexe natif donnes en argument
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
	Parsing des fichiers pdb du ligand, recepteur et complexe natif
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
		dico_compl_predit = parserPDB("complexe_predit_score1.pdb")
		
	except :
		print "EROOR : File parsing of complexe_predit_score1.pdb failed"
	
	
	try :
		dico_compl_natif = parserPDB(nativeComplex)
		
	except : 
		print "ERROR : File parsing of "+nativeComplex+" failed"
			
	'''
	Parametrisation de la fonction computeRMSD
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
		
	
	try :
		'''
		Calcul du RMSD
		dico_Lig et dico_lig_natif sont les dictionnaires correspondants aux pdb de la meilleure solution ligand et du ligand natif (-pdb4)
		dico_Rec et dico_rec_natif sont les dictionnaires correspondants aux pdb du recepteur et du recepteur natif (-pdb5)
		dico_compl_predit et dico_compl_natif sont les dictionnaires correspondants aux pdb du complexe predit et du complexe natif (-pdb3)
		'''
		print "ligand RMSD : "+str(computeRMSD(dico_Lig, dico_lig_natif, l_atomParam))
		print "receptor RMSD : "+str(computeRMSD(dico_Rec, dico_rec_natif, l_atomParam))
		print "complex RMSD : "+str(computeRMSD(dico_compl_predit, dico_compl_natif, l_atomParam))

	except :
		print "ERROR : RMSD calculation failed."


	'''
	Parametrisation de la fonction computeInterface
	'''	
	try :
		thresholdParam = float(sys.argv[sys.argv.index("-threshold")+1])
	except :
		thresholdParam = 2
		
		
	try :
		interfaceModeParam = (sys.argv[sys.argv.index("-mode")+1])
	except :
		interfaceModeParam = "center"
	
	
	try :
		'''
		Calcul des interfaces du cmplexe natif et du complexe predit
		'''
		computeInterface(dico_compl_natif, thresholdParam, interfaceModeParam)
		computeInterface(dico_compl_predit, thresholdParam, interfaceModeParam)
		
	except :
		print "ERROR : Interface calculation failed."
		
	
	try :
		'''
		Comparaison des interfaces
		'''
		nbresInter = compareInterface(dico_compl_natif, dico_compl_predit)
		
	except :
		print "ERROR : Interface comparison failed."
		
		
		
	print "Number of native contacts: "+str(nbresInter)
		
		
		
