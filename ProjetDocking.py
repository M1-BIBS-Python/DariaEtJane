from structureTools import scorelist, writePDB
from parsing import parserPDB, preparePDB, parseDirectory
from Usage import usage
from RMSD import computeRMSD
import sys
from path import path
import re



if __name__ == '__main__':

	'''
	recuperation des arguments du fichier recepteur et du dossier solutions ligand
	'''
	try:
		infile = sys.argv[1]
		indir = sys.argv[2]
		
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
		nativeComplex = sys.argv[3]
		nativeLigand = sys.argv[4]
		nativeReceptor = sys.argv[5]
		
	except :
		usage()
		print "ERROR: please enter the names of the pdb files to compare\n"
		sys.exit()
		
	'''
	Parsing des fichiers pdb du ligand, recepteur et complexe natif
	'''
	try :
		d2 = parserPDB(nativeLigand)
		
	except : 
		print "ERROR : File parsing of "+nativeLigand+" failed"
		
		
	try :
		d4 = parserPDB(nativeReceptor)
		
	except : 
		print "ERROR : File parsing of "+nativeReceptor+" failed"
	
	try :
		d5 = parserPDB("complexe_predit_score1.pdb")
	except :
		print "EROOR : File parsing of complexe_predit_score1.pdb failed"
	
	try :
		d6 = parserPDB(nativeComplex)
		
	except : 
		print "ERROR : File parsing of "+nativeComplex+" failed"
			
	'''
	Parametrisation de la fonction computeRMSD
	'''	
	try :
		atomParam = sys.argv[6]
		l_atomParam = []
		if atomParam == ALL :
			l_atomParam = ["atom"]
		else :
			for x in atomParam.split("_"):
				l_atomParam.append(x); 
	
	except : 
		l_atomParam = ["CA", "N", "C", "O"]
		print computeRMSD(dico_Lig, d2, l_atomParam)
		
	
	try :
		'''
		Calcul du RMSD
		dico_Lig et d2 sont les dictionnaires correspondants aux pdb de la meilleure solution ligand et du ligand natif (-pdb4)
		dico_Rec et d4 sont les dictionnaires correspondants aux pdb du recepteur et du recepteur natif (-pdb5)
		d5 et d6 sont les dictionnaires correspondants aux pdb du complexe predit et du complexe natif (-pdb3)
		'''
		print "ligand RMSD : "+computeRMSD(dico_Lig, d2, l_atomParam)
		print "receptor RMSD : "+computeRMSD(dico_Rec, d4, l_atomParam)
		print "complex RMSD : "+computeRMSD(d5, d6, l_atomParam)

	except :
		print "ERROR : RMSD calculation failed."
