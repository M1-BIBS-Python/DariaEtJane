from structureTools import parserPDB, preparePDB, scorelist
from RunScore import compEner
import sys
from path import path
import re



if __name__ == '__main__':

	'''
	recuperation des arguments du fichier recepteur et du dossier solutions ligand
	calcul des scores et choix du meilleur
	'''
	try:
		infile = sys.argv[sys.argv.index("-pdb1")+1]
		indir = sys.argv[sys.argv.index("-pdb2")+1]
		print "pdb files to treat:", infile1, infile2
		dico_Rec = preparePDB(infile)
		dico_Scores = {}
		n=0
		for files in path("indir").walkfiles():
			if file.endswith(".pdb") :
				dico_confH = preparePDB(files)
				dico_Scores[files]= compEner(dico_Rec, dico_confH, "B", "D")
				dico_confH = {}
		
		scorelist(dico_Scores)
		bestScore = float(max(dico_Scores.values()))

		print bestScore
		
	except:    
		usage()
		print "ERROR: please enter the names of the pdb files\n"
		sys.exit()
		
	'''
	Recuperation des fichiers pdb du ligand, recepteur et complexe natif donnes en argument
	'''
	try:
		nativeComplex = sys.argv[sys.argv.index("-pdb3")+1]
		nativeLigand = sys.argv[sys.argv.index("-pdb4")+1]
		nativeReceptor = sys.argv[sys.argv.index("-pdb5")+1]
		
	except :
		usage()
		print "ERROR: please enter the names of the pdb files\n"
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
		d6 = parserPDB(nativeComplex)
		
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
		print computeRMSD(d1, d2, l_atomParam)
		
	
	try :
		'''
		Calcul du RMSD
		d1 et d2 sont les dictionnaires correspondants aux pdb de la meilleure solution ligand et du ligand natif (-pdb4)
		dico_Rec et d4 sont les dictionnaires correspondants aux pdb du recepteur et du recepteur natif (-pdb5)
		d5 et d6 sont les dictionnaires correspondants aux pdb du complexe predit et du complexe natif (-pdb3)
		Libre a toi de les renommer bien sur
		'''
		print "ligand RMSD : "+computeRMSD(d1, d2, l_atomParam)
		print "receptor RMSD : "computeRMSD(dico_Rec, d4, l_atomParam)
		print "complex RMSD : "computeRMSD(d5, d6, l_atomParam)

	except :
		print "ERROR : RMSD calculation failed."
