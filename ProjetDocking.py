from structureTools import parserPDB, preparePDB, scorelist, writePDB
from RunScore import compEner
import sys
from path import path
import re

#pathfile = raw_input("Path of the repertory containing all the conformations : ")
dico_Rec = preparePDB(sys.argv[1])
dico_Scores = {}
n=0
listFiles=[]
listScores=[]

for files in path("/home/kazevedo/Documents/M1BIBS/S2/Python/Python2016/Docking/sample").walkfiles():
	if re.match("/home/kazevedo/Documents/M1BIBS/S2/Python/Python2016/Docking/sample/1BRS(.)*", files):
		n+=1
		dico_confH = preparePDB(files)
		listFiles.append(n)
		listScores.append(compEner(dico_Rec, dico_confH, "B", "D"))
		dico_confH = {}

print listFiles
print listScores
scorelist(listFiles, listScores)
bestScore = float(max(listScores))
print bestScore

#pathfile = raw_input("Path of the PDB file giving the best score : ")

dico_Lig = preparePDB("/home/kazevedo/Documents/M1BIBS/S2/Python/Python2016/Docking/sample/1BRS_A_1BRS_B_allatom_1_DP.pdb")
writePDB(dico_Rec, dico_Lig, prediction = "complexe_predit_score1.pdb")

