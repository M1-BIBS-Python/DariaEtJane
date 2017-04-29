from structureTools import parserPDB, preparePDB, scorelist, writePDB
from RunScore import compEner
import sys
from path import path
import re

dico_Rec = preparePDB(sys.argv[1])
dico_Scores = {}
n=0
listFiles=[]
listScores=[]

pathdir = raw_input("Directory containing all the configurations : ")

for files in path(pathdir).walkfiles():
	if files.endswith(".pdb"):
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


filename = raw_input("Name of the PDB file giving the best score : ")
pathfile=[pathdir,'/',filename]

dico_Lig = preparePDB(''.join(pathfile))
writePDB(dico_Rec, dico_Lig, prediction = "complexe_predit_score1.pdb")

