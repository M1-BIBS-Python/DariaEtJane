from structureTools import parserPDB, preparePDB, scorelist
from RunScore import compEner
import sys
from path import path
import re

#pathfile = raw_input("Entrez le chemin du dossier : ")
dico_Rec = preparePDB(sys.argv[1])
dico_Scores = {}
n=0

for files in path("/home/kazevedo/Documents/M1BIBS/S2/Python/Python2016/Docking/sample").walkfiles():
	if re.match("/home/kazevedo/Documents/M1BIBS/S2/Python/Python2016/Docking/sample/1BRS(.)*", files):
		dico_confH = preparePDB(files)
		dico_Scores[files]= compEner(dico_Rec, dico_confH, "B", "D")
		dico_confH = {}

scorelist(dico_Scores)
bestScore = float(max(dico_Scores.values()))

print bestScore
