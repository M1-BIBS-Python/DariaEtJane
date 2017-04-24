from structureTools import parserPDB, preparePDB
import sys
from path import path
import re

pathfile = raw_input("Entrez le chemin du dossier : ")

for files in path(pathfile).walkfiles():
	if re.match(pathfile+"1BRS(.)*", files):
		dicoPDB = parserPDB(files)
		newdicoPDB = preparePDB(dicoPDB)
	
print dicoPDB
