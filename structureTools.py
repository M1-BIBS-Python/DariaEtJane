import os
import string
import math
import numpy
from path import path

def distancePoints((x1,y1,z1),(x2,y2,z2)):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates 
       output: distance"""
    x = (x1-x2)
    y = (y1-y2)
    z = (z1-z2)
    return math.sqrt(x*x+y*y+z*z)

def scorelist(listFiles, listScores):
	filename="Scoring.txt"
	pathname="/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/DariaEtJane/scoring_Cornell"
	
	if not os.path.exists(pathname):
		os.makedirs(pathname)
	
	filepath = os.path.join(pathname, filename)
	fileid = open(filepath, 'w+')
	
	xarray= numpy.array(listFiles)
	yarray= numpy.array(listScores)
	
	data= numpy.array([xarray,yarray])
	data= data.T
	data= data[-data[:, 1].argsort()]
	
	numpy.savetxt(fileid,data,fmt=['%d','%d'])
	fileid.close()
	
	return data[0][0]

def writePDB(dPDBrec, dPDBlig, prediction = "complexe_predit_score1.pdb") :

    pred = open(prediction, "w")

    for chain in dPDBrec["nchaine"]:
        for res in dPDBrec[chain]["position"] :
            for atom in dPDBrec[chain][res]["atome"] :
				pred.write("ATOM  %5d  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDBrec[chain][res][atom]["ID"], atom, dPDBrec[chain][res]["residu"],chain, res,dPDBrec[chain][res][atom]["x"], dPDBrec[chain][res][atom]["y"],dPDBrec[chain][res][atom]["z"] ))
	
	for chain in dPDBlig["nchaine"]:
		for res in dPDBlig[chain]["position"] :
			for atom in dPDBlig[chain][res]["atome"] :
				pred.write("ATOM  %5d  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDBlig[chain][res][atom]["ID"], atom, dPDBlig[chain][res]["residu"],chain, res,dPDBlig[chain][res][atom]["x"], dPDBlig[chain][res][atom]["y"],dPDBlig[chain][res][atom]["z"] ))
                    
    pred.close()
