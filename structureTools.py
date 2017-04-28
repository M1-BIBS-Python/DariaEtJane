import os
import string
import math
import numpy
from path import path
from ForceField3 import chargePDB, epsilon_vdw_PDB

dcharge = chargePDB()
dvdw, depsilon = epsilon_vdw_PDB()

def parserPDB(infile):
	
	fichier = open(infile,"r")
	lines = fichier.readlines() #Chaque element est une ligne du fichier

	d_proteine= {}
	d_proteine["nchaine"]=[]

	for i in lines:
		if i[0:4]=="ATOM": #On verifie qu'on est bien sur la ligne d'un atome
			
			nchaine=i[21] #On enregistre le numero de chaine
			#On verifie que la chaine n'est pas deja enregistree dans le dictionnaire pour ne pas ecraser l'information
			if not nchaine in d_proteine["nchaine"]:
				d_proteine["nchaine"].append(nchaine)
				d_proteine[nchaine]={}
				d_proteine[nchaine]["position"]=[]
			
			#Pour chaque chaine, on enregistre les positions
			position=i[22:26].strip()
			if not position in d_proteine[nchaine]["position"]:
				d_proteine[nchaine]["position"].append(position)
				d_proteine[nchaine][position]={}
				d_proteine[nchaine][position]["residu"]=i[17:20].strip()
				d_proteine[nchaine][position]["atome"]=[]
			
			#Pour chaque position on enregistre les atomes
			atome=i[12:16].strip()
			if not atome in d_proteine[nchaine][position]["atome"]:
				d_proteine[nchaine][position]["atome"].append(atome)
				d_proteine[nchaine][position][atome]={}
			
			#Pour chaque atome, on enregistre ses coordonnees et son ID
			d_proteine[nchaine][position][atome]["x"]=float(i[30:38])
			d_proteine[nchaine][position][atome]["y"]=float(i[38:46])
			d_proteine[nchaine][position][atome]["z"]=float(i[46:54])
			d_proteine[nchaine][position][atome]["ID"]=float(i[6:11])
	fichier.close()	
	return d_proteine

def assignParams(dPDB, dcharge, dvdw, depsilon):

    first = True
    
    for chain in dPDB["nchaine"]:
		for resi in dPDB[chain]["position"] :
        
			# means this is the Nter residue
			if first and (dPDB[chain][resi].has_key("H1") or dPDB[chain][resi].has_key("H2") or dPDB[chain][resi].has_key("H3")) :

					for atomi in dPDB[chain][resi]["atome"] :
						if atomi == "H1" :
							dPDB[chain][resi]["H1"]["charge"] = 0.1984
							dPDB[chain][resi]["H1"]["vdw"] = 0.6
							dPDB[chain][resi]["H1"]["epsilon"] = 0.0157


						elif atomi == "H2" :
							dPDB[chain][resi]["H2"]["charge"] = 0.1984
							dPDB[chain][resi]["H2"]["vdw"] = 0.6
							dPDB[chain][resi]["H2"]["epsilon"] = 0.0157

						elif atomi == "H3" :
							dPDB[chain][resi]["H3"]["charge"] = 0.1984   
							dPDB[chain][resi]["H3"]["vdw"] = 0.6
							dPDB[chain][resi]["H3"]["epsilon"] = 0.0157

						elif atomi == "N" :
							dPDB[chain][resi]["N"]["charge"] = 0.1592
							dPDB[chain][resi]["N"]["vdw"] = 1.875
							dPDB[chain][resi]["N"]["epsilon"] =  0.17

						elif atomi == "CA" :
							dPDB[chain][resi]["CA"]["charge"] = 0.0221
							dPDB[chain][resi]["CA"]["vdw"] = 1.9080
							dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

						elif atomi == "HA" :
							dPDB[chain][resi]["HA"]["charge"] = 0.116
							dPDB[chain][resi]["HA"]["vdw"] = 1.1
							dPDB[chain][resi]["HA"]["epsilon"] = 0.0157

						else:
							dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
							dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
							dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]
                        
                    
					first = False
                
			# means this is the Cter residue   
			elif first == False and (dPDB[chain][resi].has_key("OXT")):

				for atomi in dPDB[chain][resi]["atome"] :

					if atomi == "CA" :
						dPDB[chain][resi]["CA"]["charge"] = -0.2493
						dPDB[chain][resi]["CA"]["vdw"] = 1.9080
						dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

					elif atomi == "C" :
						dPDB[chain][resi]["C"]["charge"] = 0.7231
						dPDB[chain][resi]["C"]["vdw"] = 1.9080
						dPDB[chain][resi]["C"]["epsilon"] = 0.0860
                    
					elif atomi == "O" :
						dPDB[chain][resi]["O"]["charge"] = -0.7855
						dPDB[chain][resi]["O"]["vdw"] = 1.6612
						dPDB[chain][resi]["O"]["epsilon"] = 0.2100

					elif atomi == "OXT" :
						dPDB[chain][resi]["OXT"]["charge"] = -0.7855
						dPDB[chain][resi]["OXT"]["vdw"] = 1.6612
						dPDB[chain][resi]["OXT"]["epsilon"] = 0.2100
                    
					else:
						dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
						dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
						dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]

        
			else :
				# for all the other residues
				for atomi in dPDB[chain][resi]["atome"] :

					# particular case of Histidine
					if dPDB[chain][resi]["residu"] == "HIS" and dPDB[chain][resi].has_key("HD1") :
						# means this is a HID
						dPDB[chain][resi]["residu"] = "HID"
                    
					elif dPDB[chain][resi]["residu"] == "HIS" and dPDB[chain][resi].has_key("HE2") :
						# means this is a HIE
						dPDB[chain][resi]["residu"] = "HIE"

					# general case                    
					dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
					dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
					dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]


def preparePDB(infile) :
    """
    assigning all the atomic params (e.g. charges, vdw and epsilon) for all the atoms
    of the infile PDB
    Input: filename of the pdb
    Output: dico dPDB with 3D coords and charges, vdw and epsilon params for each atom
    """
    
    dPDB = parserPDB(infile)
    assignParams(dPDB, dcharge, dvdw, depsilon)

    return dPDB

def distancePoints((x1,y1,z1),(x2,y2,z2)):
    """Computes the distance between the two sets of coordinates
       input: 2 tuples with the corresponding coordinates 
       output: distance"""
    x = (x1-x2)
    y = (y1-y2)
    z = (z1-z2)
    return math.sqrt(x*x+y*y+z*z)

def scorelist(dPDB):
	n=0
	listeN=[]
	print dPDB.values()
	
	filename="Scoring.txt"
	pathname="/home/kazevedo/Documents/M1BIBS/S2/Python/GitRepo/DariaEtJane/scoring_Cornell"
	
	if not os.path.exists(pathname):
		os.makedirs(pathname)
	
	filepath = os.path.join(pathname, filename)
	fileid = open(filepath, 'w+')
	
	i=0
	while i<len(dPDB.values()):
		n+=1
		listeN.append(n)
		print listeN
		i+=1
	
	xarray= numpy.array(listeN)
	yarray= numpy.array(dPDB.values())
	
	data= numpy.array([xarray,yarray])
	data= data.T
	data= numpy.sort(data, axis=0)
	
	numpy.savetxt(fileid,data,fmt=['%d','%d'])
	
	fileid.close()

def writePDB(dPDB, prediction = "complexe_predit_score1.pdb") :

    pred = open(prediction, "w")

    for chain in dPDB["nchaine"]:
        for res in dPDB[chain]["position"] :
            for atom in dPDB[chain][res]["atome"] :
				pred.write("ATOM  %5s  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDB[chain][res][atom]["id"], atom, dPDB[chain][res]["residu"],chain, res,dPDB[chain][res][atom]["x"], dPDB[chain][res][atom]["y"],dPDB[chain][res][atom]["z"] ))
                    
    pred.close()
