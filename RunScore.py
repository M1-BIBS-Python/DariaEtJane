import math
from ForceField3 import chargePDB, epsilon_vdw_PDB
from structureTools import distancePoints

def computeAij(d_atomi, d_atomj) :

    Aij = math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),8))

    return Aij

def computeBij(d_atomi, d_atomj) :
    Bij = 2*math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),6))

    return Bij

def compEner(dPDB1, dPDB2, chain1, chain2) :

    nbresInter = 0
    #structureTools.initBfactor(dPDB)
    #print dPDB["chains"]
    
    # verifie que la chaine1 existe
    if not chain1 in dPDB1["nchaine"] :
        print "la chaine %s n'existe pas"%(chain1)
        sys.exit()
    # verifie que la chaine2 existe
    if not chain2 in dPDB2["nchaine"] :
        print "la chaine %s n'existe pas"%(chain2)
        sys.exit()
        
    print "dealing chains %s and %s"%(chain1, chain2)

    # DEBUG
    cmt = 0
    cmti = 0
    cmtj = 0
    Eij = 0
    
    # computes ener
    for resi in dPDB1[chain1]["position"] :
        #print "resi ", resi 
        for atomi in dPDB1[chain1][resi]["atome"] :
            cmti+=1
            cmtj = 0
            coordi = [dPDB1[chain1][resi][atomi]["x"], dPDB1[chain1][resi][atomi]["y"], dPDB1[chain1][resi][atomi]["z"]]

            for resj in dPDB2[chain2]["position"] :
                #print "resj" , resj
                for atomj in dPDB2[chain2][resj]["atome"] :
                    cmt+=1
                    cmtj+=1
                    coordj = [dPDB2[chain2][resj][atomj]["x"], dPDB2[chain2][resj][atomj]["y"], dPDB2[chain2][resj][atomj]["z"]]

                    # computes dist, Aij & Bij
                    
                    Rij = distancePoints(coordi, coordj)
                    #print "%s %s %s VS %s %s %s"%(chain1, resi, atomi, chain2, resj, atomj)
                    #print dPDB[chain1][resi][atomi]
                    #print dPDB[chain2][resj][atomj]
                    Aij = computeAij(dPDB1[chain1][resi][atomi], dPDB2[chain2][resj][atomj])
                    Bij = computeBij(dPDB1[chain1][resi][atomi], dPDB2[chain2][resj][atomj])
                    #print "%s %s %s VS %s %s %s= %s %s %s"%(chain1, resi, atomi, chain2, resj, atomj, Rij, Aij, Bij)
                    coulij =  (332.0522*dPDB1[chain1][resi][atomi]["charge"]*dPDB2[chain2][resj][atomj]["charge"])/(20*Rij)
                    #print "%s %s %s VS %s %s %s= Rij %s Aij %s Bij %s coulij %s"%(chain1, resi, atomi, chain2, resj, atomj, Rij, Aij, Bij, coulij)

                    # computes ENER
                    Etmp = float(Aij)/pow(Rij,8) - float(Bij)/pow(Rij, 6) + coulij
                    #print "Etmp ", Etmp
                    Eij = Eij + Etmp
                    
    #print cmti, cmtj, cmt

    return Eij
