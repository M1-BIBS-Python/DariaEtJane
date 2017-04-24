import math
from ForceField3.py import chargePDB, epsilon_vdw_PDB

def computeAij(d_atomi, d_atomj) :

    Aij = math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),8))

    return Aij

def computeBij(d_atomi, d_atomj) :
    Bij = 2*math.sqrt(d_atomi["epsilon"]*d_atomj["epsilon"])*(pow((d_atomi["vdw"]+d_atomj["vdw"]),6))

    return Bij

