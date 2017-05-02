# DariaEtJane
Projet Docking de DE AZEVEDO KEVIN et DUHAMEL MARINE

Date: 2/05/2017

Used with Python 2.7.13 
OS: Linux
IDE: None used

Execution: Use command: python MainDocking.py -pdb1 Rec_natif_DP.pdb -pdb2 [directory] -pdb3 cplx_natif.pdb -pdb4 Lig_natif_DP_aligned.pdb -pdb5 Rec_natif_DP.pdb
Replace [directory] with the path of the directory containing the 948 ligand conformations.

Output: The program will create in your work directory a directory named "scoring_Cornell" which will contain two files "Scoring1.txt" and "Scoring2.txt", as well as two files "complexe_predit_score1.pdb" and "complexe_predit_score2.pdb" in the work directory

Python libraries used : before you run this program, make sure you have the following libraries installed
* itertools
math
numpy
os
path
string
sys

Note : This program is very time consuming. On a computer with an Intel(R) Celeron(R) CPU 1007U @ 1.50Hz, this program took approximately 5 hours to execute. Depending on your CPU configuration, it may take less or more time.

