def usage():
	usage = ("Computes the score involving the non-linked terms of the Cornell energy function between one receptor pdb file and a repository containing ligand solutions pdb files."
			"Then calculate the RMSD between the real complexe structure, and the best scored structure and the native contact number in the predicted structre.\n\n"
			"Input = \n"
			"\t-> -pdb1 a receptor pdb file\n"
			"\t-> -pdb2 a directory containing files of ligand solutions in pdb format\n"
			"\t-> -pdb3 the real 3D structure complexe pdb file\n"
			"\t-> -pdb4 the native ligand\n"
			"\t-> -pdb5 the native receptor\n\n"
			"Output = a text file with 3 columns and 4 row, respectively the names of the 2 pdb files used for the calculation and :\n"
			"\t->the corresponding RMSD on the entiere complexe\n"
			"\t->the corresponding RMSD on the ligand only\n"
			"\t->the calculation and the corresponding RMSD on the receptor only\n"
			"\t->the corresponding number of native contacts\n\n"
			"Obligatory argument : -pdb1,-pdb2,-pdb3, -pdb4, -pdb5 <path of the file/directory> -NbLigand <integer>\n"
			"\t-> absolute or relative path\n"	  
			"Optional argument : -atom <name of atom>\n"
			"\t-> the name of the atom used in the RMSD calculation separated with '_'\n"
			"\t-> default : CA, N, C, O\n"
			"\t-> to apply RMSD calculation on all the atoms, write 'ALL'\n")
	
	print usage
