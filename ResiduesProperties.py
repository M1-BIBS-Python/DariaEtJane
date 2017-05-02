#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-


def hydrophobicResidue(res):
	"""
	Test if the residue res is hydrophobic or not
	input : res the residue to test
	output : boolean True or False
	"""
	hydrophobicRes = ["ALA", "PHE", "GLY", "ILE", "LEU", "MET", "PRO", "VAl"]
	if res.upper() in hydrophobicRes:
		return True
	else:
		return False
