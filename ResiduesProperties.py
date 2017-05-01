#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-


def hydrophobicResidue(res):
	hydrophobicRes = ["ALA", "PHE", "GLY", "ILE", "LEU", "MET", "PRO", "VAl"]
	if res.upper() in hydrophobicRes:
		return True
	else:
		return False
