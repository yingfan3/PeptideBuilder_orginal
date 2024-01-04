"""
Simple example script demonstrating how to use the PeptideBuilder library.

The script generates a peptide consisting of six arginines in alpha-helix
conformation, and it stores the peptide under the name "example.pdb".
"""

from PeptideBuilder import Geometry
import PeptideBuilder
import random

geo = Geometry.geometry("G")
AA = ["G", "A", "S", "C", "V", "I", "L", "T", "R", "K", "D", "E", "N", "Q", "M", "H", "P", "F", "Y", "W"]

SG_chain = []
CD1_chain = []
letter = random.choice(AA)
geo1 = Geometry.geometry(letter)
SG_chain.append(letter)
letter = random.choice(AA)
geo2 = Geometry.geometry(letter)
CD1_chain.append(letter)
structure = PeptideBuilder.initialize_res(geo, geo1, geo2)
for i in range(7):
    letter = random.choice(AA)
    geo1 = Geometry.geometry(letter)
    SG_chain.append(letter)
    letter = random.choice(AA)
    geo2 = Geometry.geometry(letter)
    CD1_chain.append(letter)
    PeptideBuilder.add_residue(structure, geo, geo1, geo2)

# add terminal oxygen (OXT) to the final glycine
# PeptideBuilder.add_terminal_OXT(structure)
sidechain = []
for i in range(len(SG_chain)):
    sidechain.append((SG_chain[i] + CD1_chain[i]))

import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("%s-%s-%s-%s-%s-%s-%s-%s.pdb" % (
    sidechain[0], sidechain[1], sidechain[2], sidechain[3], sidechain[4], sidechain[5],
    sidechain[6], sidechain[7]))




