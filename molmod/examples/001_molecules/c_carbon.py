#!/usr/bin/env python

from molmod import *

# 0) Load the molecule
mol = Molecule.from_file("ibuprofen.sdf")

# 1) Print out all carbon positions
for i in range(mol.size):
    if mol.numbers[i] == 6:
        print(mol.coordinates[i]/angstrom)
