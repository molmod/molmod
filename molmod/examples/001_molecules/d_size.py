#!/usr/bin/env python

from molmod import *

# 0) Load the molecule.
mol = Molecule.from_file("ibuprofen.sdf")

# 1) Print the largest element in the distance matrix.
print("Largest interatomic distance [A]:")
print(mol.distance_matrix.max()/angstrom)
# Some comments:
#  - One can just write code that assumes the attribute distance_matrix is
#    present, but in fact it is only computed once the first time it is used.
#  - The method max is a method of the numpy array. It returns the largest value
#    in an array.
