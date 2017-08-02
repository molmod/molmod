#!/usr/bin/env python

from __future__ import print_function

from molmod import *

# 0) Create a molecule object based on the XYZ file 'caffeine.xyz'. Also
# initialize the graph.
mol = Molecule.from_file("caffeine.xyz")
mol.set_default_graph()
# The distances attribute is a square matrix of integers with a row and column
# for every atom.
print("Number of atoms:", mol.size)
print("Shape of the distances array:", mol.graph.distances.shape)

# 1) The distances array can be used to get the 'minimal' number of bonds
# between two atoms. E.g. for two atoms in a five_membered ring, this is at most
# two:
print("Distance between atom 5 and 6 (part of a 5-ring):", mol.graph.distances[5,6])

# 2) The matrix can also be used to construct a mask of atom_pairs that are
# separated by at most three bonds.
mask = (mol.graph.distances <= 3) & (mol.graph.distances != 0)
# For a nice print output, the mask is converted to integers.
print("Mask of atom pairs that are typically involved in valence interactions:")
print(mask.astype(int))
