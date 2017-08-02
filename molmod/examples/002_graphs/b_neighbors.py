#!/usr/bin/env python

from __future__ import print_function

from molmod import *

# 0) Create a molecule object based on the XYZ file 'caffeine.xyz'. Also
# initialize the graph.
mol = Molecule.from_file("caffeine.xyz")
mol.set_default_graph()
# The neighbors attribute is a dictionary with atom indexes as keys and a set
# of neighboring atom indexes as corresponding values. The dictionary will be
# constructed as soon as it is first accessed, e.g.
print("All neighbors:")
print(mol.graph.neighbors)
print()

# 1) Print the atomic number of the third atom and the atomic numbers of its
# neighbors.
print("Symbol of third atom (should be oxygen):", mol.symbols[2])
print("The number of bonds to the third atom:", len(mol.graph.neighbors[2]))
print("Indexes of neigbors of third atom:", mol.graph.neighbors[2])
print("Symbols of neigbors of third atom:", [mol.symbols[i] for i in mol.graph.neighbors[2]])

# 2) Look up all methyl groups.
for i, ns in mol.graph.neighbors.items():
    if mol.numbers[i] == 6 and len(ns) == 4:
        # Get the indexes of the hydrogen atoms.
        h_indexes = [n for n in ns if mol.numbers[n]==1]
        # We must have three hydrogens
        if len(h_indexes) != 3:
            continue
        # Print stuff
        print("C=%i H=%i H=%i H=%i" % (i, h_indexes[0], h_indexes[1], h_indexes[2]))
