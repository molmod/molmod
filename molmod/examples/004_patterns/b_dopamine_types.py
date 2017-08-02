#!/usr/bin/env python

from __future__ import print_function

from molmod import *

# 0) Load the molecule and set the default graph
mol = Molecule.from_file("dopamine.xyz")
mol.set_default_graph()

# 1) Define the atom types
atom_types = {
    "csp3": CritAnd(HasAtomNumber(6), HasNumNeighbors(4)),
    "csp2": CritAnd(HasAtomNumber(6), HasNumNeighbors(3)),
    "hoh": CritAnd(HasAtomNumber(1), HasNeighborNumbers(8)),
}

# 2) loop to detect all atom types
detected = {}
for i in range(mol.size):
    # match label is going to be the label of the matching atom type.
    match_label = None
    for label, atom_type in atom_types.items():
        if atom_type(i, mol.graph):
            # We have a match.
            if match_label is None:
                # This is how it should be.
                match_label = label
            else:
                # This should not happen after a previous match.
                raise TypeError("Atom %i matches more than one type, at least %s and %s" % (i, label, match_label))
    # Get the list of atom indexes associated with match_label. If such a list
    # is not present in the detected dictionary yet, an empty list is assigned
    # to the match_label key and returned. (Lookup setdefault in the Python
    # docs if this is not clear.)
    l = detected.setdefault(match_label, [])
    # Add the current atom index to the list.
    l.append(i)

# 3) print out the detected atom types
for label, indexes in detected.items():
    print(label, indexes)
