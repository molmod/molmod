#!/usr/bin/env python
# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
#!/usr/bin/env python

from molmod import *

# 0) Create a molecule object based on the XYZ file 'caffeine.xyz'. Also
# initialize the graph.
mol = Molecule.from_file("caffeine.xyz")
mol.set_default_graph()
# The neighbors attribute is a dictionary with atom indexes as keys and a set
# of neighboring atom indexes as corresponding values. The dictionary will be
# constructed as soon as it is first accessed, e.g.
print "All neighbors:"
print mol.graph.neighbors
print

# 1) Print the atomic number of the third atom and the atomic numbers of its
# neighbors.
print "Symbol of third atom (should be oxygen):", mol.symbols[2]
print "The number of bonds to the third atom:", len(mol.graph.neighbors[2])
print "Indexes of neigbors of third atom:", mol.graph.neighbors[2]
print "Symbols of neigbors of third atom:", [mol.symbols[i] for i in mol.graph.neighbors[2]]

# 2) Look up all methyl groups.
for i, ns in mol.graph.neighbors.iteritems():
    if mol.numbers[i] == 6 and len(ns) == 4:
        # Get the indexes of the hydrogen atoms.
        h_indexes = [n for n in ns if mol.numbers[n]==1]
        # We must have three hydrogens
        if len(h_indexes) != 3:
            continue
        # Print stuff
        print "C=%i H=%i H=%i H=%i" % (i, h_indexes[0], h_indexes[1], h_indexes[2])
