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

# 0) Create a molecule object based on the XYZ file 'caffeine.xyz'.
mol = Molecule.from_file("caffeine.xyz")

# 1) Newly created molecule objects do not have a graph, unless bond information
# is present in the data file. XYZ files only contain atomic elements and
# the corresponding coordinates. Therefore, the following line will not raise
# an error:
assert(mol.graph is None)

# 2) A molecular graph object can be derived from the geometry using a few rules
# of thumb and a database of well-known bond lengths. Such a routine is
# implemented in MolMod. It works as follows:
mol.set_default_graph()
# There are also other ways to define graphs with more control over the rules of
# thumb that detect the bonded atom pairs, e.g.
#
#   mol.graph = MolecularGraph.from_geometry(scaling=1.5)
#
# will also detect the breaking bond in a transition state.

# 3) Print all edges, i.e. bonds in the graph. The edges list is ordered and
# each item is a frozenset with two elements to stress the undericted nature of
# the molecular graphs.
print "All edges (or bonds)"
print mol.graph.edges
print
# Print the third bond.
print "The third bond:", mol.graph.edges[2]
# It is not possible to access only one of the two atom indexes of an edge. The
# following won't work because a frozenset is like unordered list.
#
#   print graph.edges[2][0]
#
# One can get both indexes of an edge at the same time:
i, j = mol.graph.edges[2]
print "The indexes of the second bond:", i, "and", j
# It is not possible to know a priori which number is assigned to i and which to
# j.

# 4) Print all the atom indexes of all C-H bonds on screen. Note that counting
# starts from zero.
# The loop relies on the 'unpack' feature in Python. It is equivalent to
#
#   for edge in mol.graph.edges:
#       i, j = edge
#
# but is a little more compact.
for i, j in mol.graph.edges:
    if mol.numbers[i] == 6 and mol.numbers[j] == 1:
        print "C-H bond:", i, j
    elif mol.numbers[j] == 6 and mol.numbers[i] == 1:
        print "C-H bond:", j, i
