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

# 0) Load the molecule and set the default graph
mol = Molecule.from_file("dopamine.xyz")
mol.set_default_graph()

# 1) Build a list of atom indexes involved in angles.
angles = []
# First loop over all atoms on the molecule.
for i1 in xrange(mol.size):
    # For each atom we will find all bending angles centered at the current
    # atom. For this we construct (an ordered!) list of all bonded neighbors.
    n = list(mol.graph.neighbors[i1])
    # The second loop iterates over all neighbors. The enumerate function is
    # used to assign a counter value to the variable index.
    for index, i0 in enumerate(n):
        # The third loop iterates over all other neighbors that came before i1.
        for i2 in n[:index]:
            # Each triple is stored as an item in the list angles.
            angles.append((i0, i1, i2))

# 2) Iterate over all angles, compute and print.
print "An overview of all bending angles in dopamine:"
for i0, i1, i2 in angles:
    # Notice again the [0] at the end.
    angle = bend_angle(mol.coordinates[[i0, i1, i2]])[0]
    # Python formatting of the indexes, symbols, and the angle in degrees.
    print "%2i %2i %2i    %2s %2s %2s    %5.1f" % (
        i0, i1, i2, mol.symbols[i0], mol.symbols[i1], mol.symbols[i2], angle/deg
    )
