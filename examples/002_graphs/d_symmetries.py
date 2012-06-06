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

# 0) Create a molecule object based on the XYZ file 'ethanol.xyz'. Also
# initialize the graph.
mol = Molecule.from_file("ethanol.xyz")
mol.set_default_graph()

# 1) Print out all the isomorphisms. For the ethanol molecule, this comes down
# to all possible ways to rotate/mirror the methyl moieties, and combination
# thereof.
print "Isomorphisms in the form of one-to-one mappings"
for symmetry in mol.graph.symmetries:
    print symmetry

# 2) One can also request isomorphisms in the form of permutations, which is
# often more convenient:
print "Isomorphisms in the form of permutations"
for cycle in mol.graph.symmetry_cycles:
    print cycle
