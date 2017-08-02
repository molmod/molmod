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
# The distances attribute is a square matrix of integers with a row and column
# for every atom.
print "Number of atoms:", mol.size
print "Shape of the distances array:", mol.graph.distances.shape

# 1) The distances array can be used to get the 'minimal' number of bonds
# between two atoms. E.g. for two atoms in a five_membered ring, this is at most
# two:
print "Distance between atom 5 and 6 (part of a 5-ring):", mol.graph.distances[5,6]

# 2) The matrix can also be used to construct a mask of atom_pairs that are
# separated by at most three bonds.
mask = (mol.graph.distances <= 3) & (mol.graph.distances != 0)
# For a nice print output, the mask is converted to integers.
print "Mask of atom pairs that are typically involved in valence interactions:"
print mask.astype(int)
