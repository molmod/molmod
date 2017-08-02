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

# 0) Load the molecule.
mol = Molecule.from_file("ibuprofen.sdf")

# 1) Print the largest element in the distance matrix.
print "Largest interatomic distance [A]:"
print mol.distance_matrix.max()/angstrom
# Some comments:
#  - One can just write code that assumes the attribute distance_matrix is
#    present, but in fact it is only computed once the first time it is used.
#  - The method max is a method of the numpy array. It returns the largest value
#    in an array.
