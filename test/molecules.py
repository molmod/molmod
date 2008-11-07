# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
# --


from molmod.molecules import *

from molmod.io.xyz import XYZFile

import unittest, numpy


__all__ = ["MoleculesTestCase"]


class MoleculesTestCase(unittest.TestCase):
    def test_distance_matrix(self):
        molecule = XYZFile("input/tpa.xyz").get_molecule()
        dm = 0
        for i in 0,1,2:
            dm += numpy.subtract.outer(molecule.coordinates[:,i],molecule.coordinates[:,i])**2
        dm = numpy.sqrt(dm)
        self.assert_((abs(molecule.distance_matrix - dm) < 1e-5).all(), "Wrong distance matrix")


