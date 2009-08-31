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


from molmod.io.fchk import FCHKFile

import numpy

import unittest


__all__ = ["FCHKTestCase"]


class FCHKTestCase(unittest.TestCase):
    def test_fchk(self):
        fchk = FCHKFile("input/1TOH.b3lyp.fchk")
        fchk = FCHKFile("input/1TOH.b3lyp.trim.fchk", ignore_errors=True)
        #print fchk.molecule.numbers
        #print fchk.molecule.coordinates
        #print fchk.optimization_coordinates()
        #print fchk.optimized_molecule()







