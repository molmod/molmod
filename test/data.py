# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of MolMod.
# 
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
# --


from molmod.data import periodic, bonds, BOND_SINGLE
from molmod.units import from_angstrom, from_unified, from_picometer

import unittest

__all__ = ["Data"]    


class Data(unittest.TestCase):
    # assert that the most important fields have been loaded
    # and that the conversions are done correctly

    def test_periodic(self):
        self.assertEqual(periodic[1].name, "Hydrogen")
        self.assertEqual(periodic[1].symbol, "H")
        self.assertEqual(periodic[1].number, 1)
        self.assertEqual(periodic[1].row, 1)
        self.assertEqual(periodic[1].col, 1)
        self.assertAlmostEqual(periodic[1].mass, from_unified(1.00794), 3)
        self.assertAlmostEqual(periodic[1].density, 0.763, 3)
        self.assertAlmostEqual(periodic[1].bond_radius, from_angstrom(0.3), 3)
        self.assertAlmostEqual(periodic[1].vdw_radius, 0.43, 3)
        self.assertEqual(periodic[1].valence, 1)
        self.assertEqual(periodic[1].artificial, False)
        self.assertAlmostEqual(periodic[1].red, 0.9, 3)
        self.assertAlmostEqual(periodic[1].green, 0.9, 3)
        self.assertAlmostEqual(periodic[1].blue, 0.9, 3)

        self.assertEqual(periodic[1], periodic["H"])
        self.assertEqual(periodic[1], periodic["h"])

    def test_periodic(self):
        self.assertAlmostEqual(bonds.lengths[BOND_SINGLE][frozenset([1,1])], from_picometer(74), 3)



