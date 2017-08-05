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
# --


from __future__ import division

import unittest

from molmod import *
from molmod.periodic import periodic
from molmod.isotopes import ame2003, nubtab03


__all__ = ["DataTestCase"]


class DataTestCase(unittest.TestCase):
    def test_periodic(self):
        # assert that the most important fields have been loaded
        # and that the conversions are done correctly
        self.assertEqual(periodic[1].name, "Hydrogen")
        self.assertEqual(periodic[1].symbol, "H")
        self.assertEqual(periodic[1].number, 1)
        self.assertEqual(periodic[1].row, 1)
        self.assertEqual(periodic[1].col, 1)
        self.assertAlmostEqual(periodic[1].mass, 1.00794*unified, 3)
        self.assertAlmostEqual(periodic[1].density, 0.12411, 3)
        self.assertAlmostEqual(periodic[1].covalent_radius/angstrom, 0.31, 3)
        self.assertAlmostEqual(periodic[1].vdw_radius/angstrom, 1.09, 3)
        self.assertEqual(periodic[1].artificial, False)
        self.assertAlmostEqual(periodic[1].red, 1.0, 3)
        self.assertAlmostEqual(periodic[1].green, 1.0, 3)
        self.assertAlmostEqual(periodic[1].blue, 1.0, 3)

        self.assertEqual(periodic[1], periodic["H"])
        self.assertEqual(periodic[1], periodic["h"])

    def test_ame2003(self):
        self.assertAlmostEqual(ame2003.masses[6][12]/12/unified, 1.0)
        self.assertAlmostEqual(ame2003.masses[1][1]/unified, 1.00782503207)

    def test_nubtab03(self):
        self.assertAlmostEqual(nubtab03.abundances[1][1], 99.9885)
