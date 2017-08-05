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


import numpy as np
import pkg_resources

from molmod.test.common import BaseTestCase
from molmod.io import *
from molmod import *


__all__ = ["FCHKTestCase"]


class FCHKTestCase(BaseTestCase):
    def test_fchk(self):
        expected_coordinates = np.array([
            -1.22166133E-01, -7.76846132E-02,  2.18269387E-02,  2.73897634E-04,
            -7.83119835E-04,  3.12916345E+00,  2.85013528E+00,  7.51246138E-02,
             4.37365470E+00, -1.32111802E+00, -2.63535370E+00,  4.12342635E+00,
            -1.40571530E+00,  2.63466892E+00,  3.99680164E+00, -2.75049837E+00,
            -3.28207067E+00,  3.21191373E+00,  7.00973573E-02,  1.48008009E+00,
            -8.88472666E-01, -1.07582390E+00,  3.28402757E+00,  5.65864900E+00,
             3.76723726E+00, -1.48292907E+00,  4.52735777E+00,
        ])

        fchk = FCHKFile(pkg_resources.resource_filename(__name__, "../../data/test/1TOH.b3lyp.fchk"))
        self.assertEqual(fchk.title, "opt")
        self.assertEqual(fchk.command, "FOpt")
        self.assertEqual(fchk.lot, "RB3LYP")
        self.assertEqual(fchk.basis, "6-311+G(d,p)")
        self.assertEqual(fchk.fields["Number of atoms"], 9)
        self.assertEqual(fchk.fields["Number of independant functions"], 142)
        self.assertAlmostEqual(fchk.fields["Virial Ratio"], 2.002408027154329)
        self.assertArraysAlmostEqual(fchk.fields["Current cartesian coordinates"], expected_coordinates)
        self.assert_(isinstance(fchk.fields["Shell to atom map"][0], np.integer))
        self.assertEqual(fchk.fields["Shell to atom map"][-1], 9)

        fchk = FCHKFile(pkg_resources.resource_filename(__name__, "../../data/test/1TOH.b3lyp.trim.fchk"), ignore_errors=True)
        self.assertEqual(fchk.title, "opt")
        self.assertEqual(fchk.command, "FOpt")
        self.assertEqual(fchk.lot, "RB3LYP")
        self.assertEqual(fchk.basis, "6-311+G(d,p)")
        self.assertEqual(fchk.fields["Number of atoms"], 9)
        self.assertEqual(fchk.fields["Number of independant functions"], 142)
        self.assertAlmostEqual(fchk.fields["Virial Ratio"], 2.002408027154329)

        fchk = FCHKFile(pkg_resources.resource_filename(__name__, "../../data/test/1TOH.b3lyp.trim.fchk"), ignore_errors=True, field_labels=["Virial Ratio"])
        self.assertAlmostEqual(fchk.fields["Virial Ratio"], 2.002408027154329)
