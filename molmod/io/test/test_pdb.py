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

import pkg_resources

from molmod.test.common import *
from molmod.io import *
from molmod import *


__all__ = ["PDBTestCase"]


class PDBTestCase(BaseTestCase):
    def test_load_pdb(self):
        mol = load_pdb(pkg_resources.resource_filename(__name__, "../../data/test/il2.pdb"))
        self.assertEqual(mol.numbers[0], 7)
        self.assertEqual(mol.numbers[5], 1)
        self.assertEqual(mol.numbers[-1], 8)
        self.assertAlmostEqual(mol.coordinates[0,0]/angstrom, 17.166)
        self.assertAlmostEqual(mol.coordinates[14,1]/angstrom, -3.746)
        self.assertAlmostEqual(mol.coordinates[-1,2]/angstrom, 3.311)
        self.assertAlmostEqual(mol.occupancies[4], 1.0)
        self.assertAlmostEqual(mol.occupancies[-1], 1.0)
        self.assertEqual(len(mol.occupancies), mol.size)
        self.assertAlmostEqual(mol.betas[4], 45.70)
        self.assertAlmostEqual(mol.betas[-1], 64.30)
        self.assertEqual(len(mol.betas), mol.size)

    def test_dump_pdb(self):
        with tmpdir(__name__, 'test_dump_pdb') as dn:
            mol = Molecule([8,1,1], [[-1,0,0],[0,0,0],[1.1,0,0]])
            dump_pdb(
                "%s/test.pdb" % dn, mol, ["OT", "HT", "HT"],
                resnames=["HOH", "HOH", "HOH"], chain_ids=["A", "A", "A"],
                occupancies=[1.0, 0.5, 1.0], betas=[41.2, 78.1, 0.25]
            )
            with open("%s/test.pdb" % dn) as f:
                content = "".join(f)
            expected_content = "ATOM      1  OT  HOH A   2      -0.529   0.000   " \
                               "0.000  1.00 41.20          O   \nATOM      2  HT " \
                               " HOH A   2       0.000   0.000   0.000  0.50 78.1" \
                               "0          H   \nATOM      3  HT  HOH A   2      " \
                               " 0.582   0.000   0.000  1.00  0.25          H   \n"
            self.assertEqual(content, expected_content)
            mol_bis = load_pdb("%s/test.pdb" % dn)
            self.assertArraysAlmostEqual(mol.coordinates, mol_bis.coordinates, 1e-3)
            self.assertArraysAlmostEqual(mol.numbers, mol_bis.numbers)
