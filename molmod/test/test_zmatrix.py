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
from molmod import *



__all__ = ["ZMatrixTestCase"]


class ZMatrixTestCase(BaseTestCase):
    def test_constency(self):
        def test_one(xyz_fn, checks, reorder=False):
            mol = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/" + xyz_fn))
            mol.set_default_graph()
            zmat_gen = ZMatrixGenerator(mol.graph)
            if reorder is False:
                self.assertArraysEqual(zmat_gen.new_index, np.arange(mol.size))
                self.assertArraysEqual(zmat_gen.old_index, np.arange(mol.size))

            zmat0 = zmat_gen.cart_to_zmat(mol.coordinates)
            for field, index, value in checks:
                self.assertAlmostEqual(zmat0[field][index], value, 2, "%s:%i %f!=%f" % (field,index,zmat0[field][index],value))

            numbers0, coordinates0 = zmat_to_cart(zmat0)
            mol0 = Molecule(numbers0, coordinates0)
            #mol0.write_to_file("zmat_%s" % os.path.basename(xyz_fn))
            mol0.set_default_graph()
            zmat_gen0 = ZMatrixGenerator(mol0.graph)
            self.assertArraysEqual(zmat_gen0.new_index, np.arange(mol.size))
            self.assertArraysEqual(zmat_gen0.old_index, np.arange(mol.size))

            zmat1 = zmat_gen0.cart_to_zmat(mol0.coordinates)
            for field, index, value in checks:
                self.assertAlmostEqual(zmat1[field][index], value, 2, "%s:%i %f!=%f" % (field,index,zmat1[field][index],value))

            numbers1, coordinates1 = zmat_to_cart(zmat1)

            self.assertArraysEqual(numbers0, numbers1)
            self.assertArraysAlmostEqual(coordinates0, coordinates1, 1e-5)

        checks = [
            ("number",   0, 7),
            ("number",   1, 6),
            ("distance", 1, 1.440*angstrom),

            ("number",   2, 6),
            ("distance", 2, 1.450*angstrom),
            ("angle",    2, 109.471*deg),

            ("number",   3, 6),
            ("distance", 3, 1.450*angstrom),
            ("angle",    3, 109.471*deg),
            ("dihed",    3, 180*deg),

            ("number",   4, 6),
            ("distance", 4, 1.440*angstrom),
            ("rel1",     4, 4),
            ("angle",    4, 109.471*deg),
            ("rel2",     4, 3),
            ("dihed",    4, 180*deg),
            ("rel3",     4, 2),

            ("number",   7, 6),
            ("distance", 7, 1.440*angstrom),
            ("rel1",     7, 7),
            ("angle",    7, 109.471*deg),
            ("rel2",     7, 3),
            ("dihed",    7, -120*deg),
            ("rel3",     7, 6),

            ("number",   19, 1),
            ("distance", 19, 1.090*angstrom),
            ("rel1",     19, 16),
            ("angle",    19, 109.471*deg),
            ("rel2",     19, 1),
            ("dihed",    19, 120*deg),
            ("rel3",     19, 2),
        ]
        test_one("tpa.xyz", checks)
        test_one("thf_single.xyz", [], reorder=True)
