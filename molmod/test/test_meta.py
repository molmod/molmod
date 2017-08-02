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


import os
import unittest

import pkg_resources

from molmod import *



__all__ = ["MetaTestCase"]


class MetaTestCase(unittest.TestCase):
    def check_example(self, dirname, fn_py):
        root = pkg_resources.resource_filename(__name__, "../examples")
        print root
        assert os.path.isdir(root)
        cwd = os.getcwd()
        command = "cd %s/%s; PYTHONPATH=%s:${PYTHONPATH} ./%s 1> /dev/null 2> /dev/null" % (root, dirname, cwd, fn_py)
        print command
        retcode = os.system(command)
        self.assertEqual(retcode, 0)

    def test_example_000(self):
        self.check_example("000_units", "a_reaction.py")
        self.check_example("000_units", "b_chbond.py")
        self.check_example("000_units", "c_h2rot.py")

    def test_example_001(self):
        self.check_example("001_molecules", "a_convert.py")
        self.check_example("001_molecules", "b_com.py")
        self.check_example("001_molecules", "c_carbon.py")
        self.check_example("001_molecules", "d_size.py")
        self.check_example("001_molecules", "e_shape.py")

    def test_example_002(self):
        self.check_example("002_graphs", "a_graphs.py")
        self.check_example("002_graphs", "b_neighbors.py")
        self.check_example("002_graphs", "c_distances.py")
        self.check_example("002_graphs", "d_symmetries.py")

    def test_example_003(self):
        self.check_example("003_internal_coordinates", "a_bond_length.py")
        self.check_example("003_internal_coordinates", "b_bending_angles.py")
        self.check_example("003_internal_coordinates", "c_ff_hessian.py")
        self.check_example("003_internal_coordinates", "d_dft_hessian.py")

    def test_example_004(self):
        self.check_example("004_patterns", "a_propane_types.py")
        self.check_example("004_patterns", "b_dopamine_types.py")
