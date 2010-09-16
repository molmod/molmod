# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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


from tamkin import *

import unittest, os, glob


__all__ = ["MetaTestCase"]


class MetaTestCase(unittest.TestCase):
    def check_example(self, dirname, fn_py):
        root = "../examples"
        self.assert_(os.path.isdir(root))
        retcode = os.system("cd %s; cd %s; ./%s 1> /dev/null 2> /dev/null" % (root, dirname, fn_py))
        self.assertEqual(retcode, 0)

    def test_example_002(self):
        self.check_example("002_graphs", "./a_graphs.py")
        self.check_example("002_graphs", "./b_neighbors.py")
        self.check_example("002_graphs", "./c_distances.py")
        self.check_example("002_graphs", "./d_symmetries.py")

    def test_code_quality(self):
        root = "../molmod"
        self.assert_(os.path.isdir(root))
        white = (" ", "\t")
        for fn in glob.glob("%s/*.py") + glob.glob("%s/io/*.py"):
            f = file(fn)
            for line in f:
                if line[-2] in white:
                    self.fail("Trailing whitespace in %s." % fn)
