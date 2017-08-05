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


import unittest

import numpy as np
import pkg_resources

from molmod.test.common import *
from molmod.io import *
from molmod import *


__all__ = ["CMLTestCase"]


class CMLTestCase(unittest.TestCase):
    def test_consistency(self):
        molecules = [
            Molecule.from_file(pkg_resources.resource_filename(__name__, "../../data/test/cyclopentane.xyz")),
            Molecule.from_file(pkg_resources.resource_filename(__name__, "../../data/test/funny.xyz")),
        ]
        for m in molecules:
            m.set_default_graph()
        with tmpdir(__name__, 'test_consistency') as dn:
            dump_cml("%s/tmp.cml" % dn, molecules)
            check = load_cml("%s/tmp.cml" % dn)
        for m1, m2 in zip(molecules, check):
            self.assertEqual(m1.title, m2.title)
            self.assert_((m1.numbers==m2.numbers).all())
            self.assert_((m1.coordinates==m2.coordinates).all())
            self.assertEqual(m1.graph.num_vertices, m2.graph.num_vertices)
            self.assertEqual(set(m1.graph.edges), set(m2.graph.edges))

    def test_load(self):
        l = load_cml(pkg_resources.resource_filename(__name__, "../../data/test/1LJL_Cys10.cml"))
        self.assertEqual(l[0].title, "Kalium [+]")
        self.assertEqual(l[2].title, "Acetic acid [-]")
        self.assert_((l[2].numbers==np.array([6,6,8,8,1,1,1])).all())
        self.assertEqual(l[2].extra["charge"],"-1.0")
