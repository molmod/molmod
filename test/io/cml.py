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


from molmod.io.cml import *
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph

import numpy, unittest


__all__ = ["CMLTestCase"]


class CMLTestCase(unittest.TestCase):
    def test_consistency(self):
        molecules = [
            Molecule.from_file("input/cyclopentane.xyz"),
            Molecule.from_file("input/funny.xyz"),
        ]
        for m in molecules:
            m.set_default_graph()
        dump_cml("output/tmp.cml", molecules)
        check = load_cml("output/tmp.cml")
        for m1, m2 in zip(molecules, check):
            self.assertEqual(m1.title, m2.title)
            self.assert_((m1.numbers==m2.numbers).all())
            self.assert_((m1.coordinates==m2.coordinates).all())
            self.assertEqual(m1.graph.num_vertices, m2.graph.num_vertices)
            self.assertEqual(set(m1.graph.edges), set(m2.graph.edges))

    def test_load(self):
        l = load_cml("input/1LJL_Cys10.cml")
        self.assertEqual(l[0].title, "Kalium [+]")
        self.assertEqual(l[2].title, "Acetic acid [-]")
        self.assert_((l[2].numbers==numpy.array([6,6,8,8,1,1,1])).all())
        self.assertEqual(l[2].extra["charge"],"-1.0")


