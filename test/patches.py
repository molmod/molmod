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


from molmod.patches import *
from molmod.io.xyz import XYZFile
from molmod.molecular_graphs import MolecularGraph
from molmod.units import deg, angstrom

import unittest, numpy


__all__ = ["PatchTestCase"]


class PatchTestCase(unittest.TestCase):
    def test_methyl(self):
        patch = MethylPatch()
        mol = XYZFile("input/water.xyz").get_molecule()
        graph = MolecularGraph.from_geometry(mol)
        new_mol, new_graph = patch.execute(mol, graph, [2], theta=109.45*deg, bond=1.54*angstrom)
        new_mol.write_to_file("output/methanol.xyz")
