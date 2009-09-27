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


from molmod.toyff import guess_geometry
from molmod.io.xyz import XYZFile
from molmod.io.sdf import SDFReader
from molmod.molecular_graphs import MolecularGraph

import unittest, numpy, os


__all__ = ["ToyFFTestCase"]


class ToyFFTestCase(unittest.TestCase):
    def load_molecule(self, xyz_fn):
        molecule = XYZFile(os.path.join("input", xyz_fn)).get_molecule()
        molecule.graph = MolecularGraph.from_geometry(molecule)
        return molecule

    def iter_molecules(self, allow_multi=False):
        xyz_fns = [
          "water.xyz", "cyclopentane.xyz", "ethene.xyz", "funny.xyz",
          "tea.xyz", "tpa.xyz", "thf_single.xyz", "precursor.xyz",
          "butane.xyz", "octane.xyz",
        ]
        for xyz_fn in xyz_fns:
            molecule = self.load_molecule(xyz_fn)
            if allow_multi or len(molecule.graph.independent_vertices) == 1:
                yield molecule
        sdf_fns = [
            "example.sdf", "CID_22898828.sdf", "SID_55127927.sdf",
            "SID_56274343.sdf", "SID_40363570.sdf", "SID_40363571.sdf",
            "SID_31646548.sdf", "SID_31646545.sdf", "SID_41893278.sdf",
            "SID_41893280.sdf", "SID_54258192.sdf", "SID_55488598.sdf",
        ]
        for sdf_fn in sdf_fns:
            for i, molecule in enumerate(SDFReader(os.path.join("input", sdf_fn))):
                if allow_multi or len(molecule.graph.independent_vertices) == 1:
                    yield molecule

    def test_guess_geometry(self):
        for input_mol in self.iter_molecules(allow_multi=False):
            output_mol = guess_geometry(input_mol.graph)
            output_mol.title = input_mol.title
            output_mol.write_to_file("output/guess_%s.xyz" % input_mol.title)

