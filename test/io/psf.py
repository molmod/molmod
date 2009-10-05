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


from molmod.io.psf import PSFFile
from molmod.molecules import Molecule

import numpy, unittest


__all__ = ["PSFTestCase"]


class PSFTestCase(unittest.TestCase):
    def test_load(self):
        psf = PSFFile("input/thf.psf")
        self.assert_(psf.bonds.shape[0] == 832)
        self.assert_(psf.bends.shape[0] == 1600)
        self.assert_(psf.dihedrals.shape[0] == 2112)
        g = psf.get_graph()

    def test_load_vmd(self):
        psf = PSFFile("input/pentapeptide.psf")
        self.assert_(psf.bonds.shape[0] == 44)
        self.assert_(psf.bends.shape[0] == 75)
        self.assert_(psf.dihedrals.shape[0] == 98)
        g = psf.get_graph()

    def test_dump(self):
        m = Molecule.from_file("input/thf.xyz")
        psf = PSFFile()
        psf.add_molecule(m)
        psf.write_to_file("output/thf.psf")

    def test_tetra(self):
        molecule = Molecule.from_file("input/tetra.xyz")
        psf = PSFFile()
        psf.add_molecule(molecule)
        self.assert_(psf.bonds.shape[0] == 4)
        self.assert_(psf.bends.shape[0] == 6)
        psf.write_to_file("output/tetra.psf")

    def test_many_separate(self):
        psf = PSFFile()
        molecule = Molecule.from_file("input/ethene.xyz")
        psf.add_molecule(molecule)
        psf.add_molecule(molecule)
        molecule = Molecule.from_file("input/tea.xyz")
        psf.add_molecule(molecule)
        psf.write_to_file("output/many_separate.psf")


