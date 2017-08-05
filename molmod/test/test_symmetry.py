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


import pkg_resources

from molmod.test.common import BaseTestCase
from molmod import *


__all__ = ["SymmetryTestCase"]


class SymmetryTestCase(BaseTestCase):
    def test_rotsym_butane(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/butane.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph)
        self.assertEqual(rotsym, 2)

    def test_rotsym_cyclopentane(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/cyclopentane.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph)
        self.assertEqual(rotsym, 1)

    def test_rotsym_octane(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/octane.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph, threshold=0.01)
        self.assertEqual(rotsym, 2)

    def test_rotsym_tetra(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/tetra.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph, threshold=0.01)
        self.assertEqual(rotsym, 12)

    def test_rotsym_water(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/water.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph, threshold=0.01)
        self.assertEqual(rotsym, 2)

    def test_rotsym_benzene(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/benzene.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph, threshold=0.01)
        self.assertEqual(rotsym, 12)

    def test_rotsym_ethane(self):
        molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/ethane.xyz"))
        molecule.set_default_graph()
        rotsym = compute_rotsym(molecule, molecule.graph, threshold=0.01)
        self.assertEqual(rotsym, 6)
