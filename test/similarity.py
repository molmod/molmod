# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from molmod.molecular_graphs import MolecularGraph
from molmod.units import angstrom
from molmod.similarity import *

from ccio.xyz import XYZFile

import unittest, numpy


__all__ = ["SimilarityTestCase"]


class SimilarityTestCase(unittest.TestCase):
    def get_molecules(self):
        tpa = XYZFile("input/tpa.xyz").get_molecule()
        tpa.title = "tpa"
        graph = MolecularGraph(tpa)
        tpa1 = graph.randomized_molecule(1000, 0.1, 0.2, 0.07)
        tpa1.title = "tpa1"
        tpa2 = graph.randomized_molecule(1000, 0.2, 0.3, 0.1)
        tpa2.title = "tpa2"
        tea = XYZFile("input/tea.xyz").get_molecule()
        tea.title = "tea"
        water = XYZFile("input/water.xyz").get_molecule()
        water.title = "water"
        cyclopentane = XYZFile("input/cyclopentane.xyz").get_molecule()
        cyclopentane.title = "cyclopentane"

        return [tpa, tpa1, tpa2, tea, water, cyclopentane]

    def test_all_distances(self):
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.distances = all_distances(molecule)
        self.do_test(molecules, margin=0.1)

    def test_cutoff_distances(self):
        radius = 4*angstrom
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.distances = cutoff_distances(molecule, radius)
        self.do_test(molecules, 0.05, radius)

    def do_test(self, molecules, margin, radius=None, verbose=False):
        for molecule in molecules:
            molecule.norm = numpy.sqrt(calculate_similarity(molecule.distances, molecule.distances, margin, radius))
        if verbose:
            print
            print " "*14, "".join("%15s" % molecule.title for molecule in molecules)
        for index1, molecule1 in enumerate(molecules):
            if verbose: print "%15s" % molecule1.title,
            for index2, molecule2 in enumerate(molecules):
                similarity = calculate_similarity(molecule1.distances, molecule2.distances, margin, radius)/molecule1.norm/molecule2.norm
                if verbose: print ("%14.3f" % similarity),
            if verbose: print


