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


from molmod.molecular_graphs import generate_molecular_graph
from molmod.units import angstrom
from molmod.similarity import DistanceDescriptor

from molmod.io.xyz import XYZFile

import unittest, numpy


__all__ = ["SimilarityTestCase"]


class SimilarityTestCase(unittest.TestCase):
    def get_molecules(self):
        tpa = XYZFile("input/tpa.xyz").get_molecule()
        tpa.title = "tpa"
        tea = XYZFile("input/tea.xyz").get_molecule()
        tea.title = "tea"
        water = XYZFile("input/water.xyz").get_molecule()
        water.title = "water"
        cyclopentane = XYZFile("input/cyclopentane.xyz").get_molecule()
        cyclopentane.title = "cyclopentane"

        return [tpa, tea, water, cyclopentane]
        #return [water, cyclopentane]

    def test_mol(self):
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.descriptor = DistanceDescriptor(molecule)
        self.do_test(molecules, margin=0.2*angstrom, cutoff=7.0*angstrom)

    def test_graph(self):
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.graph = generate_molecular_graph(molecule)
            molecule.descriptor = DistanceDescriptor(molecule.graph)
        self.do_test(molecules, margin=0.2, cutoff=10.0)

    def do_test(self, molecules, margin, cutoff, verbose=False):
        if verbose:
            print
        for molecule in molecules:
            molecule.norm = molecule.descriptor.norm(margin, cutoff)
            if verbose:
                print molecule.title, "norm:", molecule.norm
        if verbose:
            print
            print " "*14, "".join("%15s" % molecule.title for molecule in molecules)
        result = []
        for index1, molecule1 in enumerate(molecules):
            if verbose: print "%15s" % molecule1.title,
            row = []
            result.append(row)
            for index2, molecule2 in enumerate(molecules):
                similarity = (
                    molecule1.descriptor.similarity(molecule2.descriptor, margin, cutoff)/
                    (molecule1.norm*molecule2.norm)
                )
                row.append(similarity)
                if verbose: print ("%14.5f" % similarity),
            if verbose: print
        result = numpy.array(result)
        self.assert_((abs(numpy.diag(result) - 1) < 1e-5).all(), "Diagonal must be unity.")
        self.assert_((abs(result - result.transpose()) < 1e-5).all(), "Result must be symmetric.")




