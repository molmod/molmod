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


from __future__ import print_function, division

import unittest

import numpy as np
import pkg_resources

from molmod import *



__all__ = ["SimilarityTestCase"]


class SimilarityTestCase(unittest.TestCase):
    def get_molecules(self):
        tpa = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/tpa.xyz"))
        tea = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/tea.xyz"))
        water = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/water.xyz"))
        cyclopentane = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/cyclopentane.xyz"))

        return [tpa, tea, water, cyclopentane]
        #return [water, cyclopentane]

    def test_mol(self):
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.descriptor = SimilarityDescriptor.from_molecule(molecule)
        self.check(molecules, margin=0.2*angstrom, cutoff=7.0*angstrom)

    def test_graph(self):
        molecules = self.get_molecules()
        for molecule in molecules:
            molecule.set_default_graph()
            molecule.descriptor = SimilarityDescriptor.from_molecular_graph(molecule.graph)
        self.check(molecules, margin=0.2, cutoff=10.0)

    def check(self, molecules, margin, cutoff, verbose=False):
        if verbose:
            print()
        for molecule in molecules:
            molecule.norm = compute_similarity(
                molecule.descriptor,
                molecule.descriptor,
                margin,
                cutoff
             )**0.5
            if verbose:
                print(molecule.title, "norm:", molecule.norm)
        if verbose:
            print()
            print(" "*14, "".join("%15s" % molecule.title for molecule in molecules))
        result = []
        for index1, molecule1 in enumerate(molecules):
            if verbose: print("%15s" % molecule1.title, end=' ')
            row = []
            result.append(row)
            for index2, molecule2 in enumerate(molecules):
                similarity = compute_similarity(
                    molecule1.descriptor,
                    molecule2.descriptor,
                    margin, cutoff
                )/(molecule1.norm*molecule2.norm)
                row.append(similarity)
                if verbose: print(("%14.5f" % similarity), end=' ')
            if verbose: print()
        result = np.array(result)
        self.assert_((abs(np.diag(result) - 1) < 1e-5).all(), "Diagonal must be unity.")
        self.assert_((abs(result - result.transpose()) < 1e-5).all(), "Result must be symmetric.")
