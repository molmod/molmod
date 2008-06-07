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

from molmod.randomize import *
from molmod.molecular_graphs import generate_molecular_graph

from molmod.io.xyz import XYZFile

import unittest, numpy, os

__all__ = ["RandomizeTestCase"]


class RandomizeTestCase(unittest.TestCase):
    def yield_test_molecules(self):
        for filename in ["tpa.xyz", "water.xyz", "thf_single.xyz"]:
            molecule = XYZFile(os.path.join("input", filename)).get_molecule()
            molecule.filename = filename
            graph = generate_molecular_graph(molecule)
            yield molecule, graph

    def test_randomize_count(self):
        all_counts = {
            "tpa.xyz": (12+4*7, 12, 13*6, 0),
            "water.xyz": (2,0,1,0), # Stretch, Torsion, Bend, DoubleStretch
            "thf_single.xyz": (8, 10, 5*4, 5),
        }
        for molecule, graph in self.yield_test_molecules():
            manipulations = generate_manipulations(graph, molecule)
            randomized_molecule = randomize_molecule(molecule, graph, manipulations)
            randomized_molecule.write_to_file(os.path.join("output", molecule.filename))
            counts = all_counts.get(molecule.filename)
            if counts is not None:
                for cls, count in zip([RandomStretch, RandomTorsion, RandomBend, RandomDoubleStretch], counts):
                    got = sum(isinstance(mpl, cls) for mpl in manipulations)
                    self.assertEqual(got, count, "%s count problem, should be %i, got %i" % (str(cls), count, got))

