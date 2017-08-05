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


from builtins import range
import os

import pkg_resources

from molmod.test.common import *
from molmod import *

__all__ = ["RandomizeTestCase"]


# These threshold values are not good for serious applications!!!
# Just for testing the code, they are OK.
nonbond_thresholds = {
    frozenset([1,1]): 0.9*angstrom,
    frozenset([1,6]): 1.4*angstrom,
    frozenset([1,7]): 1.4*angstrom,
    frozenset([1,8]): 1.4*angstrom,
    frozenset([6,6]): 2.2*angstrom,
    frozenset([6,7]): 2.2*angstrom,
    frozenset([6,8]): 2.2*angstrom,
    frozenset([7,7]): 2.2*angstrom,
    frozenset([7,8]): 2.2*angstrom,
    frozenset([8,8]): 2.2*angstrom,
}


class RandomizeTestCase(BaseTestCase):
    def iter_test_molecules(self):
        for filename in ["tpa.xyz", "water.xyz", "thf_single.xyz"]:
            molecule = Molecule.from_file(pkg_resources.resource_filename(__name__, "../data/test/" + filename))
            molecule.filename = filename
            molecule.set_default_graph()
            yield molecule

    def test_randomize_count(self):
        all_counts = {
            "tpa.xyz": (12+4*7, 12, 13*6, 0),
            "water.xyz": (2,0,1,0), # Stretch, Torsion, Bend, DoubleStretch
            "thf_single.xyz": (8, 14, 5*4, 10),
        }
        for molecule in self.iter_test_molecules():
            manipulations = generate_manipulations(molecule)
            randomized_molecule = randomize_molecule(molecule, manipulations, nonbond_thresholds)
            with tmpdir(__name__, 'test_randomize_count') as dn:
                randomized_molecule.write_to_file(os.path.join(dn, molecule.filename))
            counts = all_counts.get(molecule.filename)
            if counts is not None:
                for cls, count in zip([RandomStretch, RandomTorsion, RandomBend, RandomDoubleStretch], counts):
                    got = sum(isinstance(mpl, cls) for mpl in manipulations)
                    self.assertEqual(got, count, "%s count problem, should be %i, got %i" % (str(cls), count, got))

    def test_random_dimer(self):
        for molecule1 in self.iter_test_molecules():
            for molecule2 in self.iter_test_molecules():
                dimer = random_dimer(molecule1, molecule2, nonbond_thresholds, 0.5*angstrom)
                self.assertEqual(dimer.coordinates.shape, (molecule1.coordinates.shape[0] + molecule2.coordinates.shape[0], 3))
                self.assertEqual(dimer.numbers.shape, (molecule1.numbers.shape[0] + molecule2.numbers.shape[0],))
                with tmpdir(__name__, 'test_random_dimer') as dn:
                    dimer.write_to_file(os.path.join(dn, "%s_%s" % (molecule1.filename, molecule2.filename)))

    def test_single_manipulation(self):
        for molecule in self.iter_test_molecules():
            manipulations = generate_manipulations(molecule)
            for i in range(100):
                randomized_molecule, mol_transformation = single_random_manipulation(molecule, manipulations, nonbond_thresholds)
                with tmpdir(__name__, 'test_single_manipulation') as dn:
                    randomized_molecule.write_to_file(os.path.join(dn, molecule.filename))
                    mol_transformation.write_to_file(os.path.join(dn, "tmp"))
                    check_transformation = MolecularDistortion.read_from_file(os.path.join(dn, "tmp"))
                self.assertEqual(mol_transformation.affected_atoms, check_transformation.affected_atoms)
                self.assertArraysAlmostEqual(mol_transformation.transformation.r, check_transformation.transformation.r, 1e-5, doabs=True)
                self.assertArraysAlmostEqual(mol_transformation.transformation.t, check_transformation.transformation.t, 1e-5, doabs=True)
