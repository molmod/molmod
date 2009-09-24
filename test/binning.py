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


from molmod.binning import *
from molmod.unit_cells import UnitCell
from molmod.units import angstrom, deg
from molmod.periodic import periodic

from molmod.io.xyz import XYZFile

import numpy, unittest


__all__ = ["BinningTestCase"]


class BinningTestCase(unittest.TestCase):
    cutoff = periodic.max_radius*2

    def verify(self, molecule, distances, unit_cell=None):
        # a few sanity checks first
        for key, value in distances:
            self.assertEqual(len(key), 2, "Singletons encountered: %s" % (key))
        count = len(distances)
        distances = dict(distances)
        self.assertEqual(len(distances), count, "Duplicate distances: %i > %i" % (count, len(distances)))

        # real check of distances
        missing_pairs = []
        wrong_distances = []
        total = 0

        def iter_pairs():
            for index1, coord1 in enumerate(molecule.coordinates):
                for index2, coord2 in enumerate(molecule.coordinates[:index1]):
                    yield (index1, coord1), (index2, coord2)

        for (id1, coord1), (id2, coord2) in iter_pairs():
            delta = coord2 - coord1
            if unit_cell is not None:
                delta = unit_cell.shortest_vector(delta)
            distance = numpy.linalg.norm(delta)
            if distance < self.cutoff:
                total += 1
                identifier = frozenset([id1, id2])
                fast_distance = distances.get(identifier)
                if fast_distance is None:
                    missing_pairs.append(tuple(identifier) + (distance,))
                elif fast_distance != distance:
                    wrong_distances.append(tuple(identifier) + (fast_distance, distance))
                else:
                    del distances[identifier]

        message  = "-"*50+"\n"
        message += "CUTOFF %s\n" % self.cutoff
        message += "MISSING PAIRS: %i\n" % len(missing_pairs)
        for missing_pair in missing_pairs:
            message += "%10s %10s: \t % 10.7f\n" % missing_pair
        message += "WRONG DISTANCES: %i\n" % len(wrong_distances)
        for wrong_distance in wrong_distances:
            message += "%10s %10s: \t % 10.7f != % 10.7f\n" % wrong_distance
        message += "DUPLICATE PAIRS: %i\n" % len(distances)
        for identifier, fast_distance in distances.iteritems():
            message += "%10s %10s: \t % 10.7f\n" % (identifier, fast_distance)
        message += "TOTAL PAIRS: %i\n" % total
        message += "-"*50+"\n"

        self.assertEqual(len(missing_pairs), 0, message)
        self.assertEqual(len(wrong_distances), 0, message)
        self.assertEqual(len(distances), 0, message)

    def test_distances(self):
        molecule = XYZFile("input/lau.xyz").get_molecule()
        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in PairSearch(molecule.coordinates, self.cutoff)
        ]
        self.verify(molecule, distances)

    def test_distances_periodic(self):
        molecule = XYZFile("input/lau.xyz").get_molecule()
        unit_cell = UnitCell.from_parameters3(
            numpy.array([14.59, 12.88, 7.61])*angstrom,
            numpy.array([ 90.0, 111.0, 90.0])*deg,
        )
        pair_search = PairSearch(molecule.coordinates, self.cutoff, unit_cell)

        for key0, bin0 in pair_search.bins.iteritems():
            encountered = set([])
            for key1, bin1 in pair_search.iter_surrounding(key0):
                self.assert_(key1 not in encountered)
                encountered.add(key1)

        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in pair_search
        ]

        self.verify(molecule, distances, unit_cell)


