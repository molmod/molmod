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
    def verify_distances(self, coordinates, cutoff, distances, unit_cell=None):
        # a few sanity checks first
        for key, value in distances:
            self.assertEqual(len(key), 2, "Singletons encountered: %s" % (key))
        count = len(distances)
        distances = dict(distances)
        self.assertEqual(len(distances), count, "Duplicate distances: %i > %i" % (count, len(distances)))

        # real check of distances
        missing_pairs = []
        wrong_distances = []
        num_total = 0
        num_correct = 0

        def iter_pairs():
            for index1, coord1 in enumerate(coordinates):
                for index2, coord2 in enumerate(coordinates[:index1]):
                    yield (index1, coord1), (index2, coord2)

        for (id1, coord1), (id2, coord2) in iter_pairs():
            delta = coord2 - coord1
            if unit_cell is not None:
                delta = unit_cell.shortest_vector(delta)
            distance = numpy.linalg.norm(delta)
            if distance < cutoff:
                num_total += 1
                identifier = frozenset([id1, id2])
                fast_distance = distances.get(identifier)
                if fast_distance is None:
                    missing_pairs.append(tuple(identifier) + (distance,))
                elif fast_distance != distance:
                    wrong_distances.append(tuple(identifier) + (fast_distance, distance))
                else:
                    num_correct += 1
                    del distances[identifier]

        message  = "-"*50+"\n"
        message += "CUTOFF %s\n" % cutoff
        message += "MISSING PAIRS: %i\n" % len(missing_pairs)
        for missing_pair in missing_pairs:
            message += "%10s %10s: \t % 10.7f\n" % missing_pair
        message += "WRONG DISTANCES: %i\n" % len(wrong_distances)
        for wrong_distance in wrong_distances:
            message += "%10s %10s: \t % 10.7f != % 10.7f\n" % wrong_distance
        message += "UNWANTED PAIRS: %i\n" % len(distances)
        for identifier, fast_distance in distances.iteritems():
            message += "%10s %10s: \t % 10.7f\n" % (identifier, fast_distance)
        message += "TOTAL PAIRS: %i\n" % num_total
        message += "CORRECT PAIRS: %i\n" % num_correct
        message += "-"*50+"\n"

        self.assertEqual(len(missing_pairs), 0, message)
        self.assertEqual(len(wrong_distances), 0, message)
        self.assertEqual(len(distances), 0, message)

    def verify_internals_periodic(self, pair_search):
        neighbor_set = set([tuple(index) for index in pair_search.neighbor_indexes])
        self.assertEqual(len(neighbor_set), len(pair_search.neighbor_indexes))

        for key0, bin0 in pair_search.bins.iteritems():
            encountered = set([])
            for key1, bin1 in pair_search.iter_surrounding(key0):
                frac_key1 = pair_search.integer_cell.to_fractional(key1)
                self.assert_(frac_key1.max() < 0.5, str(key1) + str(frac_key1))
                self.assert_(key1 not in encountered, str(key0) + str(key1))
                encountered.add(key1)


    def test_distances_lau(self):
        coordinates = XYZFile("input/lau.xyz").geometries[0]
        cutoff = periodic.max_radius*2
        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in PairSearch(coordinates, cutoff)
        ]
        self.verify_distances(coordinates, cutoff, distances)

    def test_distances_lau_periodic(self):
        coordinates = XYZFile("input/lau.xyz").geometries[0]
        cutoff = periodic.max_radius*2
        unit_cell = UnitCell.from_parameters3(
            numpy.array([14.59, 12.88, 7.61])*angstrom,
            numpy.array([ 90.0, 111.0, 90.0])*deg,
        )

        pair_search = PairSearch(coordinates, cutoff, unit_cell)
        self.verify_internals_periodic(pair_search)

        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in pair_search
        ]

        self.verify_distances(coordinates, cutoff, distances, unit_cell)

    def test_distances_random(self):
        for i in xrange(10):
            coordinates = numpy.random.uniform(0,5,(100,3))
            cutoff = 1.0
            distances = [
                (frozenset([i0, i1]), distance)
                for i0, i1, delta, distance
                in PairSearch(coordinates, cutoff)
            ]
            self.verify_distances(coordinates, cutoff, distances)

    def test_distances_random_periodic(self):
        for i in xrange(10):
            coordinates = numpy.random.uniform(0,1,(10,3))
            while True:
                unit_cell = UnitCell(
                    numpy.random.uniform(0,5,(3,3)),
                    numpy.random.randint(0,2,3).astype(bool),
                )
                if unit_cell.spacings.min() > 0.5:
                    break
            coordinates = unit_cell.to_cartesian(coordinates)*3-unit_cell.matrix.sum(axis=1)
            cutoff = numpy.random.uniform(0, 6)

            pair_search = PairSearch(coordinates, cutoff, unit_cell)
            self.verify_internals_periodic(pair_search)

            distances = [
                (frozenset([i0, i1]), distance)
                for i0, i1, delta, distance
                in pair_search
            ]
            self.verify_distances(coordinates, cutoff, distances, unit_cell)


