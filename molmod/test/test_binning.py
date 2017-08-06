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
import unittest

import numpy as np
import pkg_resources
from nose.plugins.skip import SkipTest

from molmod import *
from molmod.io import *
from molmod.periodic import periodic
from molmod.test.test_unit_cells import get_random_uc


__all__ = ["BinningTestCase"]


class BinningTestCase(unittest.TestCase):
    def verify_distances(self, iter_pairs, cutoff, distances, idcls, unit_cell=None):
        # a few sanity checks first
        for key, value in distances:
            self.assertEqual(len(key), 2, "Singletons encountered: %s" % str(key))
        count = len(distances)
        distances = dict(distances)
        self.assertEqual(len(distances), count, "Duplicate distances: %i > %i" % (count, len(distances)))

        # real check of distances
        missing_pairs = []
        wrong_distances = []
        num_total = 0
        num_correct = 0

        for (id0, coord0), (id1, coord1) in iter_pairs():
            delta = coord1 - coord0
            if unit_cell is not None:
                delta = unit_cell.shortest_vector(delta)
            distance = np.linalg.norm(delta)
            if distance < cutoff:
                num_total += 1
                identifier = idcls([id0, id1])
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
        for identifier, fast_distance in distances.items():
            message += "%10s %10s: \t % 10.7f\n" % (identifier, fast_distance)
        message += "TOTAL PAIRS: %i\n" % num_total
        message += "CORRECT PAIRS: %i\n" % num_correct
        message += "-"*50+"\n"

        self.assertEqual(len(missing_pairs), 0, message)
        self.assertEqual(len(wrong_distances), 0, message)
        self.assertEqual(len(distances), 0, message)

    def verify_bins_intra_periodic(self, bins):
        raise SkipTest
        neighbor_set = set([tuple(index) for index in bins.neighbor_indexes])
        self.assertEqual(len(neighbor_set), len(bins.neighbor_indexes))

        for key0, bin0 in bins:
            encountered = set([])
            for key1, bin1 in bins.iter_surrounding(key0):
                frac_key1 = bins.integer_cell.to_fractional(key1)
                self.assert_(frac_key1.max() < 0.5, str(key1) + str(frac_key1))
                self.assert_(key1 not in encountered, str(key0) + str(key1))
                encountered.add(key1)

    def verify_bins_inter_periodic(self, bins0, bins1):
        raise SkipTest
        for bins in bins0, bins1:
            neighbor_set = set([tuple(index) for index in bins.neighbor_indexes])
            self.assertEqual(len(neighbor_set), len(bins.neighbor_indexes))

        for key0, bin0 in bins0:
            encountered = set([])
            for key1, bin1 in bins1.iter_surrounding(key0):
                frac_key1 = bins0.integer_cell.to_fractional(key1)
                self.assert_(frac_key1.max() < 0.5, str(key1) + str(frac_key1))
                self.assert_(key1 not in encountered, str(key0) + str(key1))
                encountered.add(key1)

    def verify_distances_intra(self, coordinates, cutoff, distances, unit_cell=None):
        def iter_pairs():
            for index0, coord0 in enumerate(coordinates):
                for index1, coord1 in enumerate(coordinates[:index0]):
                    yield (index0, coord0), (index1, coord1)
        self.verify_distances(iter_pairs, cutoff, distances, frozenset, unit_cell)

    def verify_distances_inter(self, coordinates0, coordinates1, cutoff, distances, unit_cell=None):
        def iter_pairs():
            for index0, coord0 in enumerate(coordinates0):
                for index1, coord1 in enumerate(coordinates1):
                    yield (index0, coord0), (index1, coord1)
        self.verify_distances(iter_pairs, cutoff, distances, tuple, unit_cell)

    def test_distances_intra_lau(self):
        coordinates = XYZFile(pkg_resources.resource_filename(__name__, "../data/test/lau.xyz")).geometries[0]
        cutoff = periodic.max_radius*2
        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in PairSearchIntra(coordinates, cutoff)
        ]
        self.verify_distances_intra(coordinates, cutoff, distances)

    def test_distances_intra_lau_periodic(self):
        coordinates = XYZFile(pkg_resources.resource_filename(__name__, "../data/test/lau.xyz")).geometries[0]
        cutoff = periodic.max_radius*2
        unit_cell = UnitCell.from_parameters3(
            np.array([14.59, 12.88, 7.61])*angstrom,
            np.array([ 90.0, 111.0, 90.0])*deg,
        )

        pair_search = PairSearchIntra(coordinates, cutoff, unit_cell)
        self.verify_bins_intra_periodic(pair_search.bins)

        distances = [
            (frozenset([i0, i1]), distance)
            for i0, i1, delta, distance
            in pair_search
        ]

        self.verify_distances_intra(coordinates, cutoff, distances, unit_cell)

    def test_distances_intra_random(self):
        for i in range(10):
            coordinates = np.random.uniform(0,5,(20,3))
            cutoff = np.random.uniform(1, 6)
            distances = [
                (frozenset([i0, i1]), distance)
                for i0, i1, delta, distance
                in PairSearchIntra(coordinates, cutoff)
            ]
            self.verify_distances_intra(coordinates, cutoff, distances)

    def test_distances_intra_random_periodic(self):
        raise SkipTest
        for i in range(10):
            coordinates = np.random.uniform(0,1,(20,3))
            unit_cell = get_random_uc(5.0, np.random.randint(0, 4), 0.5)
            coordinates = unit_cell.to_cartesian(coordinates)*3-unit_cell.matrix.sum(axis=1)
            cutoff = np.random.uniform(1, 6)

            pair_search = PairSearchIntra(coordinates, cutoff, unit_cell)
            self.verify_bins_intra_periodic(pair_search.bins)

            distances = [
                (frozenset([i0, i1]), distance)
                for i0, i1, delta, distance
                in pair_search
            ]
            self.verify_distances_intra(coordinates, cutoff, distances, unit_cell)

    def test_distances_inter_random(self):
        for i in range(10):
            coordinates0 = np.random.uniform(0,5,(20,3))
            coordinates1 = np.random.uniform(0,5,(20,3))
            cutoff = np.random.uniform(1, 6)
            distances = [
                ((i0, i1), distance)
                for i0, i1, delta, distance
                in PairSearchInter(coordinates0, coordinates1, cutoff)
            ]
            self.verify_distances_inter(coordinates0, coordinates1, cutoff, distances)

    def test_distances_inter_random_periodic(self):
        raise SkipTest
        for i in range(10):
            fractional0 = np.random.uniform(0,1,(20,3))
            fractional1 = np.random.uniform(0,1,(20,3))
            unit_cell = get_random_uc(5.0, np.random.randint(0, 4), 0.5)
            coordinates0 = unit_cell.to_cartesian(fractional0)*3-unit_cell.matrix.sum(axis=1)
            coordinates1 = unit_cell.to_cartesian(fractional1)*3-unit_cell.matrix.sum(axis=1)
            cutoff = np.random.uniform(1, 6)

            pair_search = PairSearchInter(coordinates0, coordinates1, cutoff, unit_cell)
            self.verify_bins_inter_periodic(pair_search.bins0, pair_search.bins1)

            distances = [
                ((i0, i1), distance)
                for i0, i1, delta, distance
                in pair_search
            ]
            self.verify_distances_inter(coordinates0, coordinates1, cutoff, distances, unit_cell)
