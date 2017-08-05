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

from molmod.clusters import *



__all__ = ["ClusterTestCase"]


class ClusterTestCase(unittest.TestCase):
    def test_even_odd(self):
        cf = ClusterFactory()
        for counter in range(10000):
            a = np.random.randint(0, 200)
            b = np.random.randint(0, 200)
            if (a+b)%2 == 0:
                cf.add_related(a, b)

        counter = 0
        clusters = cf.get_clusters()
        complete = set([])
        for cluster in clusters:
            tmp = np.array(list(cluster.items)) % 2
            self.assert_((tmp == 0).all() or (tmp == 1).all())
            counter += len(cluster.items)
            complete |= cluster.items
        self.assertEqual(counter, len(complete))

    def test_rule_cluster(self):
        cf = ClusterFactory(RuleCluster)
        cf.add_related(RuleCluster(["x", "y"], ["x+y=1"]))
        cf.add_related(RuleCluster(["x", "z"], ["x*z=2"]))
        cf.add_related(RuleCluster(["u", "v"], ["u=v"]))
        clusters = list(cf.get_clusters())
        self.assertEqual(len(clusters), 2)
        clusters.sort(key=lambda x: len(x.items))
        for cluster in clusters:
            cluster.rules.sort()
        self.assertEqual(clusters[0].items, set(["u", "v"]))
        self.assertEqual(clusters[0].rules, ["u=v"])
        self.assertEqual(clusters[1].items, set(["x", "y", "z"]))
        self.assertEqual(clusters[1].rules, ["x*z=2", "x+y=1"])
