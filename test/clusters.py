# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.clusters import Cluster, ClusterFactory

import numpy

import unittest, random


__all__ = ["ClusterTestCase"]


class ClusterTestCase(unittest.TestCase):
    def test_blind(self):
        cf = ClusterFactory()
        for counter in xrange(1000):
            a = random.randint(0, 100)
            b = random.randint(0, 100)
            if (a+b)%2 == 0:
                cf.add_members([a, b])

        for cluster in cf.get_clusters():
            #print cluster.members
            tmp = numpy.array(cluster.members) % 2
            #print tmp
            self.assert_((tmp == 0).all() or (tmp == 1).all())

