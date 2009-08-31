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


from molmod.clusters import ClusterFactory

import numpy, unittest


__all__ = ["ClusterTestCase"]


class ClusterTestCase(unittest.TestCase):
    def test_even_odd(self):
        cf = ClusterFactory()
        for counter in xrange(10000):
            a = numpy.random.randint(0, 200)
            b = numpy.random.randint(0, 200)
            if (a+b)%2 == 0:
                cf.add_related(a, b)

        counter = 0
        clusters = cf.get_clusters()
        complete = set([])
        for cluster in clusters:
            tmp = numpy.array(cluster) % 2
            self.assert_((tmp == 0).all() or (tmp == 1).all())
            counter += len(cluster)
            complete |= set(cluster)
        self.assertEqual(counter, len(complete))





