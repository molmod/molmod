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


from molmod.binning import PositionedObject, SparseBinnedObjects
from molmod.environments import calculate_environments

from molmod.io.xyz import XYZFile

import unittest


__all__ = ["EnvironmentTestCase"]


class EnvironmentTestCase(unittest.TestCase):
    def test_blind(self):
        m = XYZFile("input/precursor.xyz").get_molecule()
        def yield_positioned_objects():
            for row, coordinate in enumerate(m.coordinates):
                yield PositionedObject(row, coordinate)
        sbo = SparseBinnedObjects(yield_positioned_objects(), 7.1)
        environments = calculate_environments(sbo, 7)
        #print environments
        #for id, environment in environments.iteritems():
        #    print id, environment.id
        #    print environment.deltas
        #    print environment.distances
        #    print environment.neighbors
        #    print environment.reverse_neighbors
        #    print environment.directions





