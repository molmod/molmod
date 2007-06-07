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


from molmod.binning import PositionedObject, SparseBinnedObjects
from molmod.environments import calculate_environments

from ccop.xyz import XYZFile

import unittest, random


__all__ = ["EnvironmentTestCase"]


class EnvironmentTestCase(unittest.TestCase):
    def test_blind(self):
        m = XYZFile("input/precursor.xyz").get_molecule()
        def yield_positioned_objects():
            for row, coordinate in enumerate(m.coordinates):
                yield PositionedObject(row, coordinate)
        sbo = SparseBinnedObjects(yield_positioned_objects(), 7)
        environments = calculate_environments(sbo, 7)
        #print environments
        #for id, environment in environments.iteritems():
        #    print id, environment.id
        #    print environment.deltas
        #    print environment.distances
        #    print environment.neighbors
        #    print environment.reverse_neighbors
        #    print environment.directions

