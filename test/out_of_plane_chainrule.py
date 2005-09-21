# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

from pychem.internal_coordinates import Collection
from pychem.molecular_graphs import OutOfPlaneSets, CriteriaSet
from pychem.molecules import molecule_from_xyz_filename

import unittest, math, Numeric


__all__ = ["suite"]

suite = unittest.TestSuite()   

class TestOutOfPlaneChainrule(unittest.TestCase):        
    def test_dihedral_chainrule(self):
        chlorobromoethene = molecule_from_xyz_filename("input/chlorobromoethene.xyz")
        chlorobromoethene_mod = molecule_from_xyz_filename("input/chlorobromoethene_oop.xyz")

        # Define a set of independant internal coordinates.
        collection = Collection(chlorobromoethene)
        collection.add_out_of_plane_cosines(OutOfPlaneSets([CriteriaSet("CC(HCl)-angles", ((6, 6, 1, 17), None))]))
        out_of_plane = collection["CC(HCl)-angles"][0]
        value, gradient = out_of_plane(chlorobromoethene.coordinates)
        delta = chlorobromoethene.coordinates - chlorobromoethene_mod.coordinates 
        cos_estimate = Numeric.dot(Numeric.ravel(gradient), Numeric.ravel(delta))
        cos_real = 1-math.cos(10.0/180.0*math.pi)
        self.assertAlmostEqual(
            cos_real, 
            cos_estimate, 
            1, 
            "Difference in estimate and real cosine: %f-%f=%f" % (cos_estimate, cos_real, cos_estimate-cos_real)
        )

suite.addTests([
    unittest.makeSuite(TestOutOfPlaneChainrule)
])
