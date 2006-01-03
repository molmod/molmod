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

from pychem.interfaces.base import reload_job
from pychem.interfaces.gromos96.sp import Gromos96SP
from pychem.interfaces.gromos96.file_parsers import *
from pychem.interfaces.output_parsers import OutputParser
from pychem.molecules import molecule_from_xyz_filename
from pychem.units import from_angstrom

import math, Numeric, LinearAlgebra, os
import unittest

__all__ = ["Gromos96Interface"]


class Gromos96Interface(unittest.TestCase):
    def test_single_point(self):
        def validate():
            self.assert_(job.completed)
            print job.forces1
            print job.forces2
        
        cluster = molecule_from_xyz_filename("input/2TOH.xyz")
        os.system("rm -rf output/gromos96")
        #os.mkdir("output/gromos96")
        job = Gromos96SP(
            prefix="output/gromos96/",
            title="Zeolite cluster 2TOH", 
            input_molecule=cluster, 
            parameters={}, 
            topology="2TOH", 
            box_size=from_angstrom(15.0),
            output_parser=OutputParser([
                Forces1Parser('forces1'),
                Forces2Parser('forces2')
            ])
        )
        job.run(cleanup=True)
        validate()
        job.run()
        validate()        
