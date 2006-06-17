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

from molmod.interfaces.base import reload_job
from molmod.interfaces.gromos96.sp import Gromos96SP
from molmod.interfaces.gromos96.file_parsers import *
from molmod.interfaces.output_parsers import OutputParser
from molmod.molecules import molecule_xyz_from_filename
from molmod.units import from_angstrom

import math, numpy, numpy.linalg, os
import unittest

__all__ = ["Gromos96Interface"]


class Gromos96Interface(unittest.TestCase):
    def test_single_point(self):
        def validate():
            self.assert_(job.completed)
            #print job.forces1
            #print job.forces2
            #print job.potential_energy
        
        parameters = {
            "bs_coef_51": 1.0700e+6,
            "bs_rest_51": 0.163980,
            "bs_coef_52": 1.0400e+6,
            "bs_rest_52": 0.164470,
            "bs_coef_53": 1.430e+7,
            "bs_rest_53": 0.09589,
            "ba_coef_51": 458,
            "ba_rest_51": 136.30,
            "ba_coef_52": 2206,
            "ba_rest_52": 108.09,
            "ba_coef_53": 268,
            "ba_rest_53": 116.89,
            "ba_coef_54": 897,
            "ba_rest_54": 109.06,
            "ba_coef_55": 935,
            "ba_rest_55": 106.36
        }
        
        cluster = molecule_xyz_from_filename("input/2TOH.xyz")
        os.system("rm -rf output/gromos96")
        #os.mkdir("output/gromos96")
        job = Gromos96SP(
            prefix="output/gromos96/",
            title="Zeolite cluster 2TOH", 
            input_molecule=cluster, 
            parameters=parameters, 
            topology="2TOH", 
            box_size=from_angstrom(15.0),
            output_parser=OutputParser([
                Forces1Parser('forces1'),
                Forces2Parser('forces2'),
                MDOutDataParser(
                    section=' 6. D A T A   P E R   S T E P',
                    subsection='       ENERGY        TOTAL      KINETIC    POTENTIAL',
                    row='E-TOT',
                    col=3,
                    conversion=1,
                    label='potential_energy'
                )
            ])
        )
        job.run(cleanup=True)
        validate()
        job.run()
        validate()        
