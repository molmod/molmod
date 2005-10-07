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
from pychem.interfaces.mpqc.simple import SimpleMpqcJobSinglePoint, SimpleMpqcJobOptimize
from pychem.interfaces.mpqc.oo import OOMpqcJob
from pychem.interfaces.mpqc.kvo import create_single_point_kv, create_optimize_kv
from pychem.interfaces.mpqc.keyval import KeyValObject
from pychem.interfaces.mpqc.file_parsers import MolecularEnergiesParser, OutputMoleculesParser, GradientsParser, HessianParser, OptimizationConvergedParser
from pychem.interfaces.output_parsers import OutputParser
from pychem.molecules import molecule_from_xyz_filename

import math, Numeric, LinearAlgebra
import unittest

__all__ = ["SimpleMpqcInterface", "OOMpqcInterface"]


class SimpleMpqcInterface(unittest.TestCase):
    def test_single_point(self):
        def validate():
            self.assert_(job.completed)
            self.assertAlmostEqual(job.energy, -75.9734488121, 8)
            self.assertAlmostEqual(job.gradient[0,0],  0.01174361, 6)
            self.assertAlmostEqual(job.gradient[0,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[0,2],  0.0,        6)
            self.assertAlmostEqual(job.gradient[1,0], -0.0058718,  6)
            self.assertAlmostEqual(job.gradient[1,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[1,2], -0.01381411, 6)
            self.assertAlmostEqual(job.gradient[2,0], -0.0058718,  6)
            self.assertAlmostEqual(job.gradient[2,1],  0.0,        6)
            self.assertAlmostEqual(job.gradient[2,2],  0.01381411, 6)
            
        water = molecule_from_xyz_filename("input/water.xyz")
        job = SimpleMpqcJobSinglePoint(
            prefix="output/water_sp",
            title="Water single point berekening", 
            input_molecule=water,
            method="KS (xc = B3LYP)",
            basis="3-21G*",
            memory="100MB",
            do_gradient=True
        )
        job.run(cleanup=True)
        validate()
        job.run(forcerun=True)
        validate()
        job.run()
        validate()
        
        filename = job.filename
        job = reload_job(filename + ".job")
        job.run()
        validate()
                  
    def test_optimize(self):
        def validate():
            self.assert_(job.completed)
            self.assertAlmostEqual(job.energies[-1], -75.973963163199997, 8)
            coordinates = job.output_molecule.coordinates
            delta = coordinates[0]-coordinates[1]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 6)
            delta = coordinates[0]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 6)
            delta = coordinates[1]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 2.96668446577, 6)
        
        water = molecule_from_xyz_filename("input/water.xyz")
        job = SimpleMpqcJobOptimize(
            prefix="output/water_opt",
            title="Water single point berekening", 
            input_molecule=water,
            method="KS (xc = B3LYP)",
            basis="3-21G*",
            memory="100MB"
        )
        job.run(cleanup=True)
        validate()
        job.run(forcerun=True)
        validate()
        job.run()
        validate()
        
        filename = job.filename
        job = reload_job(filename + ".job")
        job.run()
        validate()

class OOMpqcInterface(unittest.TestCase):
    def test_single_point(self):
        def validate():
            energy = job.energies[-1]
            gradient = job.gradients[-1]
            self.assert_(job.completed)
            self.assertAlmostEqual(energy, -75.9734488121, 8)
            self.assertAlmostEqual(gradient[0,0],  0.01174361, 6)
            self.assertAlmostEqual(gradient[0,1],  0.0,        6)
            self.assertAlmostEqual(gradient[0,2],  0.0,        6)
            self.assertAlmostEqual(gradient[1,0], -0.0058718,  6)
            self.assertAlmostEqual(gradient[1,1],  0.0,        6)
            self.assertAlmostEqual(gradient[1,2], -0.01381411, 6)
            self.assertAlmostEqual(gradient[2,0], -0.0058718,  6)
            self.assertAlmostEqual(gradient[2,1],  0.0,        6)
            self.assertAlmostEqual(gradient[2,2],  0.01381411, 6)
            
        water = molecule_from_xyz_filename("input/water.xyz")
        keyval = create_single_point_kv(
            molecule=water,
            charge=0,
            method="CLKS",
            basis="3-21G*",
            functional="B3LYP"
        )
        keyval['mpqc']['do_gradient'] = 'yes'
        job = OOMpqcJob(
            prefix="output/water_oo_sp",
            title="Water single point berekening", 
            keyval=keyval,
            output_parser=OutputParser([
                MolecularEnergiesParser('energies'),
                GradientsParser('gradients')
            ])
        )
        job.run(cleanup=True)
        validate()
        job.run(forcerun=True)
        validate()
        job.run()
        validate()
        
        filename = job.filename
        job = reload_job(filename + ".job")
        job.run()
        validate()

    def test_optimize(self):
        def validate():
            self.assert_(job.completed)
            self.assertAlmostEqual(job.energies[-1], -75.973963163199997, 8)
            coordinates = job.output_molecules[-1].coordinates
            delta = coordinates[0]-coordinates[1]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 3)
            delta = coordinates[0]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 1.88335259871, 3)
            delta = coordinates[1]-coordinates[2]
            self.assertAlmostEqual(math.sqrt(Numeric.dot(delta, delta)), 2.96668446577, 3)
            evals = LinearAlgebra.eigenvalues(job.hessian)
            self.assert_(min(evals) > -1e-6, "All eigenvalues of a hessian of an optimized molecule should be positive.\n%s" % str(evals))
        
        water = molecule_from_xyz_filename("input/water.xyz")
        keyval = create_optimize_kv(
            molecule=water,
            charge=0,
            method="CLKS",
            basis="3-21G*",
            functional="B3LYP"
        )
        keyval['mpqc']['freq'] = KeyValObject('MolecularFrequencies', named_items=[
            ('molecule', keyval['molecule'])
        ])
        job = OOMpqcJob(
            prefix="output/water_oo_opt",
            title="Water single point berekening", 
            keyval=keyval,
            output_parser=OutputParser([
                MolecularEnergiesParser('energies'),
                OutputMoleculesParser('output_molecules'),
                OptimizationConvergedParser('converged'),
                HessianParser('hessian')
            ])
        )
        job.run(cleanup=True)
        validate()
        if not job.converged:
            # MPQC has a problem when a converged optimization is restarted
            job.run(forcerun=True)
            validate()
        job.run()
        validate()
        
        filename = job.filename
        job = reload_job(filename + ".job")
        job.run()
        validate()

