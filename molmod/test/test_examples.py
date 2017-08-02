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
#--


from __future__ import print_function

import os
import unittest
import shutil
import subprocess
import tempfile

import pkg_resources

from molmod import *


def check_example(dirname, fn_py, fns_data):
    dntmp = tempfile.mkdtemp(dirname, fn_py)
    try:
        for fn in [fn_py] + fns_data:
            with pkg_resources.resource_stream(__name__, "../examples/{}/{}".format(dirname, fn)) as fin:
                with open(os.path.join(dntmp, fn), 'wb') as fout:
                    fout.write(fin.read())
        env = dict(os.environ)
        root_dir = os.getcwd()
        env['PYTHONPATH'] = root_dir + ':' + env.get('PYTHONPATH', '')
        command = ["python", fn_py]
        proc = subprocess.Popen(command, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=dntmp, env=env)
        outdata, errdata = proc.communicate()
        if proc.returncode != 0:
            lines = [
                'Command faild', str(command), 'Standard output', '+'*80, outdata.decode('utf-8'),
                '+'*80, 'Standard error', '+'*80, errdata.decode('utf-8'), '+'*80]
            raise AssertionError('\n'.join(lines))
    finally:
        shutil.rmtree(dntmp)


def test_example_000_a():
    check_example("000_units", "a_reaction.py", [])

def test_example_000_b():
    check_example("000_units", "b_chbond.py", [])

def test_example_000_c():
    check_example("000_units", "c_h2rot.py", [])

def test_example_001_a():
    check_example("001_molecules", "a_convert.py", ['ibuprofen.sdf'])

def test_example_001_b():
    check_example("001_molecules", "b_com.py", ['ibuprofen.sdf'])

def test_example_001_c():
    check_example("001_molecules", "c_carbon.py", ['ibuprofen.sdf'])

def test_example_001_d():
    check_example("001_molecules", "d_size.py", ['ibuprofen.sdf'])

def test_example_001_e():
    check_example("001_molecules", "e_shape.py", ['ibuprofen.sdf'])

def test_example_002_a():
    check_example("002_graphs", "a_graphs.py", ['caffeine.xyz'])

def test_example_002_b():
    check_example("002_graphs", "b_neighbors.py", ['caffeine.xyz'])

def test_example_002_c():
    check_example("002_graphs", "c_distances.py", ['caffeine.xyz'])

def test_example_002_d():
    check_example("002_graphs", "d_symmetries.py", ['ethanol.xyz'])

def test_example_003_a():
    check_example("003_internal_coordinates", "a_bond_length.py", ['dopamine.xyz'])

def test_example_003_b():
    check_example("003_internal_coordinates", "b_bending_angles.py", ['dopamine.xyz'])

def test_example_003_c():
    check_example("003_internal_coordinates", "c_ff_hessian.py", ['propane.xyz'])

def test_example_003_d():
    check_example("003_internal_coordinates", "d_dft_hessian.py", ['dopamine.fchk'])

def test_example_004_a():
    check_example("004_patterns", "a_propane_types.py", ['propane.xyz'])

def test_example_004_b():
    check_example("004_patterns", "b_dopamine_types.py", ['dopamine.xyz'])
