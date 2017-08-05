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


from __future__ import print_function

import os
import unittest
import shutil
import stat
import subprocess
import tempfile

import pkg_resources

from molmod import *
from molmod.test.common import tmpdir


def check_example(module_name, dirname, fn_script, fns_data):
    """Run an example in a temporary directory and check its exit code.

    Parameters
    ----------
    module_name : str
        You can just pass __name__.
    dirname : str
        The directory with the example, relative to the __file__ of where you call this
        function.
    fn_script : str
        The name of the script to be executed, assumed to be present in the given
        directory.
    fns_data : list of str:
        A list of data files needed by the example, which will be copied over to the
        temporary directory.
    """
    with tmpdir(module_name, dirname + fn_script) as dntmp:
        for fn in [fn_script] + fns_data:
            with pkg_resources.resource_stream(module_name, "../examples/{}/{}".format(dirname, fn)) as fin:
                # Create the directory if needed.
                if '/' in fn:
                    subdntmp = os.path.join(dntmp, os.path.dirname(fn))
                    if not os.path.isdir(subdntmp):
                        os.makedirs(subdntmp)
                # Extract the file manually.
                with open(os.path.join(dntmp, fn), 'wb') as fout:
                    fout.write(fin.read())
        env = dict(os.environ)
        root_dir = os.getcwd()
        env['PYTHONPATH'] = root_dir + ':' + env.get('PYTHONPATH', '')
        path_script = os.path.join(dntmp, fn_script)
        os.chmod(path_script, os.stat(path_script).st_mode | stat.S_IXUSR)
        command = ["python", fn_script]
        proc = subprocess.Popen(command, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=dntmp, env=env)
        outdata, errdata = proc.communicate()
        if proc.returncode != 0:
            lines = [
                'Command faild', str(command), 'Standard output', '+'*80, outdata.decode('utf-8'),
                '+'*80, 'Standard error', '+'*80, errdata.decode('utf-8'), '+'*80]
            raise AssertionError('\n'.join(lines))


def test_example_000_a():
    check_example(__name__, "000_units", "a_reaction.py", [])

def test_example_000_b():
    check_example(__name__, "000_units", "b_chbond.py", [])

def test_example_000_c():
    check_example(__name__, "000_units", "c_h2rot.py", [])

def test_example_001_a():
    check_example(__name__, "001_molecules", "a_convert.py", ['ibuprofen.sdf'])

def test_example_001_b():
    check_example(__name__, "001_molecules", "b_com.py", ['ibuprofen.sdf'])

def test_example_001_c():
    check_example(__name__, "001_molecules", "c_carbon.py", ['ibuprofen.sdf'])

def test_example_001_d():
    check_example(__name__, "001_molecules", "d_size.py", ['ibuprofen.sdf'])

def test_example_001_e():
    check_example(__name__, "001_molecules", "e_shape.py", ['ibuprofen.sdf'])

def test_example_002_a():
    check_example(__name__, "002_graphs", "a_graphs.py", ['caffeine.xyz'])

def test_example_002_b():
    check_example(__name__, "002_graphs", "b_neighbors.py", ['caffeine.xyz'])

def test_example_002_c():
    check_example(__name__, "002_graphs", "c_distances.py", ['caffeine.xyz'])

def test_example_002_d():
    check_example(__name__, "002_graphs", "d_symmetries.py", ['ethanol.xyz'])

def test_example_003_a():
    check_example(__name__, "003_internal_coordinates", "a_bond_length.py", ['dopamine.xyz'])

def test_example_003_b():
    check_example(__name__, "003_internal_coordinates", "b_bending_angles.py", ['dopamine.xyz'])

def test_example_003_c():
    check_example(__name__, "003_internal_coordinates", "c_ff_hessian.py", ['propane.xyz'])

def test_example_003_d():
    check_example(__name__, "003_internal_coordinates", "d_dft_hessian.py", ['dopamine.fchk'])

def test_example_004_a():
    check_example(__name__, "004_patterns", "a_propane_types.py", ['propane.xyz'])

def test_example_004_b():
    check_example(__name__, "004_patterns", "b_dopamine_types.py", ['dopamine.xyz'])
