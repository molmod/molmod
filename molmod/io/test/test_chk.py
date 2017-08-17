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


import os
import numpy as np

from molmod.test.common import tmpdir
from molmod.io import load_chk, dump_chk


def check_data_array(testname, data0, dtype):
    with tmpdir(__name__, testname) as dn:
        fn_test = os.path.join(dn, 'test.chk')
        dump_chk(fn_test, data0)
        data1 = load_chk(fn_test)
        assert data0.keys() == data1.keys()
        assert (data0['values'] == data1['values']).all()
        assert data1['values'].dtype == dtype


def test_strings_array():
    check_data_array('test_strings_array', {'values': ['foo', 'bar']}, np.dtype('U22'))


def test_floats_array():
    check_data_array('test_floats_array', {'values': [1.2, 3.0, 4]}, float)


def test_ints_array():
    check_data_array('test_ints_array', {'values': [1, 3, 4]}, int)


def test_bool_array():
    check_data_array('test_bool_array', {'values': [True, False, False]}, bool)


def check_data(testname, data0, dtype):
    with tmpdir(__name__, testname) as dn:
        fn_test = os.path.join(dn, 'test.chk')
        dump_chk(fn_test, data0)
        data1 = load_chk(fn_test)
        assert data0.keys() == data1.keys()
        assert data0['values'] == data1['values']
        assert isinstance(data1['values'], dtype)


def test_strings():
    check_data('test_strings_array', {'values': 'foo'}, (str, np.unicode))


def test_floats():
    check_data('test_floats_array', {'values': 1.2}, float)


def test_ints():
    check_data('test_ints_array', {'values': 42}, int)


def test_bool():
    check_data('test_bool_array', {'values': True}, bool)


def test_none():
    check_data('test_bool_array', {'values': None}, type(None))
