# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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


from molmod.utils import *

import unittest, pickle


__all__ = ["UtilsTestCase"]


class Test(ReadOnly):
    def __init__(self, a, b=None):
        ReadOnly.__init__(self)
        mandatory = {"a": a}
        optional = {"b": b}
        self._init_attributes(mandatory, optional)


class UtilsTestCase(unittest.TestCase):
    def test_pickle_read_only1(self):
        test1 = Test(5)
        s = pickle.dumps(test1)
        test2 = pickle.loads(s)
        self.assertEqual(test1.a, test2.a)

    def test_pickle_read_only2(self):
        test1 = Test(5, 3)
        s = pickle.dumps(test1)
        test2 = pickle.loads(s)
        self.assertEqual(test1.a, test2.a)
        self.assertEqual(test1.b, test2.b)


