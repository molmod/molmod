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

import unittest

import pkg_resources

from molmod.test.common import *
from molmod.io import *
from molmod import *


__all__ = ["CP2KTestCase"]


class CP2KTestCase(unittest.TestCase):
    def test_input_file(self):
        cp2k_input = CP2KInputFile.read_from_file(pkg_resources.resource_filename(__name__, "../../data/test/water_md.inp"))
        self.assert_(cp2k_input._consistent())

        # test that the read_from_file works
        self.assertEqual(cp2k_input["FORCE_EVAL"]["METHOD"].value, "Quickstep")
        self.assert_(cp2k_input._consistent())
        self.assertEqual(cp2k_input["MOTION"]["MD"]["ENSEMBLE"].value, "NVE")
        self.assert_(cp2k_input._consistent())

        # test the __getitem__, __setitem__ and __delitem__
        del cp2k_input["MOTION"]["MD"]["ENSEMBLE"]
        self.assert_(cp2k_input._consistent())
        try:
            print(cp2k_input["MOTION"]["MD"]["ENSEMBLE"].value)
            self.fail("CP2KKeyword ENSEMBLE should no longer exist.")
        except KeyError:
            pass

        cp2k_input["MOTION"]["MD"]["ENSEMBLE"] = CP2KKeyword("ENSEMBLE", "NVE")
        self.assert_(cp2k_input._consistent())
        self.assertEqual(cp2k_input["MOTION"]["MD"]["ENSEMBLE"].value, "NVE")

        try:
            cp2k_input["MOTION"]["MD"]["ENSEMBLE"] = CP2KKeyword("JOS", "NVE")
            self.fail("CP2KKeyword should have the correct name.")
        except KeyError:
            pass

        try:
            cp2k_input["MOTION"]["MD"]["ENSEMBLE"] = [CP2KKeyword("ENSEMBLE", "NVE"), CP2KKeyword("JOS", "NVE")]
            self.fail("CP2KKeyword should have the correct name.")
        except KeyError:
            pass

        cp2k_input["MOTION"]["MD"]["ENSEMBLE"] = [CP2KKeyword("ENSEMBLE", "NVE"), CP2KKeyword("ENSEMBLE", "NVT")]
        self.assert_(cp2k_input._consistent())
        cp2k_input["MOTION"]["MD"]["ENSEMBLE", 0] = CP2KKeyword("ENSEMBLE", "NVE")
        self.assert_(cp2k_input._consistent())

        self.assertEqual(len(cp2k_input["MOTION"]["MD"]["ENSEMBLE"]), 2)
        l = cp2k_input["MOTION"]["MD"]["ENSEMBLE"]
        self.assertEqual(l[0].value, "NVE")
        self.assertEqual(l[1].value, "NVT")

        del cp2k_input["MOTION"]["MD"]["ENSEMBLE", 0]
        self.assert_(cp2k_input._consistent())
        self.assertEqual(cp2k_input["MOTION"]["MD"]["ENSEMBLE"].value, "NVT")

        # test __len__
        self.assertEqual(len(cp2k_input), 3)
        self.assertEqual(len(cp2k_input["MOTION"]["MD"]), 5)
        cp2k_input["MOTION"]["MD"]["ENSEMBLE"] = [CP2KKeyword("ENSEMBLE", "NVE"), CP2KKeyword("ENSEMBLE", "NVT")]
        self.assert_(cp2k_input._consistent())
        self.assertEqual(len(cp2k_input["MOTION"]["MD"]), 6)
        del cp2k_input["MOTION"]["MD"]["ENSEMBLE", 1]
        self.assert_(cp2k_input._consistent())

        # test creating new parts
        nose = CP2KSection("NOSE", [
            CP2KKeyword("LENGTH", "3"),
            CP2KKeyword("TIMECON", "10.0")
        ])
        self.assert_(nose._consistent())

        # test iter
        tmp = [child for child in nose]
        self.assertEqual(len(tmp), 2)
        del tmp

        # test append
        cp2k_input["MOTION"]["MD"].append(nose)
        self.assert_(cp2k_input._consistent())
        self.assertEqual(cp2k_input["MOTION"]["MD"]["NOSE"]["LENGTH"].value, "3")

        # test dump, load consistency, part 1: dump a file, load it again, should be the same
        with tmpdir(__name__, 'test_input_file') as dn:
            cp2k_input.write_to_file("%s/water_md.inp" % dn)
            cp2k_input_check = CP2KInputFile.read_from_file("%s/water_md.inp" % dn)
        self.assert_(cp2k_input_check._consistent())
        self.assertEqual(cp2k_input, cp2k_input_check)

        # test dump-load consistency, part 2: no reordering of sections and keywords should be allowed
        cp2k_input = CP2KInputFile.read_from_file(pkg_resources.resource_filename(__name__, "../../data/test/water_md.inp"))
        with tmpdir(__name__, 'test_input_file') as dn:
            cp2k_input.write_to_file("%s/water_md.inp" % dn)
            with open(pkg_resources.resource_filename(__name__, "../../data/test/water_md.inp")) as f1, \
                 open("%s/water_md.inp" % dn) as f2:
                for line1, line2 in zip(f1, f2):
                    self.assertEqual(line1, line2)
