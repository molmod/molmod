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


from __future__ import division

import unittest
import tempfile
import shutil
from contextlib import contextmanager


__all__ = ["BaseTestCase", "tmpdir"]


class BaseTestCase(unittest.TestCase):
    def assertArraysEqual(self, a, b):
        self.assertEqual(a.shape, b.shape, "The array shapes do not match.")
        self.assert_((a==b).all(), "The array values do not match.")

    def assertArrayConstant(self, arr, const):
        self.assert_((arr==const).all(), "Some/All array values do not match the constant.")

    def assertArraysAlmostEqual(self, a, b, threshold=1e-6, mean=False, doabs=False):
        self.assertEqual(a.shape, b.shape, "The array shapes do not match.")
        if a.shape == (0,):
            return
        if mean:
            abserr = abs(a-b).mean()
            oom = 0.5*(abs(a).mean()+abs(b).mean())
        else:
            abserr = abs(a-b).max()
            oom = 0.5*(abs(a).max()+abs(b).max())
        if doabs:
            self.assert_(abserr <= threshold, "The absolute error is larger than the given threshold: %5.3e > %5.3e" % (abserr, threshold))
        else:
            relerr = abserr/oom
            self.assert_(relerr <= threshold, "The relative error is larger than the given threshold: %5.3e > %5.3e" % (relerr, threshold))

    def assertArrayAlmostConstant(self, arr, const, relerr_threshold):
        error = abs(arr-const).max()
        oom = const
        relerr = error/oom
        self.assert_(relerr <= relerr_threshold, "The relative error is larger than the given threshold: %5.3e > %5.3e" % (relerr, relerr_threshold))

    def assertArrayAlmostZero(self, arr, abserr_threshold):
        abserr = abs(arr).max()
        self.assert_(abserr <= abserr_threshold, "The absolute error is larger than the given threshold: %5.3e > %5.3e" % (abserr, abserr_threshold))


@contextmanager
def tmpdir(suffix, prefix):
    dn = tempfile.mkdtemp(suffix, prefix)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)
