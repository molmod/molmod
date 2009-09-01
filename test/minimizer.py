# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from molmod.molecules import *

from molmod.minimizer import *

import unittest, numpy


__all__ = ["MinimizerTestCase"]


def fun(x, do_gradient=False):
    value = 2 + numpy.sin(x[0]) + numpy.cos(x[1]) + x[0]*x[0] + x[1]*x[1] - x[0]*x[1]
    if do_gradient:
        gradient = numpy.array([
            numpy.cos(x[0]) + 2*x[0] - x[1],
            -numpy.sin(x[1]) + 2*x[1] - x[0],
        ])
        return value, gradient
    else:
        return value


class MinimizerTestCase(unittest.TestCase):
    def check_min(self, x_opt, x_tol, f_tol):
        f_opt = fun(x_opt)
        for i in xrange(100):
            x_dev = x_opt + numpy.random.uniform(-x_tol, x_tol, 2)
            f_dev = fun(x_dev)
            self.assert_(f_opt - f_dev <= f_tol)

    def test_golden(self):
        x_init = numpy.zeros(2, float)
        minimizer = Minimizer(
            x_init, fun, GoldenLineSearch, 1e-5, 1e-5, 1e-1, 1000, 50,
            do_gradient=True, verbose=False,
        )
        self.check_min(minimizer.x, 1e-5, 1e-5)

    def test_newton(self):
        x_init = numpy.zeros(2, float)
        minimizer = Minimizer(
            x_init, fun, NewtonLineSearch, 1e-5, 1e-5, 1e-1, 1000, 50,
            do_gradient=True, verbose=False,
        )
        self.check_min(minimizer.x, 1e-5, 1e-5)

    def test_newtong(self):
        x_init = numpy.zeros(2, float)
        minimizer = Minimizer(
            x_init, fun, NewtonGLineSearch, 1e-5, 1e-5, 1e-1, 1000, 50,
            do_gradient=True, verbose=False,
        )
        self.check_min(minimizer.x, 1e-5, 1e-5)



