# -*- coding: utf-8 -*-
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


from common import BaseTestCase

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


class MinimizerTestCase(BaseTestCase):
    def check_min(self, x_opt, step_rms, grad_rms):
        # check if it is really the minimum by computing small displacements
        f_opt = fun(x_opt)
        for i in xrange(100):
            delta = numpy.random.normal(0, 1, 2)
            delta /= numpy.linalg.norm(delta)
            delta *= numpy.sqrt(len(delta))
            delta *= step_rms
            x_dev = x_opt + delta
            f_dev = fun(x_dev)
            self.assert_(f_opt - f_dev <= grad_rms*step_rms)
        # check if it is really the minimum by computing the eigen values of the
        # hessian
        hessian1 = compute_fd_hessian(fun, x_opt, 1e-4, anagrad=True)
        self.assert_((numpy.linalg.eigvalsh(hessian1) > 0).all())
        # check if it is really the minimum by computing the eigen values of the
        # hessian
        hessian2 = compute_fd_hessian(fun, x_opt, 1e-4, anagrad=False)
        self.assert_((numpy.linalg.eigvalsh(hessian2) > 0).all())
        # also compare the two hessians
        self.assert_(abs(hessian1 - hessian2).max() < 1e-4)

    def test_sd_golden(self):
        x_init = numpy.zeros(2, float)
        search_direction = SteepestDescent()
        line_search = GoldenLineSearch(qtol=1e-10, qmax=1.0, max_iter=500)
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3, grad_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=False, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_cg_golden(self):
        x_init = numpy.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = GoldenLineSearch(qtol=1e-10, qmax=1.0, max_iter=500)
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3, grad_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=False, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_sd_newton(self):
        x_init = numpy.zeros(2, float)
        search_direction = SteepestDescent()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=False, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_cg_newton(self):
        x_init = numpy.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=False, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_sd_newtong(self):
        x_init = numpy.zeros(2, float)
        search_direction = SteepestDescent()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_cg_newtong(self):
        x_init = numpy.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_cg_newtong_diag_prec(self):
        x_init = numpy.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        prec_fun = DiagonalPreconditioner(fun, 3, 1e-2)
        minimizer = Minimizer(
            x_init, prec_fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.assert_(prec_fun.scales is not None)
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_cg_newtong_full_prec(self):
        x_init = numpy.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        prec_fun = FullPreconditioner(fun, 3, 1e-2)
        minimizer = Minimizer(
            x_init, prec_fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.assert_(prec_fun.scales is not None)
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_qn_newtong(self):
        x_init = numpy.zeros(2, float)
        search_direction = QuasiNewton()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_check_anagrad(self):
        x_init = numpy.zeros(2, float)
        check_anagrad(fun, x_init, 1e-5)

    def test_check_anagrad_diag_prec(self):
        x_init = numpy.zeros(2, float)
        prec_fun = DiagonalPreconditioner(fun, 20, 1e-2)
        prec_fun.scales = numpy.random.uniform(1,2,2)
        check_anagrad(prec_fun, x_init, 1e-5)

    def test_check_anagrad_full_prec(self):
        x_init = numpy.zeros(2, float)
        prec_fun = FullPreconditioner(fun, 20, 1e-2)
        A = numpy.random.normal(0,1,(2,2))
        A = 0.5*(A + A.transpose())
        evals, evecs = numpy.linalg.eigh(A)
        prec_fun.scales = abs(evals) + 1.0
        prec_fun.rotation = evecs
        check_anagrad(prec_fun, x_init, 1e-5)

    def test_full_prec_consitency(self):
        x_init = numpy.zeros(2, float)
        prec_fun = FullPreconditioner(fun, 20, 1e-2)
        A = numpy.random.normal(0,1,(2,2))
        A = 0.5*(A + A.transpose())
        evals, evecs = numpy.linalg.eigh(A)
        prec_fun.scales = abs(evals) + 1.0
        prec_fun.rotation = evecs

        for i in xrange(20):
            orig = numpy.random.normal(0, 1, 2)
            check = prec_fun.do(prec_fun.undo(orig))
            self.assertArraysAlmostEqual(orig, check)
            check = prec_fun.undo(prec_fun.do(orig))
            self.assertArraysAlmostEqual(orig, check)


