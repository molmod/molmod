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


from builtins import range
import numpy as np
from nose.plugins.skip import SkipTest

from molmod.test.common import BaseTestCase
from molmod import *


__all__ = ["MinimizerTestCase"]


def fun(x, do_gradient=False):
    value = 2 + np.sin(x[0]) + np.cos(x[1]) + x[0]*x[0] + x[1]*x[1] - x[0]*x[1]
    if do_gradient:
        gradient = np.array([
            np.cos(x[0]) + 2*x[0] - x[1],
            -np.sin(x[1]) + 2*x[1] - x[0],
        ])
        return value, gradient
    else:
        return value


def quad(x, do_gradient=False):
    value = ((x - 1)**2).sum()
    if do_gradient:
        gradient = 2*(x-1)
        return value, gradient
    else:
        return value


def circle1(x):
    return (x**2).sum()-4, 2*x


def circle2(x):
    x = x.copy()
    x[0] -= -1
    return (x**2).sum()-4, 2*x


class Half(object):
    def __init__(self, x0, normal):
        self.x0 = x0
        self.normal = normal

    def __call__(self, x):
        return np.dot(x - self.x0, self.normal), self.normal


class MinimizerTestCase(BaseTestCase):
    def check_min(self, x_opt, step_rms, grad_rms):
        # check if it is really the minimum by computing small displacements
        f_opt = fun(x_opt)
        for i in range(100):
            delta = np.random.normal(0, 1, 2)
            delta /= np.linalg.norm(delta)
            delta *= np.sqrt(len(delta))
            delta *= step_rms
            x_dev = x_opt + delta
            f_dev = fun(x_dev)
            self.assert_(f_opt - f_dev <= grad_rms*step_rms)
        # check if it is really the minimum by computing the eigen values of the
        # hessian
        hessian1 = compute_fd_hessian(fun, x_opt, 1e-4, anagrad=True)
        self.assert_((np.linalg.eigvalsh(hessian1) > 0).all())
        # check if it is really the minimum by computing the eigen values of the
        # hessian
        hessian2 = compute_fd_hessian(fun, x_opt, 1e-4, anagrad=False)
        self.assert_((np.linalg.eigvalsh(hessian2) > 0).all())
        # also compare the two hessians
        self.assert_(abs(hessian1 - hessian2).max() < 1e-4)

    def test_sd_golden(self):
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
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
        x_init = np.zeros(2, float)
        search_direction = QuasiNewton()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        self.check_min(minimizer.get_final(), 1e-6, 1e-6)

    def test_qn_newtong_rel_grad(self):
        x_init = np.zeros(2, float)
        search_direction = QuasiNewton()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(rel_grad_rms=1e-6, rel_grad_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, fun_margin=1e-3)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )

    def test_check_anagrad(self):
        x_init = np.zeros(2, float)
        check_anagrad(fun, x_init, 1e-5, 1e-4)

    def test_check_anagrad_diag_prec(self):
        x_init = np.zeros(2, float)
        prec_fun = DiagonalPreconditioner(fun, 20, 1e-2)
        prec_fun.scales = np.random.uniform(1,2,2)
        check_anagrad(prec_fun, x_init, 1e-5, 1e-4)
        dxs = random_unit((100, len(x_init)))*1e-4
        check_delta(prec_fun, x_init, dxs)

    def test_check_anagrad_full_prec(self):
        raise SkipTest
        x_init = np.zeros(2, float)
        prec_fun = FullPreconditioner(fun, 20, 1e-2)
        A = np.random.normal(0,1,(2,2))
        A = 0.5*(A + A.transpose())
        evals, evecs = np.linalg.eigh(A)
        prec_fun.scales = abs(evals) + 1.0
        prec_fun.rotation = evecs
        check_anagrad(prec_fun, x_init, 1e-5, 1e-4)
        dxs = random_unit((100, len(x_init)))*1e-4
        check_delta(prec_fun, x_init, dxs)

    def test_full_prec_consitency(self):
        x_init = np.zeros(2, float)
        prec_fun = FullPreconditioner(fun, 20, 1e-2)
        A = np.random.normal(0,1,(2,2))
        A = 0.5*(A + A.transpose())
        evals, evecs = np.linalg.eigh(A)
        prec_fun.scales = abs(evals) + 1.0
        prec_fun.rotation = evecs

        for i in range(20):
            orig = np.random.normal(0, 1, 2)
            check = prec_fun.do(prec_fun.undo(orig))
            self.assertArraysAlmostEqual(orig, check)
            check = prec_fun.undo(prec_fun.do(orig))
            self.assertArraysAlmostEqual(orig, check)

    def test_stop_loss_step(self):
        x_init = np.zeros(2, float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6, grad_max=3e-6, step_max=3e-6)
        stop_loss = StopLossCondition(max_iter=50, step_min=1e-2)
        minimizer = Minimizer(
            x_init, fun, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False,
        )
        assert np.sqrt((minimizer.step**2).mean()) < 1e-2

    def test_constraints1(self):
        x_init = np.array([0.1, 0.5], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=5)
        constraints = Constraints([(1, circle1)], 1e-10)
        minimizer = Minimizer(
            x_init, quad, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert np.sqrt((minimizer.gradient**2).mean()) < 1e-6

    def test_constraints2(self):
        x_init = np.array([0.1, 0.5], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=50)
        constraints = Constraints([(0, circle1)], 1e-10)
        minimizer = Minimizer(
            x_init, quad, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert np.sqrt((minimizer.gradient**2).mean()) < 1e-6

    def test_constraints3(self):
        x_init = np.array([0.1, 0.5], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=50)
        constraints = Constraints([(-1, circle1)], 1e-10)
        minimizer = Minimizer(
            x_init, quad, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert np.sqrt((minimizer.gradient**2).mean()) < 1e-6

    def test_constraints4(self):
        x_init = np.array([-2.0, 0.1], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=50)
        constraints = Constraints([(1, circle1), (-1, circle2)], 1e-10)
        minimizer = Minimizer(
            x_init, quad, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert np.sqrt((minimizer.gradient**2).mean()) < 1e-6

    def test_constraints5(self):
        def surf(x, do_gradient=False):
            value = x[0]*(1+0.5*x[1]*x[1])-x[1]*(1-0.5*x[0]*x[0])
            if do_gradient:
                gradient = np.array([
                    (1+0.5*x[1]*x[1]) + x[1]*x[0],
                    x[0]*x[1] - (1-0.5*x[0]*x[0]),
                ], float)
                return value, gradient
            else:
                return value
        x_init = np.array([3.5, 1.2], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=50)
        constraints = Constraints([
            (1, Half([0.0, 0.0], [1.0, 0.0])),
            (1, Half([0.0, 0.0], [0.0, 1.0])),
            (1, Half([5.0, 5.0], [-1.0, 0.0])),
            (1, Half([5.0, 5.0], [0.0, -1.0])),
        ], 1e-10)
        minimizer = Minimizer(
            x_init, surf, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert minimizer.x[0] == 0.0
        assert minimizer.x[1] == 5.0
        assert np.sqrt((minimizer.gradient**2).mean()) < 1e-6

    def test_constraints6(self):
        x_init = np.array([3.5, 1.2], float)
        search_direction = ConjugateGradient()
        line_search = NewtonLineSearch()
        convergence = ConvergenceCondition(grad_rms=1e-6)
        stop_loss = StopLossCondition(max_iter=50)
        constraints = Constraints([
            (1, Half([5.0, 0.0], [1.0, 0.0])),
            (-1, circle1),
        ], 1e-10)
        minimizer = Minimizer(
            x_init, quad, search_direction, line_search, convergence, stop_loss,
            anagrad=True, verbose=False, constraints=constraints
        )
        assert not minimizer.success
