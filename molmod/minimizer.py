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
"""General purpose minimization of smooth multidimensional functions

   The implementation is mainly concerned with robustness, rather than
   computational efficiency. Example usage::

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

       x_init = np.zeros(2, float)
       search_direction = ConjugateGradient()
       line_search = NewtonLineSearch()
       convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6)
       stop_loss = StopLossCondition(max_iter=50)
       minimizer = Minimizer(
           x_init, fun, search_direction, line_search, convergence, stop_loss,
           anagrad=True, verbose=True,
       )
       print "optimum", minimizer.x, fun(minimizer.x)

   The signature of the function ``fun`` must always be the same as in the
   example. The first argument. ``x`` is mandatory and contains a 1D numpy array
   with function arguments. The second argument, ``do_gradient`` is optional
   with default value ``False``.

   The returned values must also follow the same convention as in the example.
   When ``do_gradient==True``, two return values are given. The first one is
   the function value and the second one is a 1D numpy array with the partial
   derivatives of the function towards the arguments. When
   ``do_gradient==False``, only one value is returned, i.e. the function value.
"""


from __future__ import print_function, division

from builtins import range
import time

import numpy as np


__all__ = [
    "SearchDirection", "SteepestDescent", "ConjugateGradient", "QuasiNewton",
    "LineSearch", "GoldenLineSearch", "NewtonLineSearch",
    "Preconditioner", "DiagonalPreconditioner", "FullPreconditioner",
    "ConvergenceCondition", "StopLossCondition", "Constraints", "Minimizer",
    "check_anagrad", "check_delta", "compute_fd_hessian",
]


class SearchDirection(object):
    """Abstract base class for a search direction method"""
    def __init__(self):
        # 2 characters indicating the method used to determine the direction
        self.status = "??"

    def update(self, gradient, step):
        """Update the search direction given the latest gradient and step"""
        raise NotImplementedError

    def reset(self):
        """Reset the internal state of the search direction algorithm

           This implies that the next direction will be steepest descent.
        """
        raise NotImplementedError

    def is_sd(self):
        """Return True if the last direction was steepest descent"""
        raise NotImplementedError


class SteepestDescent(SearchDirection):
    """The steepest descent method.

       This method simply sets the search direction to minus the gradient. This
       method is the least efficient choice and becomes very inefficient for
       ill-conditioned problems.
    """
    def __init__(self):
        # the current conjugate direction
        self.direction = None
        SearchDirection.__init__(self)

    def update(self, gradient, step):
        """Update the search direction given the latest gradient and step"""
        self.direction = -gradient
        self.status = "SD"

    def reset(self):
        """Reset the internal state of the search direction algorithm

           Does nothing in this case.
        """
        pass

    def is_sd(self):
        """Return True if the last direction was steepest descent

           Always returns True in this case.
        """
        return True


class ConjugateGradient(SearchDirection):
    """The conjugate gradient method

       This method is always superior to the steepest descent method in
       practical applications. An automatic reset mechanism reverts the search
       direction to the steepest descent when beta becomes negative.
    """
    def __init__(self):
        # the current conjugate direction
        self.direction = None
        # the current gradient
        self.gradient = None
        # the previous gradient
        self.gradient_old = None
        SearchDirection.__init__(self)

    def update(self, gradient, step):
        """Update the search direction given the latest gradient and step"""
        do_sd = self.gradient_old is None
        self.gradient_old = self.gradient
        self.gradient = gradient
        if do_sd:
            self._update_sd()
        else:
            self._update_cg()

    def reset(self):
        """Reset the internal state of the search direction algorithm"""
        self.gradient_old = None

    def is_sd(self):
        """Return True if the last direction was steepest descent"""
        return self.status == "SD"

    def _update_cg(self):
        """Update the conjugate gradient"""
        beta = self._beta()
        # Automatic direction reset
        if beta < 0:
            self.direction = -self.gradient
            self.status = "SD"
        else:
            self.direction = self.direction * beta - self.gradient
            self.status = "CG"

    def _beta(self):
        # Polak-Ribiere
        return (
            np.dot(self.gradient, self.gradient - self.gradient_old) /
            np.dot(self.gradient_old, self.gradient_old)
        )

    def _update_sd(self):
        """Reset the conjugate gradient to the normal gradient"""
        self.direction = -self.gradient
        self.status = "SD"


class QuasiNewton(SearchDirection):
    """The quasi Newton method"""
    def __init__(self):
        self.inv_hessian = None
        self.gradient = None
        self.old_gradient = None
        SearchDirection.__init__(self)

    def update(self, gradient, step):
        """Update the search direction given the latest gradient and step"""
        self.old_gradient = self.gradient
        self.gradient = gradient
        N = len(self.gradient)
        if self.inv_hessian is None:
            # update the direction
            self.direction = -self.gradient
            self.status = "SD"
            # new guess of the inverse hessian
            self.inv_hessian = np.identity(N, float)
        else:
            # update the direction
            self.direction = -np.dot(self.inv_hessian, self.gradient)
            self.status = "QN"
            # new guess of the inverse hessian (BFGS)
            y = self.gradient - self.old_gradient
            s = step
            sy = abs(np.dot(s, y))+1e-5
            A = np.outer(-y/sy, s)
            A.ravel()[::N+1] += 1
            self.inv_hessian = (
                np.dot(np.dot(A.transpose(), self.inv_hessian), A) +
                np.outer(s/sy, s)
            )

    def reset(self):
        """Reset the internal state of the search direction algorithm"""
        self.inv_hessian = None
        self.gradient = None
        self.old_gradient = None


    def is_sd(self):
        """Return True if the last direction was steepest descent"""
        return self.status == "SD"


phi = 0.5*(1+np.sqrt(5))


class LineSearch(object):
    """Abstract base class for a line search"""

    def __init__(self, qmax=None):
        """
           Optional argument:
            | ``qmax``  --  The maximum step size of a line search
        """
        self.qmax = qmax

    def limit_step(self, step):
        """Clip the a step within the maximum allowed range"""
        if self.qmax is None:
            return step
        else:
            return np.clip(step, -self.qmax, self.qmax)

    def __call__(self, fun, f0, initial_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
            | ``fun``  --  function to minimize (one-dimensional)
            | ``f0``   --  the function value at the starting point q=0"
            | ``initial_step_size``  --  a guess of the order of magnitude of
                                         step size to be found.
            | ``epsilon``  --  a value that is small compared to
                               initial_step_size

           Returns:
            | ``success``  --  a boolean indicating that the line search
                               resulted in an improved solution
            | ``wolfe``  --  a boolean indicating that the new solution
                             satisfies Wolfe conditions
            | ``qopt``  --  the position of the new solution on the line
            | ``fopt``  --  the corresponding function value
        """
        raise NotImplementedError


class GoldenLineSearch(LineSearch):
    """The golden section line search algorithm"""

    def __init__(self, qtol, qmax=None, max_iter=None):
        """
           Argument:
            | ``qtol``  --  The threshold for displacements along the line.
                            (When displacements become smaller than qtol, we
                            assume convergence.)
           Optional arguments
            | ``qmax``  --  The maximum step size of a line search
            | ``max_iter``  --  the maximum number of iteration for the line
                                search (only applies to the bracketing part)
        """
        if qtol is None:
            raise ValueError("No stop condition is specified")
        self.qtol = qtol
        self.max_iter = max_iter
        LineSearch.__init__(self, qmax)
        self.num_bracket = 0
        self.num_golden = 0

    def __call__(self, fun, initial_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
            | ``fun``  --  function to minimize (one-dimensional)
            | ``f0``   --  the function value at the starting point q=0"
            | ``initial_step_size``  --  a guess of the order of magnitude of
                                         step size to be found. This is used to
                                         bracket the minimum.
            | ``epsilon``  --  a value that is small compared to
                               initial_step_size

           Returns:
            | ``success``  --  a boolean indicating that the line search
                               resulted in an improved solution
            | ``wolfe``  --  a boolean indicating that the new solution
                             satisfies Wolfe conditions
            | ``qopt``  --  the position of the new solution on the line
            | ``fopt``  --  the corresponding function value

           P.S. The wolfe parameter is always True, but this aspect is not
           guaranteed to be correct. Never use the :class:`GoldenLineSearch` in
           combination with a quasi Newton method.
        """
        # bracket the minimum
        f0 = fun(0.0)
        triplet = self._bracket(initial_step_size, f0, fun)
        if triplet is None:
            return False, False, 0.0, f0
        # do a golden section optimization
        qopt, fopt = self._golden(triplet, fun)
        qopt = self.limit_step(qopt)
        fopt = fun(qopt)
        return True, True, qopt, fopt

    def _bracket(self, qinit, f0, fun):
        """Find a bracket that does contain the minimum"""
        self.num_bracket = 0
        qa = qinit
        fa = fun(qa)
        counter = 0
        if fa >= f0:
            while True:
                self.num_bracket += 1
                #print "    bracket shrink"
                qb, fb = qa, fa
                qa /= 1+phi
                fa = fun(qa)
                if qa < self.qtol:
                    return
                if fa < f0:
                    return (0, f0), (qa, fa), (qb, fb)
                counter += 1
                if self.max_iter is not None and counter > self.max_iter:
                    return
        else:
            self.num_bracket += 1
            #print "    bracket grow1"
            qb, fb = qa, fa
            qa *= (1+phi)
            fa = fun(qa)
            if fa >= fb:
                return (0, f0), (qb, fb), (qa, fa)
            while True:
                self.num_bracket += 1
                #print "    bracket grow2"
                qc, fc = qb, fb
                qb, fb = qa, fa
                qa = qb*(1+phi) - qc
                fa = fun(qa)
                if fa >= fb:
                    return (qc, fc), (qb, fb), (qa, fa)
                counter += 1
                if self.max_iter is not None and counter > self.max_iter:
                    return

    def _golden(self, triplet, fun):
        """Reduce the size of the bracket until the minimum is found"""
        self.num_golden = 0
        (qa, fa), (qb, fb), (qc, fc) = triplet
        while True:
            self.num_golden += 1
            qd = qa + (qb-qa)*phi/(1+phi)
            fd = fun(qd)
            if fd < fb:
                #print "golden d"
                (qa, fa), (qb, fb) = (qb, fb), (qd, fd)
            else:
                #print "golden b"
                (qa, fa), (qc, fc) = (qd, fd), (qa, fa)
            if abs(qa-qb) < self.qtol:
                return qc, fc


class NewtonLineSearch(LineSearch):
    """The Newton line search algorithm

       When the curvature is negative, a steepest descent step is tried, using
       the step size from the previous step multiplied by 1.5. If the new
       function value is higher, the step size is reduced by a factor two. The
       latter is repeated at most max_iter times. If no lower value is found,
       the line search fails. (This is known is the back tracking algorithm)

       When the curvature is positive, Newton step are performed. When the
       function or the absolute value of the derivative at a new point
       increases, the procedure is interupted and the last descent point is
       used. When there is no last descent point, back tracking is used. The
       Wolfe conditions are used to determine the convergence of the line
       search. At most max_iter Newton steps are allowed.
    """
    def __init__(self, c1=1e-4, c2=1e-1, max_iter=5, qmax=None):
        """
           Optional arguments:
            | ``c1``  --  The coefficient in the first Wolfe condition
                          (sufficient decrease of the function) [default=1e-4]
            | ``c2``  --  The coefficient in the second Wolfe condition
                          (sufficient decrease of the derivative)
                          [default=1e-1]. the default is optimal for the
                          conjugate gradient method
            | ``max_iter``  --  the maximum number of iterations in the line
                                search.
            | ``qmax``  --  The maximum step size of a line search
        """
        self.c1 = c1
        self.c2 = c2
        self.max_iter = max_iter
        LineSearch.__init__(self, qmax)

    def __call__(self, fun, initial_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
            | ``fun``  --  function to minimize (one-dimensional)
            | ``f0``   --  the function value at the starting point q=0"
            | ``initial_step_size``  --  a guess of the order of magnitude of
                                         step size to be found. This is used
                                         in case the default newton line search
                                         fails and when the routine reverts to
                                         backtracking as a last resort.
            | ``epsilon``  --  a value that is small compared to
                               initial_step_size

           Returns:
            | ``success``  --  a boolean indicating that the line search
                               resulted in an improved solution
            | ``wolfe``  --  a boolean indicating that the new solution
                             satisfies Wolfe conditions
            | ``qopt``  --  the position of the new solution on the line
            | ``fopt``  --  the corresponding function value
        """
        f0, g0 = fun(0.0, do_gradient=True)
        # Approximate the second order derivative with symmetric finite differences.
        gl = fun(-0.5*epsilon, do_gradient=True)[1]
        gh = fun(+0.5*epsilon, do_gradient=True)[1]
        h0 = (gh-gl)/epsilon

        # If the line function is hollow, then try a newton step.
        if h0 > 0:
            # the initial point becomes the old point.
            q1, f1, g1, h1 = 0.0, f0, g0, h0
            counter = 0
            wolfe = False
            while True:
                q2 = q1-g1/h1
                f2, g2 = fun(q2, do_gradient=True)
                if abs(g2) > abs(g1) or f2 > f1:
                    # If the gradient or the function increases in absolute
                    # value, the newton step is clearly not working.
                    break
                if self.qmax is not None and q2 > self.qmax:
                    # This step is going to far, which we consider unsafe. This
                    # is probably because h0 is too close to zero.
                    break
                counter += 1
                if self.max_iter is not None and counter > self.max_iter:
                    break
                # The new point becomes the old point
                q1, f1, g1 = q2, f2, g2
                del q2
                del f2
                del g2
                # Check for Wolfe conditions.
                if f1 >= f0 + self.c1*abs(g0*q1):
                    break
                if f1 <= f0 - self.c1*abs(g0*q1) and abs(g1) <= abs(g0*self.c2):
                    # The Wolfe conditions are satisfied.
                    wolfe = True
                    break
                # Make a new estimate of the second order derivative at q1
                gl = fun(q1-0.5*epsilon, do_gradient=True)[1]
                gh = fun(q1+0.5*epsilon, do_gradient=True)[1]
                h1 = (gh-gl)/epsilon
            if counter > 0 and f1 <= f0:
                # If at least one step is taken in the newton procedure, we are
                # happy.
                return True, wolfe, q1, f1
            else:
                # even the first newton step failed, revert to back tracking
                pass

        # This is the fall-back part for when the Newton step did not work.
        # Simple back-tracking is carried out with tau = 0.5.
        if self.qmax is None:
            qmax = initial_step_size*1.5
        else:
            qmax = min(initial_step_size*1.5, self.qmax)
        q1 = -np.sign(g0)*qmax
        counter = 0.0
        while True:
            f1 = fun(q1)
            if f1 < f0:
                return True, False, q1, f1
            q1 *= 0.5
            counter += 1
            if self.max_iter is not None and counter > self.max_iter:
                # had enough iterations, line search fails
                return False, False, 0.0, f0


class Preconditioner(object):
    """Base class for preconditioners

       A preconditioner is a (linear) transformation of the unknowns to a new
       basis in which the Hessian of the minimum becomes a better-conditioned
       matrix. In these new coordinates the convergence of the minimizer will
       be faster. Picking the right preconditioner is a matter of experience.
       One must balance the extra computational cost of the preconditioner
       against the gains in computational cost because of the reduced number of
       iterations in the minimizer.

       The preconditioners in this package act as wrappers around the function
       to be optimized. One just replaces a function by the preconditioner in
       the constructor of the Minimizer object. E.g. ::

         >>> Minimizer(fun, ...)

       becomes::

         >>> Minimizer(SomePreconditioner(fun, ...), ...)

       Also note that the convergence and stop loss conditions are not affected
       by the preconditioner. They get the gradient and step in original
       coordinates.
    """
    def __init__(self, fun, each, grad_rms):
        """
           Arguments:
            | ``fun``  --  the function whose arguments must be transformed
            | ``each``  --  update the linear transformation after each 'each'
                            minimizer steps without updates
            | ``grad_rms``  --  only update when the rms value of the gradient
                                (in the original coordinates) is below this
                                threshold
        """
        self.fun = fun
        self.each = each
        self.grad_rms = grad_rms

        self.last_update = 0

    def __call__(self, x_prec, do_gradient=False):
        """The actual wrapper around the function call.

           Arguments:
            | ``x_prec``  --  the unknowns in preconditioned coordinates
            | ``do_gradient``  --  if True, the gradient is also computed and
                                   transformed to preconditioned coordinates

           Note that this implementation assumes that the preconditioner is a
           linear transformation.
        """
        if do_gradient:
            f, g = self.fun(self.undo(x_prec), do_gradient=True)
            return f, self.undo(g)
        else:
            return self.fun(self.undo(x_prec))

    def update(self, counter, f, x_orig, gradient_orig):
        """Perform an update of the linear transformation

           Arguments:
            | ``counter``  --  the iteration counter of the minimizer
            | ``f``  --  the function value at ``x_orig``
            | ``x_orig``  --  the unknowns in original coordinates
            | ``gradient_orig``  --  the gradient in original coordinates

           Return value:
            | ``do_update``  --  True when an update is required.

           Derived classes must call this method to test of the preconditioner
           requires updating. Derived classes must also return this boolean
           to their caller.
        """
        if counter - self.last_update > self.each:
            grad_rms = np.sqrt((gradient_orig**2).mean())
            if grad_rms < self.grad_rms:
                self.last_update = counter
                return True
        return False

    def do(self, x_orig):
        """Transform the unknowns to preconditioned coordinates

           This method also transforms the gradient to original coordinates
        """
        raise NotImplementedError

    def undo(self, x_prec):
        """Transform the unknowns to original coordinates

           This method also transforms the gradient to preconditioned coordinates
        """
        raise NotImplementedError


class DiagonalPreconditioner(Preconditioner):
    """The diagonal preconditioner

       This preconditioner derives a diagonal transformation based on a finite
       difference approximation of the diagonal elements of the Hessian. The
       trasnformation is such that these diagonal elements would become equal
       after the transformation. This type of preconditioner is especially
       usefull when the unknowns have different units. In many cases this
       preconditioner is a good trade off betweem accelerated convergence and
       extra cost. In particular when it is combined with a conjugate gradient
       minimizer, it can be more effecient that a quasi Newton method.

       (For more general info on preconditioners, read the doc string of the
       Preconditioner base class.)
    """
    def __init__(self, fun, each, grad_rms, epsilon=1e-3, scale_limit=0.0):
        """
           Arguments:
            | ``fun``  --  the function whose arguments must be transformed
            | ``each``  --  update the linear transformation after each 'each'
                            minimizer steps without updates
            | ``grad_rms``  --  only update when the rms value of the gradient
                                (in the original coordinates) is below this
                                threshold

           Optional argument:
            | ``epsilon``  --  a small scalar used for the finite differences
                               (taken in previous preconditioned coordinates)
                               [default=1e-3]
            | ``scale_limit``  --  scales smaller than scale_limit times the
                                   largest scale are fixed to scale_limit times
                                   the largest scale
        """
        self.epsilon = epsilon
        self.scale_limit = scale_limit
        Preconditioner.__init__(self, fun, each, grad_rms)
        self.scales = None

    def update(self, counter, f, x_orig, gradient_orig):
        """Perform an update of the linear transformation

           Arguments:
            | ``counter``  --  the iteration counter of the minimizer
            | ``f``  --  the function value at ``x_orig``
            | ``x_orig``  --  the unknowns in original coordinates
            | ``gradient_orig``  --  the gradient in original coordinates

           Return value:
            | ``done_update``  --  True when an update has been done

           The minimizer must reset the search direction method when an updated
           has been done.
        """
        do_update = Preconditioner.update(self, counter, f, x_orig, gradient_orig)
        if do_update:
            # determine a new preconditioner
            N = len(x_orig)
            if self.scales is None:
                self.scales = np.ones(N, float)
            for i in range(N):
                epsilon = self.epsilon/self.scales[i]
                xh = x_orig.copy()
                xh[i] += 0.5*epsilon
                fh = self.fun(xh)
                xl = x_orig.copy()
                xl[i] -= 0.5*epsilon
                fl = self.fun(xl)
                curv = (fh+fl-2*f)/epsilon**2
                self.scales[i] = np.sqrt(abs(curv))
            if self.scales.max() <= 0:
                self.scales = np.ones(N, float)
            else:
                self.scales /= self.scales.max()
                self.scales[self.scales<self.scale_limit] = self.scale_limit
        return do_update

    def do(self, x_orig):
        """Transform the unknowns to preconditioned coordinates

           This method also transforms the gradient to original coordinates
        """
        if self.scales is None:
            return x_orig
        else:
            return x_orig*self.scales

    def undo(self, x_prec):
        """Transform the unknowns to original coordinates

           This method also transforms the gradient to preconditioned coordinates
        """
        if self.scales is None:
            return x_prec
        else:
            return x_prec/self.scales


class FullPreconditioner(Preconditioner):
    """The full preconditioner

       This preconditioner is a bit experimental. The transformation is such
       that the hessian in the new coordinates becomes a constant matrix,
       i.e. diagonal with all elements the same.
    """
    def __init__(self, fun, each, grad_rms, epsilon=1e-3):
        """
           Arguments:
            | ``fun``  --  the function whose arguments must be transformed
            | ``each``  --  update the linear transformation after each 'each'
                            minimizer steps without updates
            | ``grad_rms``  --  only update when the rms value of the gradient
                                (in the original coordinates) is below this
                                threshold

           Optional argument:
            | ``epsilon``  --  a small scalar used for the finite differences
                               (taken in original coordinates) [default=1e-3]
        """
        self.epsilon = epsilon
        Preconditioner.__init__(self, fun, each, grad_rms)
        self.scales = None
        self.rotation = None

    def update(self, counter, f, x_orig, gradient_orig):
        """Perform an update of the linear transformation

           Arguments:
            | ``counter``  --  the iteration counter of the minimizer
            | ``f``  --  the function value at ``x_orig``
            | ``x_orig``  --  the unknowns in original coordinates
            | ``gradient_orig``  --  the gradient in original coordinates

           Return value:
            | ``done_update``  --  True when an update has been done

           The minimizer must reset the search direction method when an updated
           has been done.
        """
        if Preconditioner.update(self, counter, f, x_orig, gradient_orig):
            # determine a new preconditioner
            hessian = compute_fd_hessian(self.fun, x_orig, self.epsilon)
            evals, evecs = np.linalg.eigh(hessian)
            self.scales = np.sqrt(abs(evals))+self.epsilon
            self.rotation = evecs
            return True
        return False

    def do(self, x_orig):
        """Transform the unknowns to preconditioned coordinates

           This method also transforms the gradient to original coordinates
        """
        if self.scales is None:
            return x_orig
        else:
            return np.dot(self.rotation.transpose(), x_orig)*self.scales

    def undo(self, x_prec):
        """Transform the unknowns to original coordinates

           This method also transforms the gradient to preconditioned coordinates
        """
        if self.scales is None:
            return x_prec
        else:
            return np.dot(self.rotation, x_prec/self.scales)


class ConvergenceCondition(object):
    """Callable object that tests the convergence of the minimizer"""

    def __init__(self, step_rms=None, step_max=None, grad_rms=None, grad_max=None, rel_grad_rms=None, rel_grad_max=None):
        """
           Optional arguments:
            | ``step_rms``  --  threshold for the RMS value of the step vector
                                in the iterative minimizer
            | ``step_max``  --  threshold for the maximum component of the step
                                vector in the iterative minimizer
            | ``grad_rms``  --  threshold for the RMS value of the gradient
                                components of the function to be minimized
            | ``grad_max``  --  threshold for the maximum value of the gradient
                                components of the function to be minimized
            | ``rel_grad_rms``  --  threshold for the RMS value of the gradient
                                    components of the function to be minimized,
                                    divided by the function value
            | ``rel_grad_max``  --  threshold for the maximum value of the gradient
                                    components of the function to be minimized,
                                    divided by the function value

           Only the present arguments define when the minimization has
           converged. All actual values must go below the given thresholds.
        """
        if (step_rms is None and step_max is None and
            grad_rms is None and grad_max is None and
            rel_grad_rms is None and rel_grad_max is None):
            raise ValueError("Some convergence criteria must be specified")
        self.step_rms = step_rms
        self.step_max = step_max
        self.grad_rms = grad_rms
        self.grad_max = grad_max
        self.rel_grad_rms = rel_grad_rms
        self.rel_grad_max = rel_grad_max

    def get_header(self):
        """Returns the header for screen logging of the minimization"""
        result = " "
        if self.step_rms is not None:
            result += "    Step RMS"
        if self.step_max is not None:
            result += "    Step MAX"
        if self.grad_rms is not None:
            result += "    Grad RMS"
        if self.grad_max is not None:
            result += "    Grad MAX"
        if self.rel_grad_rms is not None:
            result += "  Grad/F RMS"
        if self.rel_grad_max is not None:
            result += "  Grad/F MAX"
        return result

    def __call__(self, grad, step, f):
        """Return True when the minimizer has converged

           Arguments:
            | ``grad``  --  The gradient at the current point.
            | ``step``  --  The last step vector.
            | ``f`` -- The last function value.
        """
        stop = True
        status = ""
        red = "\033[0;31m"
        green = "\033[0;32m"
        reset = "\033[m"

        def check_threshold(measure, threshold, stop, status):
            if threshold is not None:
                #measure = abs(grad).max()
                if measure > threshold:
                    color = red
                    stop = False
                else:
                    color = green
                status += "%s% 9.3e  %s" % (color, measure, reset)
            return stop, status

        def safe_divide(a, b):
            if b == 0.0:
                return 0.0
            else:
                return a/b

        stop, status = check_threshold(np.sqrt((step**2).mean()), self.step_rms, stop, status)
        stop, status = check_threshold(abs(step).max(), self.step_max, stop, status)
        stop, status = check_threshold(np.sqrt((grad**2).mean()), self.grad_rms, stop, status)
        stop, status = check_threshold(abs(grad).max(), self.grad_max, stop, status)
        stop, status = check_threshold(safe_divide(np.sqrt((grad**2).mean()), f), self.rel_grad_rms, stop, status)
        stop, status = check_threshold(safe_divide(abs(grad).max(), f), self.rel_grad_max, stop, status)

        return stop, status


class StopLossCondition(object):
    """Callable object that checks if minimizer has lost track"""
    def __init__(self, max_iter=None, fun_margin=None, grad_margin=None, step_min=None):
        """
           Optional arguments:
            | ``max_iter``  --  the maximum number of iterations allowed
            | ``fun_margin``  --  if the function to be minimized goes above the
                                  lowest value so far plus this margin, the
                                  minimization is aborted
            | ``grad_margin``  --  if the RMS value of the gradient components
                                   goes above the lowest value plus this
                                   threshold, the minimization is aborted
            | ``step_min``  --  If the RMS step size drops below this margin, the
                                optimization is interrupted.

           Only the present arguments define when the minimization has lost
           track.
        """
        self.max_iter = max_iter
        self.fun_margin = fun_margin
        self.grad_margin = grad_margin
        self.step_min = step_min

        self.reset()

    def reset(self):
        self.fn_lowest = None
        self.grad_rms_lowest = None

    def __call__(self, counter, fn, gradient, step):
        """Return True when the minimizer has lost track"""
        if self.max_iter is not None and counter >= self.max_iter:
            return True

        if self.fun_margin is not None:
            if self.fn_lowest is None or fn < self.fn_lowest:
                self.fn_lowest = fn
            elif fn > self.fn_lowest + self.fun_margin:
                return True

        if self.grad_margin is not None:
            grad_rms = np.sqrt((gradient**2).mean())
            if self.grad_rms_lowest is None or grad_rms < self.grad_rms_lowest:
                self.grad_rms_lowest = grad_rms
            elif grad_rms > self.grad_rms_lowest + self.grad_margin:
                return True

        if self.step_min is not None:
            step_rms = np.sqrt((step**2).mean())
            if self.step_min > step_rms:
                return True

        # all is fine
        return False


class LineWrapper(object):
    """A configurable line function"""
    def __init__(self, fun, anagrad, epsilon):
        """
           Argument:
            | ``fun``  --  a multivariate function, see below
            | ``anagrad``  --  boolean that indicates if fun supports analytical
                               gradients
            | ``epsilon``  --  a small scalar used for finite differences

           The function ``fun`` takes a mandatory argument ``x`` and an optional
           argument ``do_gradient``:
            | ``x``  --  the arguments of the function to be tested
            | ``do_gradient``  --  when False, only the function value is
                                   returned. when True, a 2-tuple with the
                                   function value and the gradient are returned
                                   [default=False]
        """
        self.fun = fun
        self.anagrad = anagrad
        self.epsilon = epsilon
        self.axis = None
        self.x0 = None

    def configure(self, x0, axis):
        """Configure the 1D function for a line search

           Arguments:
             x0  --  the reference point (q=0)
             axis  --  a unit vector in the direction of the line search
        """
        self.x0 = x0
        self.axis = axis

    def __call__(self, q, do_gradient=False):
        x = self.x0 + self.axis*q
        if do_gradient:
            if self.anagrad:
                f, g = self.fun(x, do_gradient=True)
                return f, np.dot(g, self.axis)
            else:
                fh = self.fun(x + (0.5*self.epsilon) * self.axis)
                fl = self.fun(x - (0.5*self.epsilon) * self.axis)
                return self.fun(x), (fh - fl)/self.epsilon
        else:
            return self.fun(x)


class FunWrapper(object):
    """Wrapper to compute the function and its gradient"""
    def __init__(self, fun, anagrad, epsilon):
        """
           Arguments:
            | ``fun``  --  a multivariate function that can also compute
                           analytical derivatives, see below
            | ``anagrad``  --  boolean that indicates if fun supports analytical
                               gradients
            | ``epsilon``  --  a small scalar used for finite differences

           The function ``fun`` takes a mandatory argument ``x`` and an optional
           argument ``do_gradient``:
            | ``x``  --  the arguments of the function to be tested
            | ``do_gradient``  --  when False, only the function value is
                                   returned. when True, a 2-tuple with the
                                   function value and the gradient are returned
                                   [default=False]
        """
        self.fun = fun
        self.anagrad = anagrad
        self.epsilon = epsilon

    def __call__(self, x, do_gradient=False):
        if do_gradient:
            if self.anagrad:
                return self.fun(x, do_gradient=True)
            else:
                g = np.zeros(x.shape)
                for j in range(len(x)):
                    xh = x.copy()
                    xh[j] += 0.5*self.epsilon
                    xl = x.copy()
                    xl[j] -= 0.5*self.epsilon
                    g[j] = (self.fun(xh) - self.fun(xl))/self.epsilon
                return self.fun(x), g
        else:
            return self.fun(x)


class ConstraintError(Exception):
    pass


class Constraints(object):
    '''Algorithm to apply half-open and convential constraints during minimization.'''
    def __init__(self, equations, threshold, rcond1=1e-10, max_iter=100):
        '''
           The constraint solver internally works with a constraint cost
           function, defined as the squared sum of the constraint functions.
           The constraints are satisfied by bringing the cost function to zero.
           This is done in iterative fashion. At each iteration, two attempts
           are made to lower the constraint cost function:

           1) Take a Levenberg-Marquardt-like step.
           2) If (1) fails or is too slow, take a step to fix only one of the
              constraints, i.e. the constraint which has the largest mismatch.

           Arguments:
            | ``equations`` -- a list of (sign,equation) pairs. sign can be +1,
                               0 of -1. equation is a function with one
                               argument: the vector of unknowns in the
                               minimizer. It returns the value of the constraint
                               function and the gradient of that function. If
                               sign is +1, the parameters will be forced in the
                               region where the constraint function is positive.
                               (Similar for -1, constraint function is forced to
                               be negative.) When the sign is 0, the constraint
                               function is forced to be zero.
            | ``threshold`` -- The acceptable allowed deviation from the
                               constraints. The deviation is defined as the
                               euclidean norm of the (active) constraint
                               functions.

           Optional arguments:
            | ``rcond1`` -- During the iterative solution of the constraint
                            equations in the shake algorithm, it may happen
                            that an ill-conditioned set of equations must be
                            solved. In that case rcond1 is the first ridge
                            parameter used to regularize these equations. If
                            needed, the ridge parameter is multiplied by 10
                            until a better fit of the constraints is found.
            | ``max_iter`` -- The maximum number of iterations in the shake
                              algorithm. This is used in several functions.
        '''
        self.equations = equations
        self.lock = np.zeros(len(equations), bool)
        self.threshold = threshold
        self.rcond1 = rcond1
        self.max_iter = max_iter

    def _compute_equations(self, x, verbose=False):
        '''Compute the values and the normals (gradients) of active constraints.

           Arguments:
            | ``x`` -- The unknowns.
        '''
        # compute the error and the normals.
        normals = []
        values = []
        signs = []
        error = 0.0
        if verbose:
            print()
            print(' '.join('% 10.3e' % val for val in x), end=' ')
            active_str = ''
        for i, (sign, equation) in enumerate(self.equations):
            value, normal = equation(x)
            if (i < len(self.lock) and self.lock[i]) or \
               (sign==-1 and value > -self.threshold) or \
               (sign==0) or (sign==1 and value < self.threshold):
                values.append(value)
                normals.append(normal)
                signs.append(sign)
                error += value**2
                if verbose:
                    active_str += 'X'
                if i < len(self.lock):
                    self.lock[i] = True
            elif verbose:
                active_str += '-'
        error = np.sqrt(error)
        normals = np.array(normals, float)
        values = np.array(values, float)
        signs = np.array(signs, int)
        if verbose:
            print('[%s]' % active_str, end=' ')
            if error < self.threshold:
                print('OK')
            else:
                print('%.5e' % error)
        return normals, values, error, signs

    def _rough_shake(self, x, normals, values, error):
        '''Take a robust, but not very efficient step towards the constraints.

           Arguments:
            | ``x`` -- The unknowns.
            | ``normals`` -- A numpy array with the gradients of the active
                             constraints. Each row is one gradient.
            | ``values`` -- A numpy array with the values of the constraint
                            functions.
            | ``error`` -- The square root of the constraint cost function.
        '''
        counter = 0
        while error > self.threshold and counter < self.max_iter:
            dxs = []
            for i in range(len(normals)):
                dx = -normals[i]*values[i]/np.dot(normals[i], normals[i])
                dxs.append(dx)
            dxs = np.array(dxs)
            dx = dxs[abs(values).argmax()]
            x = x+dx
            self.lock[:] = False
            normals, values, error = self._compute_equations(x)[:-1]
            counter += 1
        return x, normals, values, error

    def _fast_shake(self, x, normals, values, error):
        '''Take an efficient (not always robust) step towards the constraints.

           Arguments:
            | ``x`` -- The unknowns.
            | ``normals`` -- A numpy array with the gradients of the active
                             constraints. Each row is one gradient.
            | ``values`` -- A numpy array with the values of the constraint
                            functions.
            | ``error`` -- The square root of the constraint cost function.
        '''
        # filter out the degrees of freedom that do not feel the constraints.
        mask = (normals!=0).any(axis=0) > 0
        normals = normals[:,mask]
        # Take a step to lower the constraint cost function. If the step is too
        # large, it is reduced iteratively towards a small steepest descent
        # step. This is very similar to the Levenberg-Marquardt algorithm.
        # Singular Value decomposition is used to make this procedure
        # numerically more stable and efficient.
        U, S, Vt = np.linalg.svd(normals, full_matrices=False)
        rcond = None
        counter = 0
        while True:
            if rcond is None:
                rcond = 0.0
            elif rcond == 0.0:
                rcond = self.rcond1
            else:
                rcond *= 10
            # perform the least-norm correction
            Sinv = (S**2+rcond)
            if Sinv.max() == 0.0:
                continue
            Sinv = S/Sinv
            # compute the step
            dx = -np.dot(Vt.transpose(), np.dot(U.transpose(), values)*Sinv)
            new_x = x.copy()
            new_x[mask] += dx
            # try the step
            new_normals, new_values, new_error = self._compute_equations(new_x)[:-1]
            if new_error < 0.9*error:
                # Only if it decreases the constraint cost sufficiently, the
                # step is accepted. This routine is pointless of it converges
                # slowly.
                return new_x, new_normals, new_values, new_error
            elif abs(dx).sum() < self.threshold:
                # If the step becomes too small, then give up.
                break
            elif counter > self.max_iter:
                raise ConstraintError('Exceeded maximum number of shake iterations.')
            counter += 1


    def free_shake(self, x):
        '''Brings unknowns to the constraints.

           Arguments:
            | ``x`` -- The unknowns.
        '''
        self.lock[:] = False
        normals, values, error = self._compute_equations(x)[:-1]
        counter = 0
        while True:
            if error <= self.threshold:
                break
            # try a well-behaved move to the constrains
            result = self._fast_shake(x, normals, values, error)
            counter += 1
            if result is not None:
                x, normals, values, error = result
            else:
                # well-behaved move is too slow.
                # do a cumbersome move to satisfy constraints approximately.
                x, normals, values, error = self._rough_shake(x, normals, values, error)
                counter += 1
            # When too many iterations are required, just give up.
            if counter > self.max_iter:
                raise ConstraintError('Exceeded maximum number of shake iterations.')
        return x, counter, len(values)

    def safe_shake(self, x, fun, fmax):
        '''Brings unknowns to the constraints, without increasing fun above fmax.

           Arguments:
            | ``x`` -- The unknowns.
            | ``fun`` -- The function being minimized.
            | ``fmax`` -- The highest allowed value of the function being
                          minimized.

           The function ``fun`` takes a mandatory argument ``x`` and an optional
           argument ``do_gradient``:
            | ``x``  --  the arguments of the function to be tested
            | ``do_gradient``  --  when False, only the function value is
                                   returned. when True, a 2-tuple with the
                                   function value and the gradient are returned
                                   [default=False]
        '''
        self.lock[:] = False
        def extra_equation(xx):
            f, g = fun(xx, do_gradient=True)
            return (f-fmax)/abs(fmax), g/abs(fmax)
        self.equations.append((-1,extra_equation))
        x, shake_counter, constraint_couter = self.free_shake(x)
        del self.equations[-1]
        return x, shake_counter, constraint_couter

    def project(self, x, vector):
        '''Project a vector (gradient or direction) on the active constraints.

           Arguments:
            | ``x`` -- The unknowns.
            | ``vector`` -- A numpy array with a direction or a gradient.

           The return value is a gradient or direction, where the components
           that point away from the constraints are projected out. In case of
           half-open constraints, the projection is only active of the vector
           points into the infeasible region.
        '''
        scale = np.linalg.norm(vector)
        if scale == 0.0:
            return vector
        self.lock[:] = False
        normals, signs = self._compute_equations(x)[::3]
        if len(normals) == 0:
            return vector

        vector = vector/scale
        mask = signs == 0
        result = vector.copy()
        changed = True
        counter = 0
        while changed:
            changed = False
            y = np.dot(normals, result)
            for i, sign in enumerate(signs):
                if sign != 0:
                    if sign*y[i] < -self.threshold:
                        mask[i] = True
                        changed = True
                    elif mask[i] and np.dot(normals[i], result-vector) < 0:
                        mask[i] = False
                        changed = True

            if mask.any():
                normals_select = normals[mask]
                y = np.dot(normals_select, vector)
                U, S, Vt = np.linalg.svd(normals_select, full_matrices=False)
                if S.min() == 0.0:
                    Sinv = S/(S**2+self.rcond1)
                else:
                    Sinv = 1.0/S
                result = vector - np.dot(Vt.transpose(), np.dot(U.transpose(), y)*Sinv)
            else:
                result = vector.copy()

            if counter > self.max_iter:
                raise ConstraintError('Exceeded maximum number of shake iterations.')
            counter += 1

        return result*scale


class Minimizer(object):
    """A flexible multivariate minimizer

       The minimizer searches (in principle) for the 'nearest' local minimum.
    """

    def __init__(self, x_init, fun, search_direction, line_search,
                 convergence_condition, stop_loss_condition, anagrad=False,
                 epsilon=1e-6, verbose=True, callback=None,
                 initial_step_size=1.0, constraints=None, debug_line=False):
        """
           Arguments:
            | ``x_init``  --  the initial guess for the minimum
            | ``fun``  --  function to be minimized (see below)
            | ``search_direction``  --  a SearchDirection object
            | ``line_search``  --  a LineSearch object
            | ``convergence_condition``  --  a ConvergenceCondition object
            | ``stop_loss_condition``  --  a StopLossCondition object

           Optional arguments:
            | ``anagrad``  --  when set to True, analytical gradients are used
            | ``epsilon``  --  a small value compared to expected changes in the
                               unknowns [default=1e-6]. it is used to compute
                               finite differences.
            | ``verbose``  --  print progress information on screen
                               [default=True]
            | ``callback``  --  optional callback routine after each iteration.
                                the callback routine gets the minimizer as first
                                and only argument. [default=None]
            | ``initial_step_size``  --  The initial step size used in the first
                                         call to the line search. For later
                                         line searches, the actual step size
                                         found by the previous line search is
                                         used as initial step size. How the
                                         initial step size is used, depends on
                                         the line search algorithm.
            | ``constraints``  --  An instance of the Constraints class.
            | ``debug_line``  --  If True, and when the line search fails, a
                                  plot with the line function will be made with
                                  matplotlib and written as
                                  ``'line_failed_%s.png' % isodatetime``.

           The function ``fun`` takes a mandatory argument ``x`` and an optional
           argument ``do_gradient``:
            | ``x``  --  the arguments of the function to be tested
            | ``do_gradient``  --  when False, only the function value is
                                   returned. when True, a 2-tuple with the
                                   function value and the gradient are returned
                                   [default=False]
        """
        if len(x_init.shape)!=1:
            raise ValueError("The unknowns must be stored in a plain row vector.")
        # self.x always contains the current parameters
        self.x = x_init.copy()
        if isinstance(fun, Preconditioner):
            self.prec = fun
        else:
            self.prec = None
        self.fun = FunWrapper(fun, anagrad, epsilon)
        self.line = LineWrapper(fun, anagrad, epsilon)
        self.search_direction = search_direction
        self.line_search = line_search
        self.convergence_condition = convergence_condition
        self.stop_loss_condition = stop_loss_condition
        self.epsilon = epsilon
        self.verbose = verbose
        self.callback = callback
        self.initial_step_size = initial_step_size
        self.constraints = constraints
        self.debug_line = debug_line

        # perform some resets:
        self.search_direction.reset()
        if self.stop_loss_condition is not None:
            self.stop_loss_condition.reset()

        # the current function value
        self.f = None
        # the current gradient
        self.gradient = None
        # the current step size
        self.step = None

        if self.convergence_condition is not None:
            self.success = self._run()

    def get_final(self):
        """Return the final solution in the original coordinates"""
        if self.prec is None:
            return self.x
        else:
            return self.prec.undo(self.x)

    def _run(self):
        """Run the iterative optimizer"""
        success = self.initialize()
        while success is None:
            success = self.propagate()
        return success

    def initialize(self):
        self.counter = 0
        if self.constraints is not None:
            try:
                self.x = self.constraints.free_shake(self.x)[0]
            except ConstraintError:
                self._screen("SHAKE FAILED", newline=True)
                return False
        self.f, self.gradient = self.fun(self.x, do_gradient=True)
        if self.constraints is not None:
            try:
                self.gradient = -self.constraints.project(self.x, -self.gradient)
            except ConstraintError:
                self._screen("CONSTRAINT PROJECT FAILED", newline=True)
                return False
        self.last_end = time.clock()

    def propagate(self):
        # compute the new direction
        self.search_direction.update(self.gradient, self.step)
        # print some stuff on screen
        if self.counter % 20 == 0:
            self._print_header()
        self._screen("% 5i   %2s" % (self.counter, self.search_direction.status))
        # perform a line search
        if self.constraints is not None:
            # Keep a copy of the last function value to make sure there is no
            # increase.
            self.last_f = self.f
        line_success = self._line_opt()
        if not line_success:
            self._screen("Line search failed", newline=True)
            if self.search_direction.is_sd():
                self._screen("LINE FAILED", newline=True)
                return False
        if self.constraints is not None:
            try:
                x, shake_count, active_count = self.constraints.free_shake(self.x)
            except ConstraintError:
                self._screen("SHAKE FAILED", newline=True)
                return False
            # test if the function did increase. if so try harder
            if self.fun(x) > self.last_f:
                try:
                    self.x, shake_count_bis, active_count = self.constraints.safe_shake(self.x, self.fun, self.last_f)
                except ConstraintError:
                    self._screen("SHAKE FAILED", newline=True)
                    return False
                shake_count += shake_count_bis
            else:
                self.x = x
            self._screen('%3i%3i' % (shake_count, active_count))
        # compute the gradient at the new point
        self.f, self.gradient = self.fun(self.x, do_gradient=True)
        self._screen("% 15.9e  " % self.f)
        if self.constraints is not None:
            try:
                self.gradient = -self.constraints.project(self.x, -self.gradient)
            except ConstraintError:
                self._screen("CONSTRAINT PROJECT FAILED", newline=True)
                return False
        # handle the preconditioner, part 1
        if self.prec is not None:
            gradient_orig = self.prec.do(self.gradient)
            step_orig = self.prec.undo(self.step)
        else:
            gradient_orig = self.gradient
            step_orig = self.step
        if self.convergence_condition is not None:
            # check convergence on the gradient and step in original basis
            converged, status = self.convergence_condition(gradient_orig, step_orig, self.f)
            self._screen(status)
        else:
            converged = False
        # timing
        end = time.clock()
        self._screen("%5.2f" % (end - self.last_end), newline=True)
        self.last_end = end
        # check convergence, part 2
        if converged:
            self._screen("CONVERGED", newline=True)
            return True
        # check stop loss on the gradient in original basis
        if self.stop_loss_condition is not None:
            lost = self.stop_loss_condition(self.counter, self.f, gradient_orig, step_orig)
            if lost:
                self._screen("LOST", newline=True)
                return False
        # call back
        if self.callback is not None:
            self.callback(self)
        # preconditioner, part 2
        if self.prec is not None:
           x_orig = self.prec.undo(self.x)
           if self.prec.update(self.counter, self.f, x_orig, gradient_orig):
                self.x = self.prec.do(x_orig)
                self.f, self.gradient = self.fun(self.x, do_gradient=True)
                self.step = None
                self.search_direction.reset()
        self.counter += 1

    def _print_header(self):
        """Print the header for screen logging"""
        header = " Iter  Dir  "
        if self.constraints is not None:
            header += '  SC CC'
        header += "         Function"
        if self.convergence_condition is not None:
            header += self.convergence_condition.get_header()
        header += "    Time"
        self._screen("-"*(len(header)), newline=True)
        self._screen(header, newline=True)
        self._screen("-"*(len(header)), newline=True)

    def _screen(self, s, newline=False):
        """Print something on screen when self.verbose == True"""
        if self.verbose:
            if newline:
                print(s)
            else:
                print(s, end=' ')

    def _line_opt(self):
        """Perform a line search along the current direction"""
        direction = self.search_direction.direction
        if self.constraints is not None:
            try:
                direction = self.constraints.project(self.x, direction)
            except ConstraintError:
                self._screen("CONSTRAINT PROJECT FAILED", newline=True)
                return False
        direction_norm = np.linalg.norm(direction)
        if direction_norm == 0:
            return False
        self.line.configure(self.x, direction/direction_norm)

        success, wolfe, qopt, fopt = \
            self.line_search(self.line, self.initial_step_size, self.epsilon)
        if success:
            self.step = qopt*self.line.axis
            self.initial_step_size = np.linalg.norm(self.step)
            self.x = self.x + self.step
            self.f = fopt
            if wolfe:
                self._screen("W")
            else:
                self._screen(" ")
                self.search_direction.reset()
            return True
        else:
            if self.debug_line:
                import matplotlib.pyplot as pt
                import datetime
                pt.clf()
                qs = np.arange(0.0, 100.1)*(5*self.initial_step_size/100.0)
                fs = np.array([self.line(q) for q in qs])
                pt.plot(qs, fs)
                pt.xlim(qs[0], qs[-1])
                fdelta = fs.max() - fs.min()
                if fdelta == 0.0:
                    fdelta = fs.mean()
                fmargin = fdelta*0.1
                pt.ylim(fs.min() - fmargin, fs.max() + fmargin)
                pt.title('fdelta = %.2e   fmean = %.2e' % (fdelta, fs.mean()))
                pt.xlabel('Line coordinate, q')
                pt.ylabel('Function value, f')
                pt.savefig('line_failed_%s.png' % (datetime.datetime.now().isoformat()))
            self._reset_state()
            return False

    def _reset_state(self):
        """Reset of the internal state of the function"""
        self.fun(self.x) # reset the internal state of the function


def check_anagrad(fun, x0, epsilon, threshold):
    """Check the analytical gradient using finite differences

       Arguments:
        | ``fun``  --  the function to be tested, more info below
        | ``x0``  --  the reference point around which the function should be
                      tested
        | ``epsilon``  --  a small scalar used for the finite differences
        | ``threshold``  --  the maximum acceptable difference between the
                             analytical gradient and the gradient obtained by
                             finite differentiation

       The function ``fun`` takes a mandatory argument ``x`` and an optional
       argument ``do_gradient``:
        | ``x``  --  the arguments of the function to be tested
        | ``do_gradient``  --  When False, only the function value is returned.
                               When True, a 2-tuple with the function value and
                               the gradient are returned [default=False]
    """
    N = len(x0)
    f0, ana_grad = fun(x0, do_gradient=True)
    for i in range(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        num_grad_comp = (fun(xh)-fun(xl))/epsilon
        if abs(num_grad_comp - ana_grad[i]) > threshold:
            raise AssertionError("Error in the analytical gradient, component %i, got %s, should be about %s" % (i, ana_grad[i], num_grad_comp))


def check_delta(fun, x, dxs, period=None):
    """Check the difference between two function values using the analytical gradient

       Arguments:
        | ``fun``  --  The function to be tested, more info below.
        | ``x``  --  The argument vector.
        | ``dxs``  --  A matrix where each row is a vector of small differences
                       to be added to the argument vector.

       Optional argument:
        | ``period``  --  If the function value is periodic, one may provide the
                          period such that differences are computed using
                          periodic boundary conditions.

       The function ``fun`` takes a mandatory argument ``x`` and an optional
       argument ``do_gradient``:
        | ``x``  --  The arguments of the function to be tested.
        | ``do_gradient``  --  When False, only the function value is returned.
                               When True, a 2-tuple with the function value and
                               the gradient are returned. [default=False]

       For every row in dxs, the following computation is repeated:

       1) D1 = 'f(x+dx) - f(x)' is computed.
       2) D2 = '0.5 (grad f(x+dx) + grad f(x)) . dx' is computed.

       A threshold is set to the median of the D1 set. For each case where |D1|
       is larger than the threshold, |D1 - D2|, should be smaller than the
       threshold.
    """
    dn1s = []
    dn2s = []
    dnds = []
    for dx in dxs:
        f0, grad0 = fun(x, do_gradient=True)
        f1, grad1 = fun(x+dx, do_gradient=True)
        grad = 0.5*(grad0+grad1)
        d1 = f1 - f0
        if period is not None:
            d1 -= np.floor(d1/period + 0.5)*period
        if hasattr(d1, '__iter__'):
            norm = np.linalg.norm
        else:
            norm = abs
        d2 = np.dot(grad, dx)

        dn1s.append(norm(d1))
        dn2s.append(norm(d2))
        dnds.append(norm(d1-d2))
    dn1s = np.array(dn1s)
    dn2s = np.array(dn2s)
    dnds = np.array(dnds)

    # Get the threshold (and mask)
    threshold = np.median(dn1s)
    mask = dn1s > threshold
    # Make sure that all cases for which dn1 is above the treshold, dnd is below
    # the threshold
    if not (dnds[mask] < threshold).all():
        raise AssertionError((
            'The first order approximation on the difference is too wrong. The '
            'threshold is %.1e.\n\nDifferences:\n%s\n\nFirst order '
            'approximation to differences:\n%s\n\nAbsolute errors:\n%s')
            % (threshold,
            ' '.join('%.1e' % v for v in dn1s[mask]),
            ' '.join('%.1e' % v for v in dn2s[mask]),
            ' '.join('%.1e' % v for v in dnds[mask])
        ))


def compute_fd_hessian(fun, x0, epsilon, anagrad=True):
    """Compute the Hessian using the finite difference method

       Arguments:
        | ``fun``  --  the function for which the Hessian should be computed,
                       more info below
        | ``x0``  --  the point at which the Hessian must be computed
        | ``epsilon``  --  a small scalar step size used to compute the finite
                           differences

       Optional argument:
        | ``anagrad``  --  when True, analytical gradients are used
                           [default=True]

       The function ``fun`` takes a mandatory argument ``x`` and an optional
       argument ``do_gradient``:
        | ``x``  --  the arguments of the function to be tested
        | ``do_gradient``  --  When False, only the function value is returned.
                               When True, a 2-tuple with the function value and
                               the gradient are returned [default=False]
    """
    N = len(x0)

    def compute_gradient(x):
        if anagrad:
            return fun(x, do_gradient=True)[1]
        else:
            gradient = np.zeros(N, float)
            for i in range(N):
                xh = x.copy()
                xh[i] += 0.5*epsilon
                xl = x.copy()
                xl[i] -= 0.5*epsilon
                gradient[i] = (fun(xh)-fun(xl))/epsilon
            return gradient

    hessian = np.zeros((N,N), float)
    for i in range(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        hessian[i] = (compute_gradient(xh) - compute_gradient(xl))/epsilon

    return 0.5*(hessian + hessian.transpose())
