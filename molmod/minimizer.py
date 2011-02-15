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


import numpy as np, time


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
            sy = np.dot(s, y)
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
    def __init__(self, c1=1e-4, c2=1e-1, max_iter=5):
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
        """
        self.c1 = c1
        self.c2 = c2
        self.max_iter = max_iter
        LineSearch.__init__(self)

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
        gl = fun(-0.5*epsilon, do_gradient=True)[1]
        gh = fun(+0.5*epsilon, do_gradient=True)[1]
        h0 = (gh-gl)/epsilon

        if h0 > 0:
            q1, f1, g1, h1 = 0.0, f0, g0, h0
            counter = 0
            wolfe = False
            while True:
                q2 = q1-g1/h1
                f2, g2 = fun(q2, do_gradient=True)
                if abs(g2) > abs(g1):
                    break
                counter += 1
                if self.max_iter is not None and counter > self.max_iter:
                    break
                q1, f1, g1 = q2, f2, g2
                del q2
                del f2
                del g2
                if f1 >= f0 + self.c1*abs(g0*q1):
                    break
                if f1 <= f0 - self.c1*abs(g0*q1) and abs(g1) <= abs(g0*self.c2):
                    wolfe = True
                    break
                gl = fun(q1-0.5*epsilon, do_gradient=True)[1]
                gh = fun(q1+0.5*epsilon, do_gradient=True)[1]
                h1 = (gh-gl)/epsilon
            if counter > 0 and f1 <= f0:
                return True, wolfe, q1, f1
            else:
                # even the first newton step failed, revert to back tracking
                pass

        # simple back tracking with tau = 0.5, no wolfe conditions yet
        q1 = -np.sign(g0)*initial_step_size*1.5
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
            for i in xrange(N):
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

    def __init__(self, step_rms=None, step_max=None, grad_rms=None, grad_max=None):
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

           Only the present arguments define when the minimization has
           converged. All actual values must go below the given thresholds.
        """
        if (step_rms is None and step_max is None and grad_rms is None and grad_max is None):
            raise ValueError("Some convergence criteria must be specified")
        self.step_rms = step_rms
        self.step_max = step_max
        self.grad_rms = grad_rms
        self.grad_max = grad_max

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
        return result

    def __call__(self, grad, step):
        """Return True when the minimizer has converged

           Arguments:
            | ``grad``  --  The gradient at the current point
            | ``step``  --  The last step vector
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

        stop, status = check_threshold(np.sqrt((step**2).mean()), self.step_rms, stop, status)
        stop, status = check_threshold(abs(step).max(), self.step_max, stop, status)
        stop, status = check_threshold(np.sqrt((grad**2).mean()), self.grad_rms, stop, status)
        stop, status = check_threshold(abs(grad).max(), self.grad_max, stop, status)

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
           Argument:
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
                for j in xrange(len(x)):
                    xh = x.copy()
                    xh[j] += 0.5*self.epsilon
                    xl = x.copy()
                    xl[j] -= 0.5*self.epsilon
                    g[j] = (self.fun(xh) - self.fun(xl))/self.epsilon
                return self.fun(x), g
        else:
            return self.fun(x)


class Constraints(object):
    def __init__(self, equations, threshold, rcond1=1e-10):
        self.equations = equations
        self.threshold = threshold
        self.rcond1 = rcond1

    def _compute_equations(self, x, indexes=None):
        # compute the error and the normals.
        normals = []
        values = []
        signs = []
        error = 0.0
        for i, (sign, equation) in enumerate(self.equations):
            value, normal = equation(x)
            if (indexes is not None and i in indexes) or (sign==-1 and value > -2*self.threshold) or (sign==0) or (sign==1 and value < 2*self.threshold):
                values.append(value)
                normals.append(normal)
                signs.append(sign)
                error += value**2
                if indexes is not None:
                    indexes.add(i)
        error = np.sqrt(error)
        normals = np.array(normals, float)
        values = np.array(values, float)
        return normals, values, error, signs

    def shake(self, x):
        indexes = set([])
        normals, values, error = self._compute_equations(x, indexes)[:-1]
        counter = 0
        while True:
            if error < self.threshold:
                break
            U, S, Vt = np.linalg.svd(normals, full_matrices=False)
            rcond = None
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
                dx = -np.dot(Vt.transpose(), np.dot(U.transpose(), values)*Sinv)
                new_x = x + dx
                new_indexes = set(indexes)
                new_normals, new_values, new_error = self._compute_equations(new_x, new_indexes)[:-1]
                new_gradient = np.dot(new_values, new_normals)
                #print counter, new_error - error, np.linalg.norm(new_gradient), np.dot(new_gradient, dx)
                if new_error < error:
                    counter += 1
                    x = new_x
                    index = new_indexes
                    normals = new_normals
                    values = new_values
                    error = new_error
                    break
                if np.dot(new_gradient, dx)/error < self.threshold:
                    raise RuntimeError('No feasible point found.')
        return x, counter, len(values)

    def project(self, x, vector):
        normals, signs = self._compute_equations(x)[::3]
        if len(normals) == 0:
            return vector
        U, S, Vt = np.linalg.svd(normals, full_matrices=False)
        decomposition = np.dot(Vt, vector)
        for i, sign in enumerate(signs):
            if sign == 0:
                continue
            if sign*decomposition[i] > 0:
                decomposition[i] = 0.0
        result = vector - np.dot(Vt.transpose(), decomposition)
        return result


class Minimizer(object):
    """A flexible multivariate minimizer

       The minimizer searches (in principle) for the 'nearest' local minimum.
    """

    def __init__(self, x_init, fun, search_direction, line_search,
                 convergence_condition, stop_loss_condition, anagrad=False,
                 epsilon=1e-6, verbose=True, callback=None,
                 initial_step_size=1.0, constraints=None):
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

        # perform some resets:
        self.search_direction.reset()
        self.stop_loss_condition.reset()

        # the current function value
        self.f = None
        # the current gradient
        self.gradient = None
        # the current step size
        self.step = None

        self.success = self._iterate()

    def get_final(self):
        """Return the final solution in the original coordinates"""
        if self.prec is None:
            return self.x
        else:
            return self.prec.undo(self.x)

    def _iterate(self):
        """Run the iterative optimizer"""
        if self.constraints is not None:
            self.x = self.constraints.shake(self.x)[0]
        self.f, self.gradient = self.fun(self.x, do_gradient=True)
        if self.constraints is not None:
            self.gradient = -self.constraints.project(self.x, -self.gradient)
        last_end = time.clock()
        # the cg loop
        self.counter = 0
        while True:
            # compute the new direction
            self.search_direction.update(self.gradient, self.step)
            # print some stuff on screen
            if self.counter % 20 == 0:
                self._print_header()
            self._screen("% 5i   %2s" % (self.counter, self.search_direction.status))
            # perform a line search
            line_success = self._line_opt()
            if not line_success:
                self._screen("Line search failed", newline=True)
                if self.search_direction.is_sd():
                    self._screen("LINE FAILED", newline=True)
                    return False
                else:
                    self.search_direction.reset()
                    continue
            if self.constraints is not None:
                self.x, shake_count, active_count = self.constraints.shake(self.x)
                self._screen('%3i%3i' % (shake_count, active_count))
            # compute the gradient at the new point
            self.f, self.gradient = self.fun(self.x, do_gradient=True)
            self._screen("% 15.9e  " % self.f)
            if self.constraints is not None:
                self.gradient = -self.constraints.project(self.x, -self.gradient)
            # handle the preconditioner, part 1
            if self.prec is not None:
                gradient_orig = self.prec.do(self.gradient)
                step_orig = self.prec.undo(self.step)
            else:
                gradient_orig = self.gradient
                step_orig = self.step
            # check convergence on the gradient and step in original basis
            converged, status = self.convergence_condition(gradient_orig, step_orig)
            self._screen(status)
            # timing
            end = time.clock()
            self._screen("%5.2f" % (end - last_end), newline=True)
            last_end = end
            # check convergence, part 2
            if converged:
                self._screen("CONVERGED", newline=True)
                return True
            # check stop loss on the gradient in original basis
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
        header = " Iter  Dir"
        if self.constraints is not None:
            header += '  SC CC'
        header += "         Function"
        header += self.convergence_condition.get_header()
        header += "    Time"
        self._screen("-"*(len(header)), newline=True)
        self._screen(header, newline=True)
        self._screen("-"*(len(header)), newline=True)

    def _screen(self, s, newline=False):
        """Print something on screen when self.verbose == True"""
        if self.verbose:
            if newline:
                print s
            else:
                print s,

    def _line_opt(self):
        """Perform a line search along the current direction"""
        direction = self.search_direction.direction
        if self.constraints is not None:
            direction = self.constraints.project(self.x, direction)
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
            if not wolfe:
                self.search_direction.reset()
            return True
        else:
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
    for i in xrange(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        num_grad_comp = (fun(xh)-fun(xl))/epsilon
        if abs(num_grad_comp - ana_grad[i]) > threshold:
            raise ValueError("Error in the analytical gradient, component %i, got %s, should be about %s" % (i, ana_grad[i], num_grad_comp))


def check_delta(fun, x, dxs, threshold):
    """Check the difference between two function values using the analytical gradient

       Arguments:
        | ``fun``  --  The function to be tested, more info below.
        | ``x``  --  The argument vector.
        | ``dxs``  --  A matrix where each row is a vector of small differences
                       to be added to the argument vector.
        | ``threshold``  --  The maximum acceptable deviation between the
                             difference of function values in x and x+dx and the
                             approximation of the two function values using a
                             first order Taylor approximation.

       The function ``fun`` takes a mandatory argument ``x`` and an optional
       argument ``do_gradient``:
        | ``x``  --  The arguments of the function to be tested.
        | ``do_gradient``  --  When False, only the function value is returned.
                               When True, a 2-tuple with the function value and
                               the gradient are returned. [default=False]

       For every row in dxs, the following test is repeated:

       1) D1 = 'f(x+dx) - f(x)' is computed.
       2) D2 = '(grad f(x+dx) + grad f(x)) . dx' is computed.

       There should be at least 50% of dx rows for which D1 is larger than the
       threshold. If not an AssertionError is raised.

       For each case where D1 is larger than the threshold, the absolute
       difference |D1 - D2| should be larger than the threshold.
    """
    df_small = []
    for dx in dxs:
        f0, grad0 = fun(x, do_gradient=True)
        f1, grad1 = fun(x+dx, do_gradient=True)
        grad = 0.5*(grad0+grad1)
        df = f1 - f0
        if abs(df) < threshold:
            df_small.append(df)
        else:
            expected = np.dot(dx, grad)
            if abs(df - expected) > threshold:
                raise AssertionError('Error in the analytical gradient: df is %s, should be about %s.' % (df, expected))
    if len(df_small)*2 > len(dxs):
        raise AssertionError('Less than 50%% of the dx rows leads to a sufficiently large D1. (%i out of %i). small df\'s: %s' % (
            len(dxs) - len(df_small), len(dxs), '[%s]' % ' '.join('%.1e' % df for df in df_small)
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
            for i in xrange(N):
                xh = x.copy()
                xh[i] += 0.5*epsilon
                xl = x.copy()
                xl[i] -= 0.5*epsilon
                gradient[i] = (fun(xh)-fun(xl))/epsilon
            return gradient

    hessian = np.zeros((N,N), float)
    for i in xrange(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        hessian[i] = (compute_gradient(xh) - compute_gradient(xl))/epsilon

    return 0.5*(hessian + hessian.transpose())
