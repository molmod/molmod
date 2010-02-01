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
"""General purpose minimization of continuous or smooth multidimensional functions

   The implementation is mainly concerned with robustness, rather than
   computational efficiency. Example usage:

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

   x_init = numpy.zeros(2, float)
   minimizer = Minimizer(x_init, fun, NewtonGLineSearch, 1e-5, 1e-5, 1e-1, 1000, 50, do_gradient=True)
   print "optimum", minimizer.x, fun(minimizer.x)
"""

import numpy


__all__ = [
    "GoldenLineSearch", "NewtonLineSearch",
    "ConvergenceCondition", "StopLossCondition", "Minimizer"
]


phi = 0.5*(1+numpy.sqrt(5))


class LineSearch(object):
    """Abstract base class for a line search"""

    def __init__(self, qmax=None):
        """Initialize a line search

           Optional arguments:
             qmax  --  The maximum step size of a line search
        """
        self.qmax = qmax

    def limit_step(self, step):
        """Clip the a step within the maximum allowed range"""
        if self.qmax is None:
            return step
        else:
            return numpy.clip(step, -self.qmax, self.qmax)

    def __call__(self, fun, f0, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size
        """
        raise NotImplementedError


class GoldenLineSearch(LineSearch):
    """The golden section line search algorithm"""

    def __init__(self, qtol, qmax=None, max_iter=None):
        """Initialize the golden section line search

           Argument:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
           Optional arguments
             qmax  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
                           (only applies to the bracketing part)
        """
        if qtol is None:
            raise ValueError("No stop condition is specified")
        self.qtol = qtol
        self.max_iter = max_iter
        LineSearch.__init__(self, qmax)
        self.num_bracket = 0
        self.num_golden = 0

    def __call__(self, fun, f0, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size
        """
        # bracket the minimum
        triplet = self._bracket(last_step_size, f0, fun)
        if triplet is None:
            return False, 0.0, f0
        # do a golden section optimization
        qopt, fopt = self._golden(triplet, fun)
        qopt = self.limit_step(qopt)
        fopt = fun(qopt)
        return True, qopt, fopt

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
                if counter > self.max_iter:
                    raise Exception("Line search did not converge.")
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
                if counter > self.max_iter:
                    raise Exception("Line search did not converge.")

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
    """The Newton line search algorithm with numerical gradients

       Only a single Newton or steepest descent step is performed.

       When the curvature is negative, a steepest descent step is tried, using
       the step size from the previous step multiplied by 1.5. If the new
       function value is higher, the step size is reduced by a factor two. The
       latter is repeated at most max_reduce times. If no lower value is found,
       the line search fails.

       When the curvature is positive, a Newton step is performed. When the
       gradient at the new point is not lower, the step size is reduced by a
       factor two. The latter is repeated at most max_reduce times. If not lower
       gradient is found, the line search fails.
    """

    def __init__(self, qmax=None, max_reduce=None, anagrad=True):
        """Initialize the Newton line search (numerical gradients)

           Optional arguments:
             qmax  --  The maximum step size of a line search
             max_reduce  --  the maximum number of attempts to half the step
                             size in order to get a lower function or gradient
             anagrad  --  when True, use analytical gradients [default=False]
        """
        self.max_reduce = max_reduce
        self.anagrad = anagrad
        LineSearch.__init__(self, qmax)
        self.num_reduce = 0

    def _compute_derivatives(self, q0, f0, fun, epsilon):
        if self.anagrad:
            fl = fun(q0-epsilon)
            fh = fun(q0+epsilon)
            g0 = (fh-fl)/(2*epsilon)
            h0 = (fh+fl-2*f0)/epsilon
        else:
            gl = fun(q0-epsilon, do_gradient=True)[1]
            gh = fun(q0+epsilon, do_gradient=True)[1]
            g0 = 0.5*(gl+gh)
            h0 = (gh-gl)/(2*epsilon)
        return g0, h0

    def __call__(self, fun, f0, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size, used
                          for finite differences
        """
        g0, h0 = self._compute_derivatives(0.0, f0, fun, epsilon)

        # determine the initial step size
        if h0 > 0:
            qopt = -g0/h0
        else:
            qopt = numpy.sign(-g0)*last_step_size*1.5
        qopt = self.limit_step(qopt)
        self.num_reduce = 0
        while True:
            fopt = fun(qopt)
            if h0 > 0:
                gopt, hopt = self._compute_derivatives(qopt, fopt, fun, epsilon)
                if abs(gopt) < abs(g0):
                    return True, qopt, fopt
            else:
                fopt = fun(qopt)
                if fopt < f0:
                    return True, qopt, fopt
            qopt *= 0.5
            self.num_reduce += 1
            if self.max_reduce is not None and self.num_reduce > self.max_reduce:
                return False, qopt, fopt



class ConvergenceCondition(object):
    """Callable object that tests the convergence of the minimizer
    """
    def __init__(self, step_rms=None, step_max=None, grad_rms=None, grad_max=None):
        """Initialize a ConvergenceCondition object

           Optional arguments:
             step_rms  --  threshold for the RMS value of the step vector in the
                           iterative minimizer
             step_max  --  threshold for the maximum component of the step
                           vector in the iterative minimizer
             grad_rms  --  threshold for the RMS value of the gradient
                           components of the function to be minimized
             grad_max  --  threshold for the maximum value of the gradient
                           components of the function to be minimized

           The present arguments define when the minimization has converged.
           All actual values must go below the given thresholds.
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
        """Return True when the minimizer has converged"""
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

        stop, status = check_threshold(numpy.sqrt((step**2).mean()), self.step_rms, stop, status)
        stop, status = check_threshold(abs(step).max(), self.step_max, stop, status)
        stop, status = check_threshold(numpy.sqrt((grad**2).mean()), self.grad_rms, stop, status)
        stop, status = check_threshold(abs(grad).max(), self.grad_max, stop, status)

        return stop, status


class StopLossCondition(object):
    """Callable object that checks if minimizer has lost track
    """
    def __init__(self, max_iter=None, fn_margin=None, grad_margin=None):
        """Initialize a StopLossCondition object

           Optional arguments:
             max_iter  --  the maximum number of iterations allowed
             fn_margin  --  if the function to be minimized goes above the
                            lowest value so far plus this margin, the
                            minimization is aborted
             grad_margin  --  if the RMS value of the gradient components goes
                              above the lowest value plus this threshold, the
                              minimization is aborted

           The present arguments define when the minimization has lost track.
        """
        self.max_iter = max_iter
        self.fn_margin = fn_margin
        self.grad_margin = grad_margin

        self.fn_lowest = None
        self.grad_rms_lowest = None

    def __call__(self, counter, fn, gradient):
        """Return True when the minimizer has lost track"""
        if self.max_iter is not None and counter >= self.max_iter:
            return True

        if self.fn_margin is not None:
            if self.fn_lowest is None or fn < self.fn_lowest:
                self.fn_lowest = fn
            elif fn > self.fn_lowest + self.fn_margin:
                return True

        if self.grad_margin is not None:
            grad_rms = numpy.sqrt((gradient**2).mean())
            if self.grad_rms_lowest is None or grad_rms < self.grad_rms_lowest:
                self.grad_rms_lowest = grad_rms
            elif grad_rms > self.grad_rms_lowest + self.grad_margin:
                return True

        # all is fine
        return False

class Minimizer(object):
    """A flexible conjugate gradient minimizer

       The minimizer searches in principle for the 'nearest' local minimum.
    """

    def __init__(self, x_init, fun, line_search, convergence_condition,
                 stop_loss_condition, anagrad=False, epsilon_init=1e-6,
                 verbose=True, callback=None):
        """Initialize the minimizer

           Arguments:
             x_init  --  the initial guess for the minimum
             fun  --  function to be minimized (see below)
             line_search  --  a LineSearch object
             convergence_condition  --  a ConvergenceCondition object
             stop_loss_condition  --  a StopLossCondition object

          Optional arguments
             anagrad  --  when set to True, analytical gradients are used
             epsilon_init  --  a small value compared to expected changes in the
                               unknowns (default=1e-6). it is used to compute
                               finite differences. epsilon will be reduced when
                               the steps in the iterative procedure become
                               smaller
             verbose  --  print progress information on screen (default=True)
             callback  --  optional callback routine after each CG iteration.
                           the callback routine gets the minimizer as first
                           and only argument. (default=None)
        """
        if len(x_init.shape)!=1:
            raise ValueError("The unknowns must be stored in a plain row vector.")
        self.x = x_init.copy()
        self.fun = fun
        self.line_search = line_search
        self.convergence_condition = convergence_condition
        self.stop_loss_condition = stop_loss_condition
        self.anagrad = anagrad
        self.callback = callback
        self.epsilon = epsilon_init
        self.epsilon_max = epsilon_init
        self.verbose = verbose

        self.step_size = 1.0
        self.step_rms = 1.0/numpy.sqrt(len(x_init))
        # minus the (conjugate) gradient
        self.direction = numpy.zeros(self.x.shape, float)
        # minus the gradient
        self.direction_sd = numpy.zeros(self.x.shape, float)
        # previous value of minus the gradient
        self.direction_sd_old = numpy.zeros(self.x.shape, float)
        self.direction_label = "??"

        self._iterate()

    def _iterate(self):
        """Run the conjugate gradient algorithm"""
        self.f = self.fun(self.x)
        self._compute_direction_sd()
        self._update_sd()
        # the cg loop
        self.counter = 0
        while True:
            if self.counter % 20 == 0:
                self._print_header()
            self._screen("% 5i  " % self.counter)
            # perform a line search
            line_success = self._line_opt()
            if not line_success:
                self._screen("Line search failed", newline=True)
                if self.direction_label == "SD":
                    self._screen("LINE FAILED", newline=True)
                    return
                else:
                    self._update_sd()
                    continue
            # compute the direction (or -gradient) at the new point
            self._compute_direction_sd()
            # check convergence
            converged, status = self.convergence_condition(-self.direction_sd, self.step)
            self._screen(status, newline=True)
            if converged:
                self._screen("CONVERGED", newline=True)
                return
            # check stop loss
            lost = self.stop_loss_condition(self.counter, self.f, -self.direction_sd)
            if lost:
                self._screen("LOST", newline=True)
                return

            # compute the new direction
            self._update_cg()
            # update epsilon for finite differences
            if self.step_rms > 0:
                self.epsilon = min(self.step_rms*1e-2, self.epsilon_max)
            self.counter += 1
            # call back
            if self.callback is not None:
                self.callback(self)

    def _print_header(self):
        header = " Iter  Dir     Function"
        header += self.convergence_condition.get_header()
        self._screen("-"*(len(header)+2), True)
        self._screen(header, True)
        self._screen("-"*(len(header)+2), True)

    def _screen(self, s, newline=False):
        """Print something on screen when self.verbose == True"""
        if self.verbose:
            if newline:
                print s
            else:
                print s,

    def _compute_direction_sd(self):
        """Update the gradient"""
        self.direction_sd_old[:] = self.direction_sd
        if self.anagrad:
            self.f, tmp = self.fun(self.x, do_gradient=True)
            self.direction_sd[:] = -tmp
        else:
            tmp = self.x.copy()
            for j in xrange(len(self.x)):
                tmp[j] += self.epsilon
                self.direction_sd[j] = -(self.fun(tmp) - self.f)/self.epsilon
                tmp[j] = self.x[j]
            self._reset_state()

    def _line_opt(self):
        """Perform a line search along the given direction"""
        self._screen("%2s  " % self.direction_label)
        self.direction_norm = numpy.linalg.norm(self.direction)
        if self.direction_norm == 0:
            return False
        axis = self.direction/self.direction_norm

        def fun_aux(q, do_gradient=False):
            """One-dimensional cut of the function along the search direction"""
            xq = self.x + q*axis
            if do_gradient:
                fq, gq = self.fun(xq, do_gradient=True)
                return fq, numpy.dot(axis, gq)
            else:
                return self.fun(xq)

        success, qopt, fopt = self.line_search(fun_aux, self.f, self.step_size, self.epsilon)
        if success:
            self.step = qopt*axis
            self.x = self.x + self.step
            self.f = fopt
            self._screen("% 9.3e  " % self.f)
            self.step_size = numpy.linalg.norm(self.step)
            self.step_rms = self.step_size/numpy.sqrt(len(self.x))
            return True
        else:
            self._reset_state()
            return False

    def _reset_state(self):
        """Reset of the internal state of the function"""
        self.fun(self.x) # reset the internal state of the function

    def _update_cg(self):
        """Update the conjugate gradient after an iteration"""
        # Polak-Ribiere
        self.beta = (
            numpy.dot(self.direction_sd, self.direction_sd - self.direction_sd_old) /
            numpy.dot(self.direction_sd_old, self.direction_sd_old)
        )
        # Automatic direction reset
        if self.beta < 0:
            self.beta = 0
            self.direction_label = "SD"
        else:
            self.direction_label = "CG"
        self.direction[:] *= self.beta
        self.direction[:] += self.direction_sd

    def _update_sd(self):
        """Reset the conjugate gradient to the normal gradient"""
        self.direction[:] = self.direction_sd
        self.beta = 0
        self.direction_label = "SD"

