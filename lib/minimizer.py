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
    "GoldenLineSearch", "NewtonLineSearch", "NewtonGLineSearch",
    "StopCondition", "Minimizer"
]


phi = 0.5*(1+numpy.sqrt(5))


class LineSearch(object):
    """Abstract base class for a line search implementation"""

    def __init__(self, qtol, qmax, max_iter):
        """Initialize a line search

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             qmax  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        self.qtol = qtol
        self.qmax = qmax
        self.max_iter = max_iter

    def limit_step(self, step):
        """Clip the a step within the maximum allowed range"""
        return numpy.clip(step, -self.qmax, self.qmax)

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return []

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        return ()

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
    """An implementation of the golden section line search algorithm"""

    def __init__(self, qtol, qmax, max_iter):
        """Initialize the golden section line search

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             qmax  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, qmax, max_iter)
        self.num_bracket = 0
        self.num_golden = 0

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return [("num_bracket", int), ("num_golden", int)]

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        result = (self.num_bracket, self.num_golden)
        self.num_bracket = 0
        self.num_golden = 0
        return result

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
    """An implementation of the Newton line search algorithm with numerical gradients"""

    def __init__(self, qtol, qmax, max_iter):
        """Initialize the Newton line search (numerical gradients)

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             qmax  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, qmax, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return [("num_reduce", int)]

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        result = (self.num_reduce, )
        self.num_reduce = 0
        return result

    def __call__(self, fun, f0, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size
        """
        # determine the curvature:
        fl = fun(-epsilon)
        fh = fun(+epsilon)

        f_d1 = (fh-fl)/(2*epsilon)
        f_d2 = (fh+fl-2*f0)/epsilon
        if f_d2 > 0:
            qopt = -f_d1/f_d2
        else:
            qopt = last_step_size*1.5
        qopt = self.limit_step(qopt)
        self.num_reduce = 0
        while True:
            if qopt < self.qtol:
                return False, 0.0, f0
            fopt = fun(qopt)
            if fopt < f0:
                break
            qopt *= 0.5
            self.num_reduce += 1
            if self.num_reduce > self.max_iter:
                raise Exception("Line search did not converge.")
        return True, qopt, fopt


class NewtonGLineSearch(LineSearch):
    """An implementation of the Newton line search algorithm with analytical gradients"""

    def __init__(self, qtol, qmax, max_iter):
        """Initialize the Newton line search (analytical gradients)

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             qmax  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, qmax, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return [("num_reduce", int)]

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        result = (self.num_reduce, )
        self.num_reduce = 0
        return result

    def __call__(self, fun, f0, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size
        """
        # determine the curvature:
        dl = fun(-epsilon, do_gradient=True)[1]
        dh = fun(+epsilon, do_gradient=True)[1]

        f_d1 = 0.5*(dl+dh)
        f_d2 = (dh-dl)/(2*epsilon)
        if f_d2 > 0:
            qopt = -f_d1/f_d2
        else:
            qopt = last_step_size*1.5
        qopt = self.limit_step(qopt)
        self.num_reduce = 0
        while True:
            if qopt < self.qtol:
                return False, 0.0, f0
            fopt = fun(qopt)
            if fopt < f0:
                break
            qopt *= 0.5
            self.num_reduce += 1
            if self.num_reduce > self.max_iter:
                raise Exception("Line search did not converge.")
        return True, qopt, fopt


class StopCondition(object):
    def __init__(self, step_rms=None, step_max=None, grad_rms=None, grad_max=None, max_iter=None):
        self.step_rms = step_rms
        self.step_max = step_max
        self.grad_rms = grad_rms
        self.grad_max = grad_max
        self.max_iter = max_iter

    def get_header(self):
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

    def __call__(self, grad, step, counter):
        """Determine if the minimizer has converged"""
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

        if counter >= self.max_iter:
            stop = True
        return stop, status


class Minimizer(object):
    """A flexible conjugate gradient optimizer with configurable line search

       The minimizer searches in principle for the 'nearest' local minimum.
    """

    def __init__(self, x_init, fun, line_search, stop_condition, anagrad=False,
                 epsilon_init=1e-6, verbose=True, callback=None,
                 extra_log_dtypes=[]):
        """Initialize the minimizer

           Arguments:
             x_init  --  the initial guess for the minimum
             fun  --  function to be minimized (see below)
             line_search  --  a LineSearch object
             stop_condition  --  a StopCondition object

          Optional arguments
             anagrad  --  when set to True, analytical gradients are used
             epsilon_init  --  a small value compared to expected changes in the
                               unknowns (default=1e-6). it is used to compute
                               finite differences. epsilon will be reduced when
                               the steps in the iterative procedure become
                               smaller
             verbose  --  print progress information on screen (default=True)
             callback  --  optional callback routine after each CG iteration.
                           the callback routine gets the current unknowns as
                           first and only argument. The callback routine must
                           return a list of additional log data. (default=None)
             extra_log_dtypes  --  dtype descriptions of the extra log fields
                                   returned by the callback routine
                                   (default=[])
        """
        if len(x_init.shape)!=1:
            raise ValueError("The unknowns must be stored in a plain row vector.")
        self.x = x_init.copy()
        self.fun = fun
        self.line_search = line_search
        self.stop_condition = stop_condition
        self.anagrad = anagrad
        self.callback = callback
        self.epsilon = epsilon_init
        self.epsilon_max = epsilon_init
        self.verbose = verbose

        self.extra_log_dtypes = self.line_search.get_extra_log_dtypes()
        self.extra_log_dtypes.extend(extra_log_dtypes)

        self.log = []
        self.step_size = 1.0
        self.step_rms = 1.0/numpy.sqrt(len(x_init))
        # minus the (conjugate) gradient
        self.direction = numpy.zeros(self.x.shape, float)
        # minus the gradient
        self.direction_sd = numpy.zeros(self.x.shape, float)
        # previous value of minus the gradient
        self.direction_sd_old = numpy.zeros(self.x.shape, float)

        self._iterate()

    def _iterate(self):
        """Run the conjugate gradient algorithm"""
        self.f = self.fun(self.x)
        self._update_sd()
        last_reset = True
        # the cg loop
        self.counter = 0
        while True:
            if self.counter % 20 == 0:
                self._print_header()
            self._screen("% 5i  " % self.counter)
            if last_reset:
                self._line_opt("SD")
            else:
                self._line_opt("CG")
            self._append_log()
            # decide what to do in the next loop
            if self.decrease > 0:
                # check stop condition
                stop, status = self.stop_condition(-self.direction_sd, self.step, self.counter)
                self._screen(status, newline=True)
                if stop:
                    # Stopping because the stop condition was reached
                    break
                self._update_cg()
                last_reset = False
            else:
                self._screen("", newline=True)
                if last_reset or self.direction_norm == 0:
                    # Stopping because the conjugate gradient in combination
                    # with the line search is not finding a better solution,
                    # even after a plain steepest descent step. Most of the time
                    # this due to:
                    #  (i) ill-defined cost functions
                    #  (ii) errors in the gradient inherent to finite
                    #       differences
                    #  (iii) errors in the implementation of the analytical
                    #        gradient provided by the user
                    #  (iv) convergence criteria that go beyond double precision
                    # Note that the last issue should be seen in the relative
                    # sense. If the function value itself is about 10e10, it is
                    # not possible to see the function decrease with amounts
                    # smaller than 10e-5.
                    break
                self._update_sd()
                last_reset = True
            # update epsilon for finite differences
            if self.step_rms > 0:
                self.epsilon = min(self.step_rms*1e-2, self.epsilon_max)
            self.counter += 1

        self._screen("Done", True)

    def _print_header(self):
        header = " Iter  Dir     Decrease    Function"
        header += self.stop_condition.get_header()
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

    def _line_opt(self, label):
        """Perform a line search along the given direction"""
        self._screen("%2s  " % label)
        self.direction_norm = numpy.linalg.norm(self.direction)
        if self.direction_norm == 0:
            self._screen("Direction zero")
            self.decrease = 0.0
            return
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
            return self._handle_new_solution(self.x + qopt*axis, fopt)
        else:
            # returning the old x and f, means that the optimization failed
            self._screen("Line search failed")
            self._reset_state()
            self.decrease = 0.0

    def _handle_new_solution(self, xnew, fnew):
        """Analyze the new solution generated by the line search

           Returns True when the function lowered, false otherwise
        """
        self.step = self.x - xnew
        self.step_size = numpy.linalg.norm(self.step)
        self.step_rms = self.step_size/numpy.sqrt(len(self.x))
        self.decrease = self.f - fnew
        self._screen("% 9.3e  % 9.3e  " % (self.decrease, fnew))
        if self.decrease > 0:
            self.x = xnew
            self.f = fnew
            return True
        else:
            self._screen("Function increased")
            self._reset_state()
            return False

    def _reset_state(self):
        """Reset of the internal state of the function"""
        self.fun(self.x) # reset the internal state of the function

    def _update_cg(self):
        """Update the conjugate gradient after an iteration"""
        self.direction_sd_old[:] = self.direction_sd
        self._compute_direction_sd()
        # Polak-Ribiere
        self.beta = (
            numpy.dot(self.direction_sd, self.direction_sd - self.direction_sd_old) /
            numpy.dot(self.direction_sd_old, self.direction_sd_old)
        )
        # Automatic direction reset
        if self.beta < 0:
            self.beta = 0
        self.direction[:] *= self.beta
        self.direction[:] += self.direction_sd

    def _update_sd(self):
        """Reset the conjugate gradient to the normal gradient"""
        self._compute_direction_sd()
        self.direction[:] = self.direction_sd
        self.beta = 0

    def _append_log(self):
        """Collect various data into the log"""
        extra_fields = self.line_search.get_extra_log()
        if self.callback is not None:
            extra_fields = extra_fields + tuple(self.callback(self.x))
        self.log.append((
            self.f, self.step_rms, self.beta,
            numpy.linalg.norm(self.direction_sd),
            numpy.linalg.norm(self.direction),
            self.x.copy(),
        ) + extra_fields)

    def get_log(self):
        """Return the log"""
        return numpy.array(self.log, [
            ("f", float), ("step_rms", float), ("beta", float),
            ("norm_sd", float), ("norm_cg", float), ("x", float, self.x.shape),
        ] + self.extra_log_dtypes)


