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
    "GoldenLineSearch", "NewtonLineSearch", "NewtonGLineSearch", "Minimizer"
]


phi = 0.5*(1+numpy.sqrt(5))


class LineSearch(object):
    """Abstract base class for a line search implementation"""

    def __init__(self, qtol, max_step, max_iter):
        """Initialize a line search

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             max_step  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        self.qtol = qtol
        self.max_step = max_step
        self.max_iter = max_iter

    def limit_step(self, step):
        """Clip the a step within the maximum allowed range"""
        return numpy.clip(step, -self.max_step, self.max_step)

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

    def __init__(self, qtol, max_step, max_iter):
        """Initialize the golden section line search

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             max_step  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, max_step, max_iter)
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

    def __init__(self, qtol, max_step, max_iter):
        """Initialize the Newton line search (numerical gradients)

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             max_step  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, max_step, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return [("num_reduce", int)]

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        result = (self.num_reduce,)
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

    def __init__(self, qtol, max_step, max_iter):
        """Initialize the Newton line search (analytical gradients)

           Arguments:
             qtol  --  The threshold for displacements along the line. (When
                       displacements become smaller than qtol, we assume
                       convergence)
             max_step  --  The maximum step size of a line search
             max_iter  --  the maximum number of iteration for the line search
        """
        LineSearch.__init__(self, qtol, max_step, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        """Get a description of the extra fields from the line search for the log

           The format of the returned list is the same as a dtype of a record
           array in numpy.
        """
        return [("num_reduce", int)]

    def get_extra_log(self):
        """The values of the extra log fields from the line search"""
        result = (self.num_reduce,)
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
        fl, dl = fun(-epsilon, do_gradient=True)
        fh, dh = fun(+epsilon, do_gradient=True)

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


class Minimizer(object):
    """A flexible conjugate gradient optimizer with configurable line search

       The minimizer searches in principle for the 'nearest' local minimum.
    """

    def __init__(
        self, x_init, fun, LineSearchCls, ftol, xtol, max_step_rms, max_iter, max_line_iter,
        do_gradient=False, epsilon_init=1e-6, absftol=False, verbose=True, callback=None,
        min_iter=0, extra_log_dtypes=None
    ):
        """Initialize the minimizer

           Arguments:
             x_init  --  the initial guess for the minimum
             fun  --  function to be minimized (see below)
             LineSearchCls  --  The class of the line search implementation
             ftol  --  when changes in the function value drop below this
                       threshold, the algorithm stops
             xtol  --  when the norm of the changes in the unknowns drops below
                       this threshold, the algorithm stops
             max_step_rms  --  the maximum step size of one line search
             max_iter  --  the maximum number of CG iterations
             max_line_iter  --  the maximum number of iterations in the line
                                search
             do_gradient  --  use analytical gradients (default=False)
             epsilon_init  --  a small value compared to expected changes in the
                               unknowns (default=1e-6)
             absftol  --  consider the ftol threshold as an absolute change
                          instead of a relative change (default=False)
             verbose  --  print progress information on screen (default=True)
             callback  --  optional callback routine after each CG iteration.
                           the callback routine gets the current unknowns as
                           first and only argument. The callback routine must
                           return a (empty) list of additional log data.
                           (default=None)
             min_iter  --  the minimum number of CG iterations (default=0)
             extra_log_dtypes  --  dtype descriptions of the extra log fields
                                   returned by the callback routine
                                   (default=None)
        """
        if len(x_init.shape)!=1:
            raise ValueError("The unknowns must be stored in a plain row vector.")
        max_step = max_step_rms*numpy.sqrt(len(x_init))
        self.x = x_init.copy()
        self.fun = fun
        self.line_search = LineSearchCls(xtol, max_step, max_line_iter)
        self.ftol = ftol
        self.xtol = xtol
        self.max_step = max_step
        self.do_gradient = do_gradient
        self.callback = callback
        self.epsilon = epsilon_init
        self.epsilon_max = epsilon_init
        self.absftol = absftol
        self.verbose = verbose
        self.min_iter = min_iter

        self.extra_log_dtypes = self.line_search.get_extra_log_dtypes()
        if extra_log_dtypes is not None:
            self.extra_log_dtypes.extend(extra_log_dtypes)

        self.log = []
        self.step_size = 1.0
        self.step_rms = 1.0/numpy.sqrt(len(x_init))
        self.direction_cg = numpy.zeros(self.x.shape, float)
        self.direction_sd = numpy.zeros(self.x.shape, float)
        self.direction_sd_old = numpy.zeros(self.x.shape, float)

        self._iterate(max_iter)

    def _iterate(self, max_iter):
        """Run the conjugate gradient algorithm"""
        self.f = self.fun(self.x)

        self.beta = 0
        self._update_direction_sd()
        self.direction_cg[:] = self.direction_sd
        last_reset = True
        self.lower = False # do have one iteration were the function lowered?
        # the cg loop
        for self.counter in xrange(max_iter):
            self._screen("Iter % 5i of % 5i" % (self.counter, max_iter), False)
            if last_reset:
                lower = self._line_opt("SD")
            else:
                lower = self._line_opt("CG")
            self._append_log()
            if lower is True:
                self._update_cg()
                last_reset = False
                self.lower = True
            else:
                if last_reset and (self.counter > self.min_iter or lower is None):
                    break
                self._update_sd()
                last_reset = True
            if self.step_rms > 0:
                self.epsilon = min(self.step_rms*1e-2, self.epsilon_max)

        self._screen("Done")

    def _screen(self, s, newline=True):
        """Print something on screen when self.verbose == True"""
        if self.verbose:
            if newline:
                print s
            else:
                print s,

    def _update_direction_sd(self):
        """Update the gradient"""
        if self.do_gradient:
            self.f, tmp = self.fun(self.x, do_gradient=True)
            self.direction_sd[:] = -tmp
            #print self.direction_sd
        else:
        #if True:
            tmp = self.x.copy()
            for j in xrange(len(self.x)):
                tmp[j] += self.epsilon
                self.direction_sd[j] = -(self.fun(tmp) - self.f)/self.epsilon
                tmp[j] = self.x[j]
            self._reset_state()
            #print self.direction_sd
        #sys.exit()

    def _line_opt(self, label):
        """Perform a line search along the given direction"""
        direction = self.direction_cg
        if numpy.linalg.norm(direction) == 0:
            self._screen("  'Opt% 3s'              Direction zero" % label)
            return False
        direction = direction/numpy.linalg.norm(direction)
        def fun_aux(q, do_gradient=False):
            xq = self.x + q*direction
            if do_gradient:
                fq, gq = self.fun(xq, do_gradient=True)
                return fq, numpy.dot(direction, gq)
            else:
                return self.fun(xq)

        success, qopt, fopt = self.line_search(fun_aux, self.f, self.step_size, self.epsilon)
        if success:
            self._screen("  'Opt% 3s'             " % label, False)
            return self._handle_new_solution(self.x + qopt*direction, fopt)
        else:
            # returning the old x and f, means that the optimization failed
            self._screen("  'Opt% 3s'              Line search failed" % label)
            self._reset_state()
            return None

    def _handle_new_solution(self, xnew, fnew):
        """Decide what to do with a new set of unknowns generated by the line search"""
        self.step_size = numpy.linalg.norm(self.x - xnew)
        self.step_rms = self.step_size/numpy.sqrt(len(self.x))
        if self.absftol:
            self.decrease = self.f - fnew
            self._screen("step_rms=% 9.7e   absdecr=% 9.7e   fnew=% 9.7e" % (self.step_rms, self.decrease, fnew))
        else:
            self.decrease = (self.f - fnew)/self.f
            self._screen("step_rms=% 9.7e   reldecr=% 9.7e   fnew=% 9.7e" % (self.step_rms, self.decrease, fnew))
        if self.step_rms < self.xtol:
            result = False
        if self.decrease < self.ftol:
            result = False
        else:
            result = True
        if self.decrease > 0:
            self.x = xnew
            self.f = fnew
        else:
            self._reset_state()
        return result

    def _reset_state(self):
        """Reset of the internal state of the function"""
        self.fun(self.x) # reset the internal state of the function
        self.decrease = -1.0

    def _update_cg(self):
        """Update the conjugate gradient after an iteration"""
        self.direction_sd_old[:] = self.direction_sd
        self._update_direction_sd()
        self.beta = (
            numpy.dot(self.direction_sd, self.direction_sd - self.direction_sd_old) /
            numpy.dot(self.direction_sd_old, self.direction_sd_old)
        )
        if self.beta < 0:
            self.beta = 0
        self.direction_cg[:] *= self.beta
        self.direction_cg[:] += self.direction_sd

    def _update_sd(self):
        """Reset the conjugate gradient to the normal gradient"""
        self._update_direction_sd()
        self.direction_cg[:] = self.direction_sd
        self.beta = 0

    def _append_log(self):
        """Collect various data into the log"""
        extra_fields = self.line_search.get_extra_log()
        if self.callback is not None:
            extra_fields = extra_fields + tuple(self.callback(self.x))
        self.log.append((
            self.f, self.step_rms, self.beta,
            numpy.linalg.norm(self.direction_sd),
            numpy.linalg.norm(self.direction_cg),
            self.x.copy(),
        ) + extra_fields)

    def get_log(self):
        """Return the log"""
        return numpy.array(self.log, [
            ("f", float), ("step_rms", float), ("beta", float),
            ("norm_sd", float), ("norm_cg", float), ("x", float, self.x.shape),
        ] + self.extra_log_dtypes)


