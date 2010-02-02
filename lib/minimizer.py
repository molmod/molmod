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
   search_direction = ConjugateGradient()
   line_search = NewtonLineSearch()
   convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6)
   stop_loss = StopLossCondition(max_iter=50)
   minimizer = Minimizer(
       x_init, fun, search_direction, line_search, convergence, stop_loss,
       anagrad=True, verbose=True,
   )
   print "optimum", minimizer.x, fun(minimizer.x)
"""


import numpy, time


__all__ = [
    "SearchDirection", "SteepestDescent", "ConjugateGradient", "QuasiNewton",
    "LineSearch", "GoldenLineSearch", "NewtonLineSearch",
    "ConvergenceCondition", "StopLossCondition", "Minimizer",
    "check_anagrad", "compute_fd_hessian",
]


class SearchDirection(object):
    """Abstract base class for a search direction method"""
    def __init__(self):
        """Initialize a SearchDirection object"""
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
        """Initialize a SearchDirection object"""
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
    def __init__(self, variant='Polak-Ribiere'):
        """Initialize a ConjugateGradient object

           Optional argument:
             variant  --  the variant of the Conjugate Gradient method to use:
                          'Fletcher-Reeves' or 'Polak-Ribiere' (default)
        """
        if variant == 'Polak-Ribiere':
            self._beta = self._beta_pr
        elif variant == 'Fletcher-Reeves':
            self._beta = self._beta_fr
        else:
            raise ValueError("Unknown conjugate gradient variant: %s" % variant)
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

    def _beta_pr(self):
        # Polak-Ribiere
        return (
            numpy.dot(self.gradient, self.gradient - self.gradient_old) /
            numpy.dot(self.gradient_old, self.gradient_old)
        )

    def _beta_fr(self):
        # Fletcher-Reeves
        return (
            numpy.dot(self.gradient, self.gradient) /
            numpy.dot(self.gradient_old, self.gradient_old)
        )

    def _update_sd(self):
        """Reset the conjugate gradient to the normal gradient"""
        self.direction = -self.gradient
        self.status = "SD"


class QuasiNewton(SearchDirection):
    """The quasi Newton method

       bla bla bla
    """
    def __init__(self):
        """Initialize a QuasiNewton object"""
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
            self.inv_hessian = numpy.identity(N, float)
        else:
            # update the direction
            self.direction = -numpy.dot(self.inv_hessian, self.gradient)
            self.status = "QN"
            # new guess of the inverse hessian (BFGS)
            y = self.gradient - self.old_gradient
            s = step
            sy = numpy.dot(s, y)
            A = numpy.outer(-y/sy, s)
            A.ravel()[::N+1] += 1
            self.inv_hessian = (
                numpy.dot(numpy.dot(A.transpose(), self.inv_hessian), A) +
                numpy.outer(s/sy, s)
            )

    def reset(self):
        """Reset the internal state of the search direction algorithm"""
        self.inv_hessian = None
        self.gradient = None
        self.old_gradient = None


    def is_sd(self):
        """Return True if the last direction was steepest descent"""
        return self.status == "SD"


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

           Returns:
             success  --  a boolean indicating that the line search resulted
                          in an improved solution
             wolfe  --  a boolean indicating that the new solution satisfies
                        wolfe conditions
             qopt  --  the position of the new solution on the line
             fopt  --  the corresponding function value
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

    def __call__(self, fun, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size

           Returns:
             success  --  a boolean indicating that the line search resulted
                          in an improved solution
             wolfe  --  a boolean indicating that the new solution satisfies
                        wolfe conditions
             qopt  --  the position of the new solution on the line
             fopt  --  the corresponding function value

           P.S. The wolfe parameter is always True, but this aspect is not
           guaranteed to be correct. Never use the GoldenLineSearch in
           combination with a quasi Newton method.
        """
        # bracket the minimum
        f0 = fun(0.0)
        triplet = self._bracket(last_step_size, f0, fun)
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
                if counter > self.max_iter:
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
                if counter > self.max_iter:
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
        """Initialize the Newton line search

           Optional arguments:
             c1  --  The coefficient in the first Wolfe condition (sufficient
                     decrease of the function) [default=1e-4]
             c2  --  The coefficient in the second Wolfe condition (sufficient
                     decrease of the derivative) [default=1e-1]. the default
                     is optimal for the conjugate gradient method
             max_iter  --  the maximum number of iterations in the line search.
        """
        self.c1 = c1
        self.c2 = c2
        self.max_iter = max_iter
        LineSearch.__init__(self)

    def __call__(self, fun, last_step_size, epsilon):
        """Return the value that minimizes the one-dimensional function 'fun'

           Arguments:
             fun  --  function to minimize (one-dimensional)
             f0   --  the function value at the starting point q=0"
             last_step_size  --  the norm of step from the previous line search
             epsilon  --  a value that is small compared to last_step_size, used
                          for finite differences

           Returns:
             success  --  a boolean indicating that the line search resulted
                          in an improved solution
             wolfe  --  a boolean indicating that the new solution satisfies
                        wolfe conditions
             qopt  --  the position of the new solution on the line
             fopt  --  the corresponding function value
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
                if counter > self.max_iter:
                    break
                q1, f1, g1 = q2, f2, g2
                del q2
                del f2
                del g2
                if f1 <= f0 - self.c1*abs(g0*q1) and abs(g1) <= abs(g0*self.c2):
                    wolfe = True
                    break
                gl = fun(q1-0.5*epsilon, do_gradient=True)[1]
                gh = fun(q1+0.5*epsilon, do_gradient=True)[1]
                h1 = (gh-gl)/epsilon
            if counter > 0:
                return True, wolfe, q1, f1
            else:
                # even the first newton step failed, revert to back tracking
                pass

        # simple back tracking with tau = 0.5, no wolfe conditions yet
        q1 = -numpy.sign(g0)*last_step_size*1.5
        counter = 0.0
        while True:
            f1 = fun(q1)
            if f1 < f0:
                return True, False, q1, f1
            q1 *= 0.5
            counter += 1
            if counter > self.max_iter:
                # had enough iterations, line search fails
                return False, False, 0.0, f0


class ConvergenceCondition(object):
    """Callable object that tests the convergence of the minimizer"""
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
    def __init__(self, max_iter=None, fun_margin=None, grad_margin=None):
        """Initialize a StopLossCondition object

           Optional arguments:
             max_iter  --  the maximum number of iterations allowed
             fun_margin  --  if the function to be minimized goes above the
                             lowest value so far plus this margin, the
                             minimization is aborted
             grad_margin  --  if the RMS value of the gradient components goes
                              above the lowest value plus this threshold, the
                              minimization is aborted

           The present arguments define when the minimization has lost track.
        """
        self.max_iter = max_iter
        self.fun_margin = fun_margin
        self.grad_margin = grad_margin

        self.reset()

    def reset(self):
        self.fn_lowest = None
        self.grad_rms_lowest = None

    def __call__(self, counter, fn, gradient):
        """Return True when the minimizer has lost track"""
        if self.max_iter is not None and counter >= self.max_iter:
            return True

        if self.fun_margin is not None:
            if self.fn_lowest is None or fn < self.fn_lowest:
                self.fn_lowest = fn
            elif fn > self.fn_lowest + self.fun_margin:
                return True

        if self.grad_margin is not None:
            grad_rms = numpy.sqrt((gradient**2).mean())
            if self.grad_rms_lowest is None or grad_rms < self.grad_rms_lowest:
                self.grad_rms_lowest = grad_rms
            elif grad_rms > self.grad_rms_lowest + self.grad_margin:
                return True

        # all is fine
        return False


class FunctionWrapper(object):
    """Generic implementation of function and its analytical derivates

       Also supports function evaluation and derivative along a line.
    """
    def __init__(self, fun):
        """Initialize a FunctionWrapper object

           Argument:
             fun  --  a multivariate function
        """
        self._fun = fun
        self.axis = None
        self.x0 = None

    def set_line(self, x0, axis):
        """Configure the 1D function for a line search

           Arguments:
             x0  --  the reference point (q=0)
             axis  --  a unit vector in the direction of the line search
        """
        self.x0 = x0
        self.axis = axis
        self.line_cache = {}



class AnaGradWrapper(FunctionWrapper):
    """Wrapper of a function that also supports evaluation along a line

       This version assumes that the underlying function implements an
       analytical gradient.
    """
    def __init__(self, fun):
        """Initialize a AnaGradWrapper object

           Argument:
             fun  --  a multivariate function that can also compute analytical
                      derivatives, see below

           The function fun takes a mandatory argument x and an optional argument
           do_gradient:
             x  --  the arguments of the function to be tested
             do_gradient  --  when False, only the function value is returned.
                              when True, a 2-tuple with the function value and
                              the gradient are returned [default=False]
        """
        FunctionWrapper.__init__(self, fun)

    def __call__(self, x, do_gradient=False):
        """Just call the underlying function with the same arguments"""
        return self._fun(x, do_gradient=do_gradient)

    def line(self, q, do_gradient=False):
        """Evaluate the function along a predefined line

           The derivative is computed as the dot product of the axis and the
           analytical gradient.
        """
        x = self.x0 + self.axis*q
        if do_gradient:
            f, g = self(x, do_gradient=True)
            return f, numpy.dot(g, self.axis)
        else:
            return self(x)


class NumGradWrapper(FunctionWrapper):
    def __init__(self, fun, epsilon):
        """Initialize a NumGradWrapper object

           Argument:
             fun  --  a multivariate function, no analytical derivatives
                      required
             epsilon  --  small scalar used for finite differences

           The function fun only takes a mandatory argument x:
             x  --  the arguments of the function to be tested
        """
        self.epsilon = epsilon
        FunctionWrapper.__init__(self, fun)

    def __call__(self, x, do_gradient=False):
        """Call the underlying for function evaluates, use finite differences for gradient"""
        f = self._fun(x)
        if do_gradient:
            g = numpy.zeros(x.shape)
            for j in xrange(len(x)):
                xh = x.copy()
                xh[j] += 0.5*self.epsilon
                xl = x.copy()
                xl[j] -= 0.5*self.epsilon
                g[j] = (self._fun(xh) - self._fun(xl))/self.epsilon
            return f, g
        else:
            return f

    def line(self, q, do_gradient=False):
        """Evaluate the function along a predefined line

           The derivative is computed using finite differences
        """
        x = self.x0 + self.axis*q
        f = self._fun(x)
        if do_gradient:
            fh = self._fun(x + (0.5*self.epsilon) * self.axis)
            fl = self._fun(x - (0.5*self.epsilon) * self.axis)
            return f, (fh - fl)/self.epsilon
        else:
            return f


class Minimizer(object):
    """A flexible multivariate minimizer

       The minimizer searches (in principle) for the 'nearest' local minimum.
    """

    def __init__(self, x_init, fun, search_direction, line_search,
                 convergence_condition, stop_loss_condition, anagrad=False,
                 epsilon=1e-6, verbose=True, callback=None):
        """Initialize the minimizer

           Arguments:
             x_init  --  the initial guess for the minimum
             fun  --  function to be minimized (see below)
             search_direction  --  a SearchDirection object
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
             verbose  --  print progress information on screen [default=True]
             callback  --  optional callback routine after each CG iteration.
                           the callback routine gets the minimizer as first
                           and only argument. (default=None)

           The function fun takes a mandatory argument x and an optional argument
           do_gradient:
             x  --  the arguments of the function to be tested
             do_gradient  --  when False, only the function value is returned.
                              when True, a 2-tuple with the function value and
                              the gradient are returned [default=False]
        """
        if len(x_init.shape)!=1:
            raise ValueError("The unknowns must be stored in a plain row vector.")
        # self.x always contains the current parameters
        self.x = x_init.copy()
        if anagrad:
            self.fun = AnaGradWrapper(fun)
        else:
            self.fun = NumGradWrapper(fun, epsilon)
        self.search_direction = search_direction
        self.line_search = line_search
        self.convergence_condition = convergence_condition
        self.stop_loss_condition = stop_loss_condition
        self.callback = callback
        self.epsilon = epsilon
        self.verbose = verbose

        # perform some resets:
        self.search_direction.reset()
        self.stop_loss_condition.reset()

        # the current function value
        self.f = None
        # the current gradient
        self.gradient = None
        # the current step size
        self.step = None

        self._iterate()

    def _iterate(self):
        """Run the iterative optimizer"""
        self.f, self.gradient = self.fun(self.x, do_gradient=True)
        last_end = time.clock()
        # the cg loop
        self.counter = 0
        while True:
            # compute the new direction
            self.search_direction.update(self.gradient, self.step)
            # print some stuff on screen
            if self.counter % 20 == 0:
                self._print_header()
            self._screen("% 5i   %2s  " % (self.counter, self.search_direction.status))
            # perform a line search
            line_success = self._line_opt()
            if not line_success:
                self._screen("Line search failed", newline=True)
                if self.search_direction.is_sd():
                    self._screen("LINE FAILED", newline=True)
                    return
                else:
                    self.search_direction.reset()
                    continue
            # compute the gradient at the new point
            self.f, self.gradient = self.fun(self.x, do_gradient=True)
            # check convergence
            converged, status = self.convergence_condition(self.gradient, self.step)
            self._screen(status)
            # timing
            end = time.clock()
            self._screen("%5.2f" % (end - last_end), newline=True)
            last_end = end
            # check convergence, part 2
            if converged:
                self._screen("CONVERGED", newline=True)
                return
            # check stop loss
            lost = self.stop_loss_condition(self.counter, self.f, self.gradient)
            if lost:
                self._screen("LOST", newline=True)
                return
            # call back
            if self.callback is not None:
                self.callback(self)
            self.counter += 1

    def _print_header(self):
        header = " Iter  Dir     Function"
        header += self.convergence_condition.get_header()
        header += "    Time"
        self._screen("-"*(len(header)+2), newline=True)
        self._screen(header, newline=True)
        self._screen("-"*(len(header)+2), newline=True)

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
        direction_norm = numpy.linalg.norm(direction)
        if direction_norm == 0:
            return False
        self.fun.set_line(self.x, direction/direction_norm)

        if self.step is None:
            last_step_size = 1.0
        else:
            last_step_size = numpy.linalg.norm(self.step)
        success, wolfe, qopt, fopt = self.line_search(self.fun.line, last_step_size, self.epsilon)
        if success:
            self.step = qopt*self.fun.axis
            self.x = self.x + self.step
            self.f = fopt
            self._screen("% 9.3e  " % self.f)
            if not wolfe:
                self.search_direction.reset()
            return True
        else:
            self._reset_state()
            return False

    def _reset_state(self):
        """Reset of the internal state of the function"""
        self.fun(self.x) # reset the internal state of the function


def check_anagrad(fun, x0, epsilon, scale=10):
    """Check the analytical gradient using finite differences

       Arguments:
         fun  --  the function to be tested, more info below
         x0  --  the reference point around which the function should be tested
         epsilon  --  a small scalar used for the finite differences

       Optional argument
         scale  --  scale*epsilon is the threshold for an error between the
                    analytical and the numerical gradient

       The function fun takes a mandatory argument x and an optional argument
       do_gradient:
         x  --  the arguments of the function to be tested
         do_gradient  --  when False, only the function value is returned. when
                          True, a 2-tuple with the function value and the
                          gradient are returned [default=False]
    """
    N = len(x0)
    f0, ana_grad = fun(x0, do_gradient=True)
    for i in xrange(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        num_grad_comp = (fun(xh)-fun(xl))/epsilon
        if abs(num_grad_comp - ana_grad[i]) > scale*epsilon:
            raise ValueError("Error in the analytical gradient, component %i, got %s, should be about %s" % (i, ana_grad[i], num_grad_comp))


def compute_fd_hessian(fun, x0, epsilon, anagrad=True):
    """Compute the Hessian using the finite difference method

       Arguments:
         fun  --  the function for which the Hessian should be computed, more
                  info below
         x0  --  the point at which the hessian must be computed
         epsilon  --  a small scalar step size used to compute the finite
                      differences

       Optional argument:
         anagrad  --  when True, analytical gradients are used [default=True]

       The function fun takes a mandatory argument x and an optional argument
       do_gradient:
         x  --  the arguments of the function to be tested
         do_gradient  --  when False, only the function value is returned. when
                          True, a 2-tuple with the function value and the
                          gradient are returned [default=False]
    """
    N = len(x0)

    def compute_gradient(x):
        if anagrad:
            return fun(x, do_gradient=True)[1]
        else:
            gradient = numpy.zeros(N, float)
            for i in xrange(N):
                xh = x.copy()
                xh[i] += 0.5*epsilon
                xl = x.copy()
                xl[i] -= 0.5*epsilon
                gradient[i] = (fun(xh)-fun(xl))/epsilon
            return gradient

    hessian = numpy.zeros((N,N), float)
    for i in xrange(N):
        xh = x0.copy()
        xh[i] += 0.5*epsilon
        xl = x0.copy()
        xl[i] -= 0.5*epsilon
        hessian[i] = (compute_gradient(xh) - compute_gradient(xl))/epsilon

    return 0.5*(hessian + hessian.transpose())

