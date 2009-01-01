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


import numpy


__all__ = [
    "GoldenLineSearch", "NewtonLineSearch", "NewtonGLineSearch", "Minimizer"
]


phi = 0.5*(1+numpy.sqrt(5))


class LineSearch(object):
    def __init__(self, qtol, max_step, max_iter):
        self.qtol = qtol
        self.max_step = max_step
        self.max_iter = max_iter

    def limit_step(self, step):
        if step > self.max_step: return self.max_step
        if step < -self.max_step: return -self.max_step
        return step

    def get_extra_log_dtypes(self):
        return []

    def get_extra_log(self):
        return ()



class GoldenLineSearch(LineSearch):
    def __init__(self, qtol, max_step, max_iter):
        LineSearch.__init__(self, qtol, max_step, max_iter)
        self.num_bracket = 0
        self.num_golden = 0

    def get_extra_log_dtypes(self):
        return [("num_bracket", int), ("num_golden", int)]

    def get_extra_log(self):
        result = (self.num_bracket, self.num_golden)
        self.num_bracket = 0
        self.num_golden = 0
        return result

    def __call__(self, fun, f0, last_step_size, epsilon):
        # bracket the minimum
        triplet = self.bracket(last_step_size, f0, fun)
        if triplet is None:
            return False, 0.0, f0
        # do a golden section optimization
        qopt, fopt = self.golden(triplet, fun)
        qopt = self.limit_step(qopt)
        fopt = fun(qopt)
        return True, qopt, fopt

    def bracket(self, qinit, f0, fun):
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

    def golden(self, triplet, fun):
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
    def __init__(self, qtol, max_step, max_iter):
        LineSearch.__init__(self, qtol, max_step, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        return [("num_reduce", int)]

    def get_extra_log(self):
        result = (self.num_reduce,)
        self.num_reduce = 0
        return result

    def __call__(self, fun, f0, last_step_size, epsilon):
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
    def __init__(self, qtol, max_step, max_iter):
        LineSearch.__init__(self, qtol, max_step, max_iter)
        self.num_reduce = 0

    def get_extra_log_dtypes(self):
        return [("num_reduce", int)]

    def get_extra_log(self):
        result = (self.num_reduce,)
        self.num_reduce = 0
        return result

    def __call__(self, fun, f0, last_step_size, epsilon):
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
    def __init__(
        self, x_init, fun, LineSearchCls, ftol, xtol, max_step, max_iter, max_line_iter,
        do_gradient=False, epsilon_init=1e-6, absftol=False, verbose=True, callback=None,
        min_iter=0, extra_log_dtypes=None
    ):
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
        self.step_size = 1
        self.direction_cg = numpy.zeros(self.x.shape, float)
        self.direction_sd = numpy.zeros(self.x.shape, float)
        self.direction_sd_old = numpy.zeros(self.x.shape, float)

        self.iterate(max_iter)

    def screen(self, s, newline=True):
        if self.verbose:
            if newline:
                print s
            else:
                print s,

    def step_cg(self, label):
        return self.line_opt(self.direction_cg, label)

    def update_cg(self):
        self.direction_sd_old[:] = self.direction_sd
        self.update_direction_sd()
        self.beta = (
            numpy.dot(self.direction_sd, self.direction_sd - self.direction_sd_old) /
            numpy.dot(self.direction_sd_old, self.direction_sd_old)
        )
        if self.beta < 0:
            self.beta = 0
        self.direction_cg[:] *= self.beta
        self.direction_cg[:] += self.direction_sd

    def update_sd(self):
        self.update_direction_sd()
        self.direction_cg[:] = self.direction_sd
        self.beta = 0

    #def line_scan_global(self):
    #    delta_xs = numpy.arange(self.num_line_scan, dtype=float)/(self.num_line_scan-1)*2*self.max_step - self.max_step
    #    fs = []
    #    for delta_x in delta_xs:
    #        fs.append(self.fun(self.x + delta_x))
    #    fs = numpy.array(fs)
    #    m = fs.argmin()
    #    if (self.f-fs[m])/self.f < 1e-16:
    #        print "  'Global scan'         No Improvement"
    #        self.fun(self.x) # reset the internal state of the function
    #        return False
    #    print "  'Global scan'        ",
    #    return self.handle_new_solution(self.x + delta_xs[m], fs[m])

    #def line_scan_parameter(self, index=None):
    #    delta_xs = numpy.arange(self.num_line_scan, dtype=float)/(self.num_line_scan-1)*2*self.max_step - self.max_step
    #    xnew = self.x.copy()
    #    fs = []
    #    for delta_x in delta_xs:
    #        xnew[:] = self.x
    #        xnew[index] += delta_x
    #        fs.append(self.fun(xnew))
    #    fs = numpy.array(fs)
    #    m = fs.argmin()
    #    if (self.f-fs[m])/self.f < 1e-16:
    #        print "  'Parameter %03i scan'  No Improvement" % index
    #        self.fun(self.x) # reset the internal state of the function
    #        return False
    #    print "  'Parameter %03i scan' " % index,
    #    xnew[:] = self.x
    #    xnew[index] += delta_xs[m]
    #    return self.handle_new_solution(xnew, fs[m])

    def line_opt(self, direction, label):
        if numpy.linalg.norm(direction) == 0:
            self.screen("  'Opt% 3s'              Direction zero" % label)
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
            self.screen("  'Opt% 3s'             " % label, False)
            return self.handle_new_solution(self.x + qopt*direction, fopt)
        else:
            # returning the old x and f, means that the optimization failed
            self.screen("  'Opt% 3s'              Line search failed" % label)
            self.fun(self.x) # reset the internal state of the function
            return False


    def update_direction_sd(self):
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
            self.fun(self.x) # reset the internal state of the function
            #print self.direction_sd
        #sys.exit()

    def handle_new_solution(self, xnew, fnew):
        self.step_size = numpy.linalg.norm(self.x - xnew)
        if self.absftol:
            self.decrease = self.f - fnew
            self.screen("step=% 9.7e   absdecr=% 9.7e   fnew=% 9.7e" % (self.step_size, self.decrease, fnew))
        else:
            self.decrease = (self.f - fnew)/self.f
            self.screen("step=% 9.7e   reldecr=% 9.7e   fnew=% 9.7e" % (self.step_size, self.decrease, fnew))
        if self.step_size < self.xtol:
            result = False
        if self.decrease < self.ftol:
            result = False
        else:
            result = True
        if self.decrease > 0:
            self.x = xnew
            self.f = fnew
        else:
            self.fun(self.x) # reset the internal state of the function
        return result

    def append_log(self):
        extra_fields = self.line_search.get_extra_log()
        if self.callback is not None:
            extra_fields = extra_fields + tuple(self.callback(self.x))
        self.log.append((
            self.f, self.step_size, self.beta,
            numpy.linalg.norm(self.direction_sd),
            numpy.linalg.norm(self.direction_cg),
            self.x.copy(),
        ) + extra_fields)

    def iterate(self, max_iter):
        self.f = self.fun(self.x)

        self.beta = 0
        self.update_direction_sd()
        self.direction_cg[:] = self.direction_sd
        last_reset = True
        cg_lower = False
        # the cg loop
        for counter in xrange(max_iter):
            self.screen("Iter % 5i of % 5i" % (counter, max_iter), False)
            if last_reset:
                lower = self.step_cg("SD")
            else:
                lower = self.step_cg("CG")
            self.append_log()
            lower = lower or (counter < self.min_iter and self.decrease > 0)
            if lower:
                self.update_cg()
                last_reset = False
            else:
                if last_reset:
                    break
                self.update_sd()
                last_reset = True
            if self.step_size > 0:
                self.epsilon = min(self.step_size*1e-2, self.epsilon_max)

            #if self.do_scan:
            #    # the parameter scan loop
            #    self.num_bracket = 0
            #    self.num_golden = 0
            #    scan_lower = False
            #    print "Iter % 5i of % 5i" % (counter, max_iter),
            #    scan_lower |= self.line_scan_global()
            #    self.append_log(False)
            #    counter += 1
            #    for parameter_index in xrange(len(self.x)):
            #        print "Iter % 5i of % 5i" % (counter, max_iter),
            #        scan_lower |= self.line_scan_parameter(parameter_index)
            #        self.append_log(False)
            #        counter += 1

            #    if not scan_lower:
            #        self.set_coarse(False)

        self.screen("Done")
        return self.x

    def get_log(self):
        return numpy.array(self.log, [
            ("f", float), ("step_size", float), ("beta", float),
            ("norm_sd", float), ("norm_cg", float), ("x", float, self.x.shape),
        ] + self.extra_log_dtypes)


