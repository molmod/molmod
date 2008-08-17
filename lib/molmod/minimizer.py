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


__all__ = ["Minimizer", "minimize"]


phi = 0.5*(1+numpy.sqrt(5))


class Minimizer(object):
    def __init__(
        self, x_init, fun, ftol_fine, ftol_coarse, xtol, epsilon_init=1e-3, 
        num_line_scan=61, max_step=5.0, callback=None, max_iter_cg=100, 
        do_scan=False
    ):
        self.x = x_init.copy()
        self.fun = fun
        self.ftol_coarse = ftol_coarse
        self.ftol_fine = ftol_fine
        self.xtol = xtol
        self.epsilon = epsilon_init
        self.step_size = 1

        self.num_line_scan = num_line_scan
        self.max_step = max_step
        self.callback = callback
        self.max_iter_cg = max_iter_cg
        self.do_scan = do_scan

        self.log = []
        self.direction_cg = numpy.zeros(self.x.shape, float)
        self.direction_sd = numpy.zeros(self.x.shape, float)
        self.direction_sd_old = numpy.zeros(self.x.shape, float)

    def set_coarse(self, coarse):
        self.coarse = coarse
        if coarse:
            self.ftol = self.ftol_coarse
        else:
            self.ftol = self.ftol_fine

    def opt_cg(self, label):
        return self.line_opt(self.direction_cg, label)

    def update_cg(self):
        self.direction_sd_old[:] = self.direction_sd
        self.update_direction_sd()
        self.beta = (
            numpy.dot(self.direction_sd, self.direction_sd) /# - self.direction_sd_old) /
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

    def line_scan_global(self):
        delta_xs = numpy.arange(self.num_line_scan, dtype=float)/(self.num_line_scan-1)*2*self.max_step - self.max_step
        fs = []
        for delta_x in delta_xs:
            fs.append(self.fun(self.x + delta_x))
        fs = numpy.array(fs)
        m = fs.argmin()
        if (self.f-fs[m])/self.f < 1e-16:
            print "  'Global scan'         No Improvement"
            self.fun(self.x) # reset the internal state of the function
            return False
        print "  'Global scan'        ",
        return self.handle_new_solution(self.x + delta_xs[m], fs[m])

    def line_scan_parameter(self, index=None):
        delta_xs = numpy.arange(self.num_line_scan, dtype=float)/(self.num_line_scan-1)*2*self.max_step - self.max_step
        xnew = self.x.copy()
        fs = []
        for delta_x in delta_xs:
            xnew[:] = self.x
            xnew[index] += delta_x
            fs.append(self.fun(xnew))
        fs = numpy.array(fs)
        m = fs.argmin()
        if (self.f-fs[m])/self.f < 1e-16:
            print "  'Parameter %03i scan'  No Improvement" % index
            self.fun(self.x) # reset the internal state of the function
            return False
        print "  'Parameter %03i scan' " % index,
        xnew[:] = self.x
        xnew[index] += delta_xs[m]
        return self.handle_new_solution(xnew, fs[m])

    def update_direction_sd(self):
        tmp = self.x.copy()
        for j in xrange(len(self.x)):
            tmp[j] += self.epsilon
            self.direction_sd[j] = -(self.fun(tmp) - self.f)/self.epsilon
            tmp[j] = self.x[j]
        self.fun(self.x) # reset the internal state of the function

    def line_opt(self, direction, label):
        if numpy.linalg.norm(direction) == 0:
            print "  'Opt% 3s'              Direction zero" % label
            return False 
        direction = direction/numpy.linalg.norm(direction)
        def fun_aux(q):
            xq = self.x + q*direction
            return self.fun(xq)

        # bracket the minimum
        triplet = self.bracket(self.step_size, self.f, fun_aux)
        if triplet is None:
            # returning the old x and f, means that the optimization failed
            print "  'Opt% 3s'              No bracket found" % label
            self.fun(self.x) # reset the internal state of the function
            return False
        # do a golden section optimization
        qopt, fopt = self.golden(triplet, fun_aux)
        if qopt > self.max_step: qopt = self.max_step
        if qopt < -self.max_step: qopt = -self.max_step
        fopt = fun_aux(qopt)
        print "  'Opt% 3s'             " % label,
        return self.handle_new_solution(self.x + qopt*direction, fopt)

    def handle_new_solution(self, xnew, fnew):
        self.step_size = numpy.linalg.norm(self.x - xnew)
        relative_decrease = (self.f - fnew)/self.f
        print "step=% 9.7e   reldecr=% 9.7e   fnew=% 9.7e" % (self.step_size, relative_decrease, fnew)
        if self.step_size < self.xtol:
            result = False
        if relative_decrease < self.ftol:
            result = False
        else:
            result = True
        if relative_decrease > 0:
            self.x = xnew
            self.f = fnew
        else:
            self.fun(self.x) # reset the internal state of the function
        return result

    def bracket(self, qinit, f0, fun):
        self.num_bracket = 0
        qa = qinit
        fa = fun(qa)
        if fa >= f0:
            while True:
                self.num_bracket += 1
                #print "    bracket shrink"
                qb, fb = qa, fa
                qa /= 1+phi
                fa = fun(qa)
                if qa < self.xtol:
                    return
                if fa < f0:
                    return (0, f0), (qa, fa), (qb, fb)
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
            if abs(qa-qb) < self.xtol:
                return qc, fc

    def append_log(self, cg):
        if self.callback is not None:
            extra_fields = tuple(self.callback())
        else:
            extra_fields = ()
        self.log.append((
            self.f, self.step_size, self.beta,
            numpy.linalg.norm(self.direction_sd),
            numpy.linalg.norm(self.direction_cg),
            cg, self.num_bracket, self.num_golden,
            self.x.copy(),
        ) + extra_fields)

    def iterate(self):
        self.f = self.fun(self.x)
        print "Init      ",
        self.set_coarse(True)
        self.line_scan_global()
        callback()

        counter = 0
        while True:
            # the cg loop
            self.beta = 0
            self.num_golden = 0
            self.num_bracket = 0
            self.update_direction_sd()
            self.direction_cg[:] = self.direction_sd
            last_sd = True
            cg_lower = False
            for i in xrange(self.max_iter_cg):
                print "Iter % 5i" % counter,
                counter += 1
                if last_sd:
                    lower = self.opt_cg("SD")
                else:
                    lower = self.opt_cg("CG")
                cg_lower |= lower
                self.append_log(True)
                if lower:
                    self.update_cg()
                    last_sd = False
                else:
                    if last_sd:
                        break
                    self.update_sd()
                    last_sd = True
                callback()
                if self.step_size > 0:
                    self.epsilon = self.step_size*1e-2
            
            if not cg_lower and not self.coarse:
                break

            if self.coarse and last_sd:
                self.set_coarse(False)
            
            if self.do_scan:
                # the parameter scan loop
                self.num_bracket = 0
                self.num_golden = 0
                scan_lower = False
                print "Iter % 5i" % counter,
                scan_lower |= self.line_scan_global()
                self.append_log(False)
                counter += 1
                for parameter_index in xrange(len(self.x)):
                    print "Iter % 5i" % counter,
                    scan_lower |= self.line_scan_parameter(parameter_index)
                    self.append_log(False)
                    callback()
                    counter += 1
                
                if not scan_lower:
                    self.set_coarse(False)

        print "Done"
        return self.x

    def get_log(self):
        return numpy.array(self.log, [
            ("f", float), ("step_size", float), ("beta", float),
            ("norm_sd", float), ("norm_cg", float), ("cg", bool),
            ("num_bracket", int), ("num_golden", int),
            ("x", float, self.x.shape), ("solution", float, self.fun.solution.shape),
        ])


def minmize(*args, **kwargs):
    minmizer = Minimizer(*args, **kwargs)
    x_opt = minmizer.iterate
    return x_opt, minimizer


