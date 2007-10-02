# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from molmod.constants import boltzman
import numpy, pickle


class Error(Exception):
    pass


class TrajectoryMaker(object):
    def __init__(self, calculate_energy, calculate_gradient, dof, log):
        self.initialized = False

        self.calculate_energy = calculate_energy
        self.calculate_gradient = calculate_gradient
        self.dof = dof
        self.log = log

        self.current_state = numpy.zeros(self.dof, float)

        self.trajectory_vars = []
        self.init_trajectory_vars(["states", "potential_energies", "gradients"])

    def init_trajectory_vars(self, names):
        self.trajectory_vars.extend(names)

    def clear_trajectory_vars(self):
        for var in self.trajectory_vars:
            self.__dict__[var] = []

    def dump_trajectory(self, filename):
        f = file(filename, "w")
        pickle.dump(dict((var, numpy.array(self.__dict__[var])) for var in self.trajectory_vars), f)
        f.close()

    def verbose(self):
        if self.log is not None:
            self.log(self)

    def run(self):
        if not self.initialized:
            raise Error("Please first initialize the state of the integrator.")


class Fix(object):
    def __init__(self, skip):
        self.skip = skip

    def active(self, vi):
        return len(vi.states) % self.skip == 0

    def run(self, vi):
        raise NotImplementedError


class VelocityScaler(Fix):
    def __init__(self, skip, temperature):
        self.ke = temperature*boltzman
        Fix.__init__(self, skip)

    def run(self, vi):
        if self.active(vi):
            factor = numpy.sqrt(self.ke*vi.dof/(0.5*(vi.masses*vi.tmp**2).sum()))
            print "Velocity rescale", factor
            vi.tmp *= factor



class VerletIntegrator(TrajectoryMaker):
    def __init__(self, calculate_energy, calculate_gradient, log, masses, time_step, fixes=[]):
        TrajectoryMaker.__init__(self, calculate_energy, calculate_gradient, len(masses), log)

        self.masses = masses
        self.time_step = time_step
        self.fixes = fixes

        self.current_velocities = numpy.zeros(self.dof, float)
        self.accelerations = numpy.zeros(self.dof, float)
        self.tmp = numpy.zeros(self.dof, float)
        self.gradient = numpy.zeros(self.dof, float)

        self.init_trajectory_vars(["velocities", "energies", "kinetic_energies", "times", "temperatures"])

    def initialize_state(self, initial_state, initial_velocities=0):
        self.clear_trajectory_vars()

        self.current_state[:] = initial_state
        self.current_velocities[:] = initial_velocities
        self.gradient[:] = self.calculate_gradient(self.current_state)
        self.accelerations[:] = -self.gradient/self.masses
        self.tmp[:] = self.current_velocities + 0.5*self.accelerations*self.time_step
        self.initialized = True

        self.log_state()

    def log_state(self):
        self.states.append(self.current_state.copy())
        self.velocities.append(self.current_velocities.copy())

        potential_energy = self.calculate_energy(self.current_state)
        kinetic_energy = 0.5*(self.masses*self.current_velocities**2).sum()
        self.potential_energies.append(potential_energy)
        self.kinetic_energies.append(kinetic_energy)
        self.energies.append(potential_energy + kinetic_energy)
        step = len(self.times)
        self.times.append(step*self.time_step)
        self.temperatures.append(kinetic_energy/boltzman/self.dof)
        self.gradients.append(self.gradient)

        self.verbose()

    def run(self, steps):
        TrajectoryMaker.run(self)
        try:
            for step in xrange(steps):
                self.current_state[:] += self.tmp*self.time_step
                self.gradient[:] = self.calculate_gradient(self.current_state)
                self.accelerations[:] = -self.gradient/self.masses
                self.current_velocities[:] = self.tmp + 0.5*self.accelerations[:]*self.time_step
                self.tmp[:] = self.current_velocities + 0.5*self.accelerations[:]*self.time_step

                self.log_state()

                for fix in self.fixes:
                    fix.run(self)
        except KeyboardInterrupt:
            print "The iteration has been interupted."

    def initialize_kinetic_energy(self, initial_state, kinetic_energy, exact=False):
        initial_velocities = numpy.random.normal(0.0, 1.0, self.dof)*numpy.sqrt(2*kinetic_energy/self.masses)
        if exact:
            factor = numpy.sqrt(kinetic_energy*self.dof/(0.5*(self.masses*initial_velocities**2).sum()))
            initial_velocities *= factor
        self.initialize_state(initial_state, initial_velocities)

    def initialize_temperature(self, initial_state, temperature, exact=False):
        self.initialize_kinetic_energy(initial_state, temperature*boltzman, exact)


class ConjugateGradientOptimizer(TrajectoryMaker):
    def __init__(self, calculate_energy, calculate_gradient, dof, log, ofm, epsilon):
        TrajectoryMaker.__init__(self, calculate_energy, calculate_gradient, dof, log)
        self.ofm = ofm
        self.epsilon = epsilon

        self.init_trajectory_vars(["gradient_norms"])

    def initialize_state(self, initial_state):
        self.clear_trajectory_vars()

        self.current_state[:] = initial_state
        self.states.append(initial_state)
        self.initialized = True

    def run(self, steps, gradient_norm_threshold):
        TrajectoryMaker.run(self)
        gradient = self.calculate_gradient(self.current_state)
        direction = -gradient.copy()
        gradient_norm = numpy.sqrt(numpy.dot(gradient, gradient))
        try:
            for step in xrange(steps):
                state_delta = self.current_state + self.epsilon*direction/gradient_norm
                gradient_delta = self.calculate_gradient(state_delta)
                hessian_applied = (gradient_delta - gradient)/self.epsilon
                second_order_deriv_in_direction = numpy.dot(hessian_applied, direction)
                if second_order_deriv_in_direction > 0:
                    # hier een CG stap
                    self.last_method = "CG"
                    step_size = gradient_norm**2/second_order_deriv_in_direction
                    self.current_state[:] += step_size*direction/gradient_norm
                    new_gradient = self.calculate_gradient(self.current_state)
                    new_gradient_norm = numpy.sqrt(numpy.dot(new_gradient, new_gradient))
                    #beta = (new_gradient_norm*new_gradient_norm)/(gradient_norm*gradient_norm)
                    beta = numpy.dot(new_gradient - gradient, new_gradient)/(gradient_norm*gradient_norm)
                    direction *= beta
                    direction -= new_gradient
                    gradient = new_gradient
                    gradient_norm = new_gradient_norm
                else:
                    # hier een ruwe SD stap
                    self.last_method = "SD"
                    self.current_state[:] += self.ofm*direction/gradient_norm
                    gradient = self.calculate_gradient(self.current_state)
                    gradient_norm = numpy.sqrt(numpy.dot(gradient, gradient))
                self.gradient_norms.append(gradient_norm)
                self.potential_energies.append(self.calculate_energy(self.current_state))
                self.gradients.append(gradient)
                self.states.append(self.current_state.copy())
                self.verbose()

                if gradient_norm < gradient_norm_threshold and self.last_method == "CG":
                    break
        except KeyboardInterrupt:
            print "The iteration has been interupted."


