# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.constants import boltzman
import numpy


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

        self.states = []
        self.potential_energies = []

    def verbose(self):
        if self.log is not None:
            self.log(self)

    def run(self):
        if not self.initialized:
            raise Error("Please first initialize the state of the integrator.")


class VerletIntegrator(TrajectoryMaker):
    def __init__(self, calculate_energy, calculate_gradient, log, masses, time_step):
        TrajectoryMaker.__init__(self, calculate_energy, calculate_gradient, len(masses), log)

        self.masses = masses
        self.time_step = time_step

        self.current_velocities = numpy.zeros(self.dof, float)
        self.accelerations = numpy.zeros(self.dof, float)
        self.tmp = numpy.zeros(self.dof, float)

        self.velocities = []
        self.energies = []
        self.kinetic_energies = []

    def initialize_state(self, initial_state, initial_velocities=0):
        self.current_state[:] = initial_state
        self.current_velocities[:] = initial_velocities[:]
        self.accelerations[:] = -self.calculate_gradient(self.current_state)/self.masses
        self.tmp[:] = self.current_velocities + 0.5*self.accelerations*self.time_step
        self.initialized = True

    def run(self, steps):
        TrajectoryMaker.run(self)
        for step in xrange(steps):
            self.current_state[:] += self.tmp*self.time_step
            self.accelerations[:] = -self.calculate_gradient(self.current_state)/self.masses
            self.current_velocities[:] = self.tmp + 0.5*self.accelerations[:]*self.time_step
            self.tmp[:] = self.current_velocities + 0.5*self.accelerations[:]*self.time_step

            self.states.append(self.current_state.copy())
            self.velocities.append(self.current_velocities.copy())

            potential_energy = self.calculate_energy(self.current_state)
            kinetic_energy = 0.5*(self.masses*self.current_velocities**2).sum()
            self.energies.append(potential_energy + kinetic_energy)
            self.potential_energies.append(potential_energy)
            self.kinetic_energies.append(kinetic_energy)
            self.energies.append(potential_energy + kinetic_energy)

            self.verbose()

    def initialize_kinetic_energy(self, initial_state, kinetic_energy):
        initial_velocities = numpy.random.normal(0.0, 1.0)*numpy.sqrt(2*kinetic_energy/self.masses)
        self.initialize_state(initial_state, initial_velocities)

    def initialize_temperature(self, initial_state, temperature):
        self.initialize_kinetic_energy(initial_state, temperature*boltzman)


class ConjugateGradientOptimizer(TrajectoryMaker):
    def __init__(self, calculate_energy, calculate_gradient, dof, log, ofm, epsilon):
        TrajectoryMaker.__init__(self, calculate_energy, calculate_gradient, dof, log)
        self.ofm = ofm
        self.epsilon = epsilon

        self.gradient_norms = []

    def initialize_state(self, initial_state):
        self.current_state[:] = initial_state
        self.initialized = True

    def run(self, steps):
        TrajectoryMaker.run(self)
        gradient = self.calculate_gradient(self.current_state)
        direction = -gradient.copy()
        gradient_norm = numpy.sqrt(numpy.dot(gradient, gradient))
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
                beta = (new_gradient_norm*new_gradient_norm)/(gradient_norm*gradient_norm)
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
            self.states.append(self.current_state.copy())
            self.verbose()
