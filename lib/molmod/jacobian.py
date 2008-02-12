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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --

# Tools related to solving the jacobian system

import numpy, numpy.linalg, math

def jacobian_analysis(configuration, internal_coordinates):
    configuration.internal_coordinates = internal_coordinates

    jacobian = []
    values = []
    for internal_coordinate in internal_coordinates:
        value, tangent = internal_coordinate.value_tangent(configuration.cartesian_values)
        values.append(value)
        jacobian.append(numpy.ravel(tangent))

    configuration.internal_values = numpy.array(values)
    configuration.jacobian = numpy.transpose(numpy.array(jacobian))


def energy_analysis(configuration, internal_coordinates):
    jacobian_analysis(configuration, internal_coordinates)
    configuration.energy_error = configuration.energy_accuracy

    V, S, Wt = numpy.linalg.svd(configuration.jacobian, True)
    W = numpy.transpose(Wt)
    rank = sum(abs(S)>(max(abs(S))*1e-7))
    configuration.rank = rank
    S = S[:rank]
    V = V[:,:rank]

    particular_transform = numpy.dot(W[:,:rank], numpy.transpose(V/S))

    configuration.gradient = numpy.ravel(configuration.gradient)
    configuration.particular = numpy.dot(particular_transform, configuration.gradient)
    configuration.gradient_error = numpy.ones(configuration.gradient.shape, float)*configuration.gradient_accuracy
    configuration.particular_error = numpy.sqrt(numpy.dot(particular_transform**2, configuration.gradient_error**2))

    if W.shape[1] > rank:
        configuration.nullspace = W[:,rank:]


