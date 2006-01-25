# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

# Tools related to solving the jacobian system

import Numeric, LinearAlgebra, math

def jacobian_analysis(configuration, internal_coordinates):
    configuration.internal_coordinates = internal_coordinates
    
    jacobian = []
    values = []
    for internal_coordinate in internal_coordinates:
        value, tangent = internal_coordinate.value_tangent(configuration.carthesian_values)
        values.append(value)
        jacobian.append(Numeric.ravel(tangent))
        
    configuration.internal_values = Numeric.array(values)
    configuration.jacobian = Numeric.transpose(Numeric.array(jacobian))


def energy_analysis(configuration, internal_coordinates):
    jacobian_analysis(configuration, internal_coordinates)
    configuration.energy_error = configuration.energy_accuracy

    V, S, Wt = LinearAlgebra.singular_value_decomposition(configuration.jacobian, True)
    W = Numeric.transpose(Wt)
    rank = sum(abs(S)>(max(abs(S))*1e-7))
    configuration.rank = rank
    S = S[:rank]
    V = V[:,:rank]
    
    particular_transform = Numeric.dot(W[:,:rank], Numeric.transpose(V/S))
    
    configuration.gradient = Numeric.ravel(configuration.gradient)
    configuration.particular = Numeric.dot(particular_transform, configuration.gradient)
    configuration.gradient_error = Numeric.ones(configuration.gradient.shape, Numeric.Float)*configuration.gradient_accuracy
    configuration.particular_error = Numeric.sqrt(Numeric.dot(particular_transform**2, configuration.gradient_error**2))

    if W.shape[1] > rank:
        configuration.nullspace = W[:,rank:]
