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

class JacobianAnalysis(object):
    def __init__(self, internal_coordinates, carthesian_values):
        self.internal_coordinates = internal_coordinates
        self.carthesian_values = carthesian_values
        
        jacobian = []
        values = []
        for internal_coordinate in self.internal_coordinates:
            value, derivates = internal_coordinate(carthesian_values)
            values.append(value)
            jacobian.append(Numeric.ravel(derivates))
            
        self.values = Numeric.array(values)
        self.jacobian = Numeric.transpose(Numeric.array(jacobian))


class GradientAnalysis(JacobianAnalysis):
    def __init__(self, internal_coordinates, carthesian_values, job):
        JacobianAnalysis.__init__(self, internal_coordinates, carthesian_values)
        self.job = job
        self.energy = job.energy
        self.energy_error = job.energy_accuracy

        self.V, self.S, self.Wt = LinearAlgebra.singular_value_decomposition(self.jacobian, True)
        self.W = Numeric.transpose(self.Wt)
        self.rank = sum(self.S>(max(self.S)*1e-7))
        self.S = self.S[:self.rank]
        self.V = self.V[:,:self.rank]
        
        particular_transform = Numeric.dot(self.W[:,:self.rank], Numeric.transpose(self.V/self.S))
        
        self.gradient = Numeric.ravel(job.gradient)
        self.particular = Numeric.dot(particular_transform, self.gradient)
        self.gradient_error = Numeric.ones(self.gradient.shape, Numeric.Float)*job.gradient_accuracy
        self.particular_error = Numeric.sqrt(Numeric.dot(particular_transform**2, self.gradient_error**2))

        if self.W.shape[1] > self.rank:
            self.nullspace = self.W[:,self.rank:]



class HessianAnalysis(JacobianAnalysis):
    def __init__(self, internal_coordinates, carthesian_values, job):
        JacobianAnalysis.__init__(self, internal_coordinates, carthesian_values)
        self.job = job
        self.hessian = job.hessian
        
        self.normalized_jacobian = self.jacobian.copy()
        for col in Numeric.transpose(self.normalized_jacobian):
            col[:] /= math.sqrt(Numeric.dot(col, col))
        self.overlap = Numeric.dot(Numeric.transpose(self.normalized_jacobian), self.normalized_jacobian)

        self.diag = []
        for col in Numeric.transpose(self.normalized_jacobian):
            temp = Numeric.dot(self.hessian, col)
            norm_transf = math.sqrt(Numeric.dot(temp, temp))
            self.diag.append(Numeric.dot(col, temp)/norm_transf)
    
