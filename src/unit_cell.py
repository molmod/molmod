# Zeobuilder is an extensible GUI-toolkit for molecular model construction.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of Zeobuilder.
# 
# Zeobuilder is free software; you can redistribute it and/or
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


from molmod.units import angstrom
from molmod.vectors import random_orthonormal

import numpy

import math, copy


__all__ = ["check_cell", "UnitCell"]


def check_cell(cell, norm_threshold=1e-6, volume_threshold=1e-6):
    cell = cell.copy()
    for col, name in enumerate(["A", "B", "C"]):
        norm = math.sqrt(numpy.dot(cell[:,col], cell[:,col]))
        if norm < 1e-6:
            raise ValueError("The length of ridge %s is (nearly) zero." % name)
        cell[:,col] /= norm
    if abs(numpy.linalg.det(cell)) < 1e-6:
        raise ValueError("The ridges of the unit cell are (nearly) linearly dependent vectors.")


class UnitCell(object):
    def __init__(self):
        self.cell = numpy.array([
            [10.0,  0.0,  0.0], 
            [ 0.0, 10.0,  0.0], 
            [ 0.0,  0.0, 10.0]]
        )*angstrom,
        self.cell_active = numpy.array([False, False, False])
        self.update_reciproke()
        
    def set_cell(self, cell, norm_threshold=1e-6, volume_threshold=1e-6):
        check_cell(cell)
        self.cell = cell
        self.update_reciproke()
    
    def set_cell_active(self, cell_active):
        self.cell_active = cell_active
        self.update_reciproke()
    
    def get_active_inactive(self):
        active_indices = []
        inactive_indices = []
        for index, active in enumerate(self.cell_active):
            if active:
                active_indices.append(index)
            else:
                inactive_indices.append(index)
        return active_indices, inactive_indices
    
    def update_reciproke(self):
        active, inactive = self.get_active_inactive()
        if len(active) == 0:
            self.cell_reciproke = numpy.zeros((3, 3), float)
            return
        elif len(active) == 1:
            temp = copy.deepcopy(self.cell)
            temp[:, inactive[0]] = random_orthonormal(temp[:, active[0]])
            temp[:, inactive[1]] = numpy.cross(temp[:, inactive[0]], temp[:, active[0]])
        elif len(active) == 2:
            temp = copy.deepcopy(self.cell)
            temp[:, inactive[0]] = numpy.cross(temp[:, active[0]], temp[:, active[1]])
        elif len(active) == 3:
            temp = self.cell
        self.cell_reciproke = numpy.transpose(self.cell_active*numpy.transpose(numpy.linalg.inv(temp)))

    def to_fractional(self, coordinate):
        return numpy.dot(self.cell_reciproke, coordinate)
        
    def to_index(self, coordinate):
        return self.to_fractional(coordinate).round().astype(int)

    def to_coordinate(self, fractional):
        return numpy.dot(self.cell, fractional)

    def move_to_cell(self, delta, index):
        return delta + numpy.dot(self.cell, index)
        
    def shortest_vector(self, delta):
        return self.move_to_cell(delta, -self.to_index(delta))
        
    def add_cell_vector(self, vector, norm_threshold=1e-6, volume_threshold=1e-6):
        active, inactive = self.get_active_inactive()
        if len(active) == 3:
            raise ValueError("The unit cell already has three axes.")
        norm_vector = math.sqrt(numpy.dot(vector, vector))
        if norm_vector < norm_threshold:
            raise ValueError("The norm of the proposed vector must be significantly larger than zero.")
        if len(active) == 0:
            # Add the vector
            self.cell[:,0] = vector
            self.cell_active[0] = True
            # Make sure that the unused vectors are not linearly dependent
            normal = vector/norm_vector
            self.cell[:,1] = random_orthonormal(normal)
            self.cell[:,2] = numpy.cross(self.cell[:,0], self.cell[:,1])
            self.cell[:,1] *= length
            # update
            self.update_reciproke()
        a = self.cell[:,active[0]]
        norm_a = math.sqrt(numpy.dot(a, a))
        if len(active) == 1:
            if 1 - abs(numpy.dot(a, vector) / (norm_vector * norm_a)) < volume_threshold:
                raise ValueError("Can not add the vector since it is colinear with the existing unit cell axis.")
            else:
                # Add the vector
                self.cell[:,inactive[0]] = vector
                self.cell_active[inactive[0]] = True
                # Make sure that the unused vector is not linearly dependent
                self.cell[:,2] = numpy.cross(self.cell[:,0], self.cell[:,1])
                # update
                self.update_reciproke()
        b = self.cell[:,active[1]]
        norm_b = math.sqrt(numpy.dot(b, b))
        if len(active) == 2:
            backup = self.cell[:,inactive[0]].copy()
            self.cell[:,inactive[0]] = vector
            if abs(numpy.linalg.det(self.cell) / (norm_vector * norm_a * norm_b)) < volume_threshold:
                self.cell[:,inactive[0]] = backup
                raise ValueError("Can not add the vector since it is linearly dependent on the existing unit cell axes.")
            else:
                self.cell_active[inactive[0]] = True
                self.update_reciproke()
                return True

    def get_parameters(self):
        length_a = math.sqrt(numpy.dot(self.cell[:,0], self.cell[:,0]))
        length_b = math.sqrt(numpy.dot(self.cell[:,1], self.cell[:,1]))
        length_c = math.sqrt(numpy.dot(self.cell[:,2], self.cell[:,2]))
        alpha = math.acos(numpy.dot(self.cell[:,1], self.cell[:,2]) / (length_b * length_c))
        beta = math.acos(numpy.dot(self.cell[:,2], self.cell[:,0]) / (length_c * length_a))
        gamma = math.acos(numpy.dot(self.cell[:,0], self.cell[:,1]) / (length_a * length_b))
        return (numpy.array([length_a, length_b, length_c], float), numpy.array([alpha, beta, gamma], float))

    def generalized_volume(self):
        active, inactive = self.get_active_inactive()
        if len(active) == 0:
            return -1
        elif len(active) == 1:
            return math.sqrt(numpy.dot(self.cell[:,active[0]], self.cell[:,active[0]]))
        elif len(active) == 2:
            return abs(numpy.dot(self.cell[:,active[0]], self.cell[:,active[1]]))
        elif len(active) == 3:
            return abs(numpy.linalg.det(self.cell))
