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


import numpy, math

class Transformation(object):
    def __init__(self):
        self.a = numpy.zeros((3, 3), float)
        self.b = numpy.zeros(3, float)
        
    def apply(self, vector):
        return numpy.dot(self.a, vector) + self.b


class Rotation(Transformation):
    def __init__(self, angle, axis, center=numpy.zeros(3, float)):
        Transformation.__init__(self)
        # Rodriguez
        norm = math.sqrt(numpy.dot(axis, axis))
        assert norm > 0
        x = axis[0] / norm
        y = axis[1] / norm
        z = axis[2] / norm
        c = math.cos(angle)
        s = math.sin(angle)
        self.a = numpy.array([
            [x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s],
            [x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s],
            [x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  ]
        ])
        
        self.b = center - numpy.dot(self.a, center)
