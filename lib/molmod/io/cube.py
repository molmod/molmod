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


__all__ = ["CubeReader"]


class CubeReader(object):
    def __init__(self, filename):
        self.f = file(filename)
        # skip the first two lines
        self.f.readline()
        self.f.readline()

        def read_header_line(line):
          words = line.split()
          return (
              int(words[0]),
              numpy.array([float(words[1]), float(words[2]), float(words[3])], float)
          )

        # number of atoms and origin of the grid
        self.num_atoms, self.origin = read_header_line(self.f.readline())
        # numer of grid points in A direction and step vector A, and so on
        self.num_a, self.vector_a = read_header_line(self.f.readline())
        self.num_b, self.vector_b = read_header_line(self.f.readline())
        self.num_c, self.vector_c = read_header_line(self.f.readline())

        def read_header_line(line):
          words = line.split()
          return (
              int(words[0]),
              numpy.array([float(words[2]), float(words[3]), float(words[4])], float)
          )

        self.numbers = numpy.zeros(self.num_atoms, int)
        self.coordinates = numpy.zeros((self.num_atoms, 3), float)
        for i in xrange(self.num_atoms):
            self.numbers[i], self.coordinates[i] = read_header_line(self.f.readline())

        self.counter_a = 0
        self.counter_b = 0
        self.counter_c = 0
        self.values = []

    def __del__(self):
        self.f.close()

    def __iter__(self):
        return self

    def next(self):
        if len(self.values) == 0:
            line = self.f.readline()
            if len(line) == 0:
                raise StopIteration
            self.values = [float(word) for word in line.split()]
        value = self.values.pop(0)
        vector = self.origin + self.counter_a*self.vector_a + self.counter_b*self.vector_b + self.counter_c*self.vector_c
        self.counter_c += 1
        if self.counter_c >= self.num_c:
            self.counter_c = 0
            self.counter_b += 1
            if self.counter_b >= self.num_b:
                self.counter_b = 0
                self.counter_a += 1
                if self.counter_a >= self.num_a:
                    raise StopIteration
        return value, vector


