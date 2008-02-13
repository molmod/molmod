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


from molmod.molecules import Molecule
from molmod.units import angstrom
from molmod.grids import Grid

import numpy


class ReadError(Exception):
    pass


class Cube(Grid):
    def __init__(self, filename, label, grid_conversion=1, volume_filename=None):
        f = file(filename, 'r')

        # 1) skip the first two lines
        self.label1 = f.readline()[:-1]
        if self.label1 == "": raise ReadError("Unexpected EOF, expected first comment line.")
        self.label2 = f.readline()[:-1]
        if self.label2 == "": raise ReadError("Unexpected EOF, Expected second comment line.")

        # 2) read grid specification

        # 2a) origin line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with origin.")
        words = line.split()
        num_atoms = int(words[0])
        origin = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)*grid_conversion

        # 2b) x line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with x-axis.")
        words = line.split()
        nx = int(words[0])
        x = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)*grid_conversion

        # 2c) y line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with y-axis.")
        words = line.split()
        ny = int(words[0])
        y = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)*grid_conversion

        # 2d) z line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with z-axis.")
        words = line.split()
        nz = int(words[0])
        z = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)*grid_conversion

        # 3) read atoms xyz data
        self.molecule = Molecule()
        self.molecule.numbers = numpy.zeros(num_atoms, int)
        self.molecule.coordinates = numpy.zeros((num_atoms,3), float)
        for index in xrange(num_atoms):
            line = f.readline()
            if line == "": raise ReadError("Unexpected EOF, expected line with atom (%i)." % index)
            words = line.split()
            self.molecule.numbers[index] = int(words[0])
            self.molecule.coordinates[index,0] = float(words[2])
            self.molecule.coordinates[index,1] = float(words[3])
            self.molecule.coordinates[index,2] = float(words[4])

        user_points = (sum(origin**2) + sum(x**2) + sum(y**2) + sum(z**2) == 0.0)

        if user_points:
            # 4a) read user points
            points = numpy.zeros((nx, 3), float)
            data = numpy.zeros(nx, float)
            for index in xrange(nx):
                line = f.readline()
                if line == "": raise ReadError("Unexpected EOF, expected line with grid data (index=%i)" % index)
                words = line.split()
                points[index] = [float(word) for word in words[:3]]
                points[index] *= angstrom
                data[index] = float(words[3])
        else:
            # 4) read grid values
            n = nx*ny*nz
            points = numpy.zeros((n, 3), float)
            data = numpy.zeros(n, float)
            buffer = []
            counter = 0
            for i in xrange(nx):
                for j in xrange(ny):
                    for k in xrange(nz):
                        if len(buffer) == 0:
                            line = f.readline()
                            if line == "": raise ReadError("Unexpected EOF, expected line with grid data (i=%i,j=%i,k=%i)." % (i,j,k))
                            buffer = [float(word) for word in line.split()]
                        points[counter] = origin + i*x + j*y + k*z
                        data[counter] = buffer.pop(0)
                        counter += 1

        f.close()

        if volume_filename is not None:
            volumes = numpy.zeros(len(points), float)
            f = file(volume_filename)
            for index, line in enumerate(f):
                volumes[index] = float(line)
            f.close()
        else:
            if user_points:
                volumes = None
            else:
                volumes = numpy.zeros(len(points), float)
                volumes += abs(numpy.linalg.det(numpy.array([x, y, z], float)))

        Grid.__init__(self, points, volumes, {label: data})





