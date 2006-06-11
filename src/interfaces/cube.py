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


import numpy

from molmod.molecules import Molecule


class ReadError(Exception):
    pass


class Cube(object):
    def __init__(self, filename, grid_conversion=1):
        self.load_from_file(filename, grid_conversion)

    def load_from_file(self, filename, grid_conversion):
        f = file(filename, 'r')
        
        # 1) skip the first two lines
        self.label1 = f.readline()
        if self.label1 == "": raise ReadError("Unexpected EOF, expected first comment line.")
        self.label2 = f.readline()
        if self.label2 == "": raise ReadError("Unexpected EOF, Expected second comment line.")

        # 2) read grid specification
        
        # 2a) origin line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with origin.")
        words = line.split()
        num_atoms = int(words[0])
        self.origin = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)*grid_conversion
        
        # 2b) x line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with x-axis.")
        words = line.split()
        self.nx = int(words[0])
        self.x = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)
        
        # 2c) y line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with y-axis.")
        words = line.split()
        self.ny = int(words[0])
        self.y = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)
        
        # 2d) z line
        line = f.readline()
        if line == "": raise ReadError("Unexpected EOF, expected line with z-axis.")
        words = line.split()
        self.nz = int(words[0])
        self.z = numpy.array([float(words[1]), float(words[2]), float(words[3])], float)
        
        # 3) read atoms xyz data
        self.molecule = Molecule()
        self.molecule.numbers = numpy.zeros(num_atoms, int)
        self.molecule.coordinates = numpy.zeros((num_atoms,3), float)
        for index in xrange(num_atoms):
            line = f.readline()
            if line == "": raise ReadError("Unexpected EOF, expected line with atom (%i)." % index)
            words = line.split()
            self.molecule.numbers[index] = int(words[0])
            self.molecule.coordinates[index,0] = float(words[1])
            self.molecule.coordinates[index,1] = float(words[2])
            self.molecule.coordinates[index,2] = float(words[3])
            
        # 4) read grid values
        self.grid_data = numpy.zeros((self.nx, self.ny, self.nz), float)
        buffer = []
        for i in xrange(self.nx):
            for j in xrange(self.ny):
                for k in xrange(self.nz):
                    if len(buffer) == 0:
                        line = f.readline()
                        if line == "": raise ReadError("Unexpected EOF, expected line with grid data (i=%i,j=%i,k=%i)." % (i,j,k))
                        buffer = [float(word) for word in line.split()]
                    self.grid_data[i, j, k] = buffer.pop(0)
        
        f.close()

    def points_and_values(self):
        result = numpy.zeros((self.nx*self.ny*self.nz, 4), float)
        counter = 0
        for i in xrange(self.nx):
            for j in xrange(self.ny):
                for k in xrange(self.nz):
                    result[counter,0:3] = self.origin + i*self.x + j*self.y + k*self.z
                    result[counter,3] = self.grid_data[i, j, k]
                    counter += 1
        return result
        
    
