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


class FormattedCheckpoint(object):
    def __init__(self, filename):
        self.read_filename(filename)
        
    def read_filename(self, filename):
        f = file(filename, 'r')
        self.title = f.readline()[:-1]
        self.command, self.lot, self.basis = f.readline().split()
        self.fields = {}
        
        def read_field(f):
            line = f.readline()
            if line == "":
                return False
            
            label = line[:43].strip()
            line = line[43:]
            words = line.split()
            
            if words[0] == 'I':
                datatype = int
            elif words[0] == 'R':
                datatype = float
            else:
                raise ReadError("Unexpected datatype in formatted checkpoint file %s\n%s" % (filename, words[1])) 
            
            if len(words) == 2:
                value = datatype(words[1])
            elif len(words) == 3:
                if words[1] != "N=":
                    raise ReadError("Unexpected line in formatted checkpoint file %s\n%s" % (filename, line[:-1])) 
                length = int(words[2])
                value = numpy.zeros(length, datatype)
                counter = 0
                while counter < length:
                    line = f.readline()
                    if line == "":
                        raise ReadError("Unexpected end of formatted checkpoint file %s" % filename) 
                    for word in line.split():
                        value[counter] = datatype(word)
                        counter += 1
            else:
                raise ReadError("Unexpected line in formatted checkpoint file %s\n%s" % (filename, line[:-1])) 
            
            self.fields[label] = value
            return True
        
        while read_field(f):
            pass
            
        self.analyze()
        
    def analyze(self):
        print self.fields.keys()
        self.molecule = Molecule()
        self.molecule.numbers = self.fields["Atomic numbers"]
        self.molecule.coordinates = numpy.reshape(self.fields["Current cartesian coordinates"], (-1,3))
        
    def optimization_coordinates(self):
        geom_array = self.fields.get("Opt point       1 Geometries")
        assert geom_array is not None, "No optimization geometries found."
        return numpy.reshape(geom_array, (-1, len(self.molecule.numbers), 3))
