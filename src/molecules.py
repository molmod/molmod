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

from pychem.moldata import periodic
from pychem.units import from_angstrom

from StringIO import StringIO

import Numeric, math


__all__ = [
    "Molecule", 
    "molecule_from_xyz_stream", "molecule_from_xyz_str", 
    "molecule_from_xyz_filename"
]


class Molecule:
    """
    A Molecule instance describes a molecule in the following representation:
    - cartesian coordinates
    """

    def __init__(self, atoms=[]):
        """
        Initialiaze a Molecule instance.
        
        arguments:
        atoms -- [[number, x, y, z], ...]
        """
        self.numbers = Numeric.zeros(len(atoms), Numeric.Int)
        self.coordinates = Numeric.zeros((len(atoms), 3), Numeric.Float)
        for index, line in enumerate(atoms):
            self.numbers[index] = line[0]
            self.coordinates[index] = line[1:4]
                
    def normalize(self):
        """
        Bring the molecule in a normalized frame. This only works if the
        first three atoms are not colinear which is the case for general
        molecules.
        """
        # first translate the first atom to the center
        self.coordinates -= self.coordinates[0].copy()
        # then rotate the molecule so that the second atom lies on the positive x-axis
        # and the third atom lies in the xy-plane with positive y.
        new_x = self.coordinates[1].copy()
        new_x /= math.sqrt(Numeric.dot(new_x, new_x))
        third = self.coordinates[2].copy()
        new_z = Numeric.array([
            new_x[1]*third[2]-third[1]*new_x[2],
            new_x[2]*third[0]-third[2]*new_x[0],
            new_x[0]*third[1]-third[0]*new_x[1]
        ])
        new_z /= math.sqrt(Numeric.dot(new_z, new_z))
        new_y = Numeric.array([
            new_z[1]*new_x[2]-new_x[1]*new_z[2],
            new_z[2]*new_x[0]-new_x[2]*new_z[0],
            new_z[0]*new_x[1]-new_x[0]*new_z[1]
        ])
        rotation = Numeric.transpose(Numeric.array([new_x, new_y, new_z]))
        self.coordinates = Numeric.dot(self.coordinates, rotation)


def molecule_from_xyz_stream(stream):
    num_atoms = None
    atoms = []
    for line in stream:
        line = line[0:line.find("#")]
        words = line.split()
        if len(words) == 1 and num_atoms == None:
            num_atoms = int(words[0])
        elif len(words) == 4:
            atoms.append([
                periodic.symbol_lookup(words[0]), 
                from_angstrom(float(words[1])), 
                from_angstrom(float(words[2])),
                from_angstrom(float(words[3]))
            ])
    assert len(atoms) == num_atoms, "Inconsistent number of atoms."
    return Molecule(atoms)

def molecule_from_xyz_filename(filename):
    """Load an xyz file and return a Molecule instance."""
    f = file(filename)
    result = molecule_from_xyz_stream(f)
    f.close()
    return result

def molecule_from_xyz_string(string):
    stream = StringIO(string)
    result = molecule_from_xyz_stream(stream)
    stream.close()
    return result

