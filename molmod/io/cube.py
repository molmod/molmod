# -*- coding: utf-8 -*-
# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
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
"""Reader for the cube format"""


from __future__ import print_function

from builtins import range, object
import numpy as np

from molmod.molecules import Molecule


__all__ = ['get_cube_points', 'CubeReader', 'Cube']


def get_cube_points(origin, axes, nrep):
    '''Generate the Cartesian coordinates of the points in a cube file

       *Arguemnts:*

       origin
            The cartesian coordinate for the origin of the grid.

       axes
            The 3 by 3 array with the grid spacings as rows.

       nrep
            The number of grid points along each axis.
    '''
    points = np.zeros((nrep[0], nrep[1], nrep[2], 3), float)
    points[:] = origin
    points[:] += np.outer(np.arange(nrep[0], dtype=float), axes[0]).reshape((-1,1,1,3))
    points[:] += np.outer(np.arange(nrep[1], dtype=float), axes[1]).reshape((1,-1,1,3))
    points[:] += np.outer(np.arange(nrep[2], dtype=float), axes[2]).reshape((1,1,-1,3))
    return points


def read_cube_header(f):
    # skip the first two lines
    title = f.readline().strip()
    subtitle = f.readline().strip()

    def read_grid_line(line):
        """Read a grid line from the cube file"""
        words = line.split()
        return (
            int(words[0]),
            np.array([float(words[1]), float(words[2]), float(words[3])], float)
            # all coordinates in a cube file are in atomic units
        )

    # number of atoms and origin of the grid
    natom, origin = read_grid_line(f.readline())
    # numer of grid points in A direction and step vector A, and so on
    nrep0, axis0 = read_grid_line(f.readline())
    nrep1, axis1 = read_grid_line(f.readline())
    nrep2, axis2 = read_grid_line(f.readline())
    nrep = np.array([nrep0, nrep1, nrep2], int)
    axes = np.array([axis0, axis1, axis2])

    def read_coordinate_line(line):
        """Read an atom number and coordinate from the cube file"""
        words = line.split()
        return (
            int(words[0]), float(words[1]),
            np.array([float(words[2]), float(words[3]), float(words[4])], float)
            # all coordinates in a cube file are in atomic units
        )

    numbers = np.zeros(natom, int)
    nuclear_charges = np.zeros(natom, float)
    coordinates = np.zeros((natom, 3), float)
    for i in range(natom):
        numbers[i], nuclear_charges[i], coordinates[i] = read_coordinate_line(f.readline())

    molecule = Molecule(numbers, coordinates, title=title)
    return molecule, origin, axes, nrep, subtitle, nuclear_charges


class CubeReader(object):
    """Iterator that reads cube files. See the cubegen manual for more
       information about cube files:
       http://www.gaussian.com/g_tech/g_ur/u_cubegen.htm

       Use this object as an iterator::

         >>> cr = CubeReader("test.cube")
         >>> print cr.numbers
         >>> for vector, value in cr:
         ...     print vector, value
    """
    def __init__(self, filename):
        """
           Argument:
            | ``filename``  --  the filename with the formatted cube data
        """
        self.f = open(filename)

        self.molecule, self.origin, self.axes, self.nrep, self.subtitle, \
            self.nuclear_charges = read_cube_header(self.f)

        self._counter0 = 0
        self._counter1 = 0
        self._counter2 = 0
        self._done = False
        self._values = []

    def __del__(self):
        self.f.close()

    def __iter__(self):
        return self

    def __next__(self):
        """Read the next datapoint from the cube file

           This method is part of the iterator protocol.
        """
        if self._done:
            raise StopIteration
        if len(self._values) == 0:
            line = self.f.readline()
            if len(line) == 0:
                raise StopIteration
            self._values = [float(word) for word in line.split()]
        value = self._values.pop(0)
        vector = self.origin + self._counter0*self.axes[0] \
                             + self._counter1*self.axes[1] \
                             + self._counter2*self.axes[2]
        self._counter2 += 1
        if self._counter2 >= self.nrep[2]:
            self._counter2 = 0
            self._counter1 += 1
            if self._counter1 >= self.nrep[1]:
                self._counter1 = 0
                self._counter0 += 1
                if self._counter0 >= self.nrep[0]:
                    self._done = True
        return vector, value


class Cube(object):
    '''A data structure for cube file data.
    '''
    @classmethod
    def from_file(cls, filename):
        '''Create a cube object by loading data from a file.

           *Arguemnts:*

           filename
                The file to load. It must contain the header with the
                description of the grid and the molecule.
        '''
        with open(filename) as f:
            molecule, origin, axes, nrep, subtitle, nuclear_charges = \
                read_cube_header(f)
            data = np.zeros(tuple(nrep), float)
            tmp = data.ravel()
            counter = 0
            while True:
                line = f.readline()
                if len(line) == 0:
                    break
                words = line.split()
                for word in words:
                    tmp[counter] = float(word)
                    counter += 1
        return cls(molecule, origin, axes, nrep, data, subtitle, nuclear_charges)

    def __init__(self, molecule, origin, axes, nrep, data, subtitle='', nuclear_charges=None):
        '''
           *Arguments:*

           molecule
                A Molecule instance.

           origin
                The cartesian coordinate for the origin of the grid.

           axes
                The 3 by 3 array with the grid spacings as rows.

           nrep
                The number of grid points along each axis.

           subtitle
                The title on the second line in the cube file.

           nuclear_charges
                The nuclear charges, can be different from the atomic numbers
                in case of effective core potentials.
        '''
        self.molecule = molecule
        self.origin = np.array(origin, copy=False)
        self.axes = np.array(axes, copy=False)
        self.nrep = np.array(nrep, copy=False)
        self.data = np.array(data, copy=False)
        self.subtitle = subtitle
        if nuclear_charges is None:
            self.nuclear_charges = self.molecule.numbers.astype(float)
        else:
            self.nuclear_charges = nuclear_charges

    def write_to_file(self, fn):
        '''Write the cube to a file in the Gaussian cube format.'''
        with open(fn, 'w') as f:
            f.write(' {}\n'.format(self.molecule.title))
            f.write(' {}\n'.format(self.subtitle))

            def write_grid_line(n, v):
                f.write('%5i % 11.6f % 11.6f % 11.6f\n' % (n, v[0], v[1], v[2]))

            write_grid_line(self.molecule.size, self.origin)
            write_grid_line(self.data.shape[0], self.axes[0])
            write_grid_line(self.data.shape[1], self.axes[1])
            write_grid_line(self.data.shape[2], self.axes[2])

            def write_atom_line(n, nc, v):
                f.write('%5i % 11.6f % 11.6f % 11.6f % 11.6f\n' % (n, nc, v[0], v[1], v[2]))

            for i in range(self.molecule.size):
                write_atom_line(self.molecule.numbers[i], self.nuclear_charges[i],
                                self.molecule.coordinates[i])

            for i0 in range(self.data.shape[0]):
                for i1 in range(self.data.shape[1]):
                    col = 0
                    for i2 in range(self.data.shape[2]):
                        value = self.data[i0, i1, i2]
                        if col % 6 == 5:
                            f.write(' % 12.5e\n' % value)
                        else:
                            f.write(' % 12.5e' % value)
                        col += 1
                    if col % 6 != 5:
                        f.write('\n')

    def copy(self, newdata=None):
        '''Return a copy of the cube with optionally new data.'''
        if newdata is None:
            newdata = self.data.copy()
        return self.__class__(
            self.molecule, self.origin.copy(), self.axes.copy(),
            self.nrep.copy(), newdata, self.subtitle, self.nuclear_charges
        )

    def get_points(self):
        '''Return a Nz*Nb*Nc*3 array with all Cartesian coordinates of the
           points in the cube.
        '''
        return get_cube_points(self.origin, self.axes, self.nrep)
