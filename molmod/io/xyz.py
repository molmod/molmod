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
"""Tools for reading and writing XYZ trajectory files"""


from __future__ import print_function


from __future__ import division

from builtins import range
import numpy as np

from molmod.io.common import SlicedReader, FileFormatError
from molmod.periodic import periodic
from molmod.molecules import Molecule
from molmod.units import angstrom


__all__ = ["XYZReader", "XYZWriter", "XYZFile"]


class XYZReader(SlicedReader):
    """A reader for XYZ trajectory files

       Use this reader as an iterator::

         >>> xr = XYZReader("somefile.xyz")
         >>> for title, coordinates in xr:
                 print title
    """

    def __init__(self, f, sub=slice(None), file_unit=angstrom):
        """Initialize an XYZ reader

           Arguments:
            | ``f``  --  a filename or a file-like object

           Optional arguments:
            | ``sub``  --  a slice indicating which frames to read/skip
            | ``file_unit``  --  the conversion constant to convert data into atomic
                                  units [default=angstrom]

           After initialization, the following attributes are defined:
            | ``symbols``  --  The atom symbols
            | ``numbers``  --  The atom numbers
        """
        SlicedReader.__init__(self, f, sub)
        self.file_unit = file_unit

        try:
            self.symbols = None
            self._first = self._read_frame()
            self.numbers = np.zeros(len(self.symbols), int)
            for index, symbol in enumerate(self.symbols):
                try:
                    number = int(symbol)
                    symbol = periodic[number].symbol
                    self.symbols[index] = symbol
                except ValueError:
                    atom_info = periodic[symbol]
                    if atom_info is not None:
                        number = atom_info.number
                    else:
                        number = 0
                self.numbers[index] = number
            self.symbols = tuple(self.symbols)
            self._f.seek(0)
        except StopIteration:
            raise FileFormatError("Could not read first frame from XYZ file. Incorrect file format.")

    def read_size(self):
        """Read the number of atoms"""
        try:
            return int(self._f.readline().strip())
        except ValueError:
            raise StopIteration

    def _read_frame(self):
        """Read a frame from the XYZ file"""

        size = self.read_size()
        title = self._f.readline()[:-1]
        if self.symbols is None:
            symbols = []
        coordinates = np.zeros((size, 3), float)
        for counter in range(size):
            line = self._f.readline()
            if len(line) == 0:
                raise StopIteration
            words = line.split()
            if len(words) < 4:
                raise StopIteration
            if self.symbols is None:
                symbols.append(words[0])
            try:
                coordinates[counter, 0] = float(words[1])
                coordinates[counter, 1] = float(words[2])
                coordinates[counter, 2] = float(words[3])
            except ValueError:
                raise StopIteration
        coordinates *= self.file_unit
        if self.symbols is None:
            self.symbols = symbols
        return title, coordinates

    def _skip_frame(self):
        """Skip a single frame from the trajectory"""
        size = self.read_size()
        for i in range(size+1):
            line = self._f.readline()
            if len(line) == 0:
                raise StopIteration

    def get_first_molecule(self):
        """Get the first molecule from the trajectory

           This can be useful to configure your program before handeling the
           actual trajectory.
        """
        title, coordinates = self._first
        molecule = Molecule(self.numbers, coordinates, title, symbols=self.symbols)
        return molecule


class XYZWriter(object):
    """A XYZ trajectory file writer

       This writer is designed to be used in conjunction with an iterator that
       generates coordinates in one way or the other. Example::

         >>> xr = XYZReader("somefile.xyz")
         >>> xw = XYZWriter("otherfile.xyz", xr.symbols[5:10])
         >>> for title, coordinates in xr:
         ...    xw.dump(title, -coordinates[5:10])
    """
    def __init__(self, f, symbols, file_unit=angstrom):
        """
           Arguments:
            | ``f``  -- a filename or a file-like object to write to
            | ``symbols``  --  the atom symbols

           Optional argument
            | ``file_unit``  --  the unit of the values written to file
                                 [default=angstrom]
        """
        if isinstance(f, str):
            self._auto_close = True
            self._f = open(f, 'w')
        else:
            self._auto_close = False
            self._f = f
        self.symbols = symbols
        self.file_unit = file_unit

    def __del__(self):
        if self._auto_close:
            self._f.close()

    def dump(self, title, coordinates):
        """Dump a frame to the trajectory file

           Arguments:
            | ``title``  --  the title of the frame
            | ``coordinates``  --  a numpy array with coordinates in atomic units
        """
        print("% 8i" % len(self.symbols), file=self._f)
        print(str(title), file=self._f)
        for symbol, coordinate in zip(self.symbols, coordinates):
            print("% 2s % 12.9f % 12.9f % 12.9f" % ((symbol, ) + tuple(coordinate/self.file_unit)), file=self._f)


class XYZFile(object):
    """Data structure representing an XYZ File

       This implementation is extra layer on top the XYZReader and XYZWriter,
       and is somewhat easier to use. Example::

         >>> xyz_file = XYZFile("some.xyz")
         >>> mol = xyz_file.get_molecule(3)
         >>> print mol.title
         >>> xyz_file.geometries[0, 4, 2] = 5.0 # frame 0, atom 4, Z-coordinate
         >>> xyz_file.write_to_file("other.xyz")
    """
    def __init__(self, f, sub=slice(None), file_unit=angstrom):
        """Initialize an XYZFile object

           Argument:
            | ``f``  --  a filename or a file-like object to read from

           Optional arguments:
            | ``sub``  --  a slice indicating which frames to read/skip
            | ``file_unit``  --  the conversion constant to convert data into atomic
                                 units [default=angstrom]

           XYZFile instances always have to following attriubtes:
            | ``numbers``  --  The atom numbers (of one frame)
            | ``symbols``  --  The atom symbols (of one frame)
            | ``titles``   --  The titles of all the frames
            | ``geometries``  --  A MxNx3 array with all the atom coordinates
        """
        xyz_reader = XYZReader(f, sub, file_unit=file_unit)
        self.file_unit = file_unit

        self.numbers = xyz_reader.numbers
        self.symbols = xyz_reader.symbols
        self.titles = []
        if sub.stop is not None:
            start = sub.start
            if start is None: start = 0
            step = sub.step
            if step is None: step = 1
            count = (sub.stop - start + step - 1)//step
            self.geometries = np.zeros((count, len(self.numbers), 3), float)

            for counter, (title, coordinates) in enumerate(xyz_reader):
                self.titles.append(title)
                self.geometries[counter] = coordinates
            assert counter+1 == count
        else:
            geometries = []
            for title, coordinates in xyz_reader:
                self.titles.append(title)
                geometries.append(coordinates)
            self.geometries = np.array(geometries, float)

    def get_molecule(self, index=0):
        """Get a molecule from the trajectory

           Optional argument:
            | ``index``  --  The frame index [default=0]
        """
        return Molecule(self.numbers, self.geometries[index], self.titles[index], symbols=self.symbols)

    def write_to_file(self, f, file_unit=angstrom):
        """Write the trajectory to a file

           Argument:
            | ``f``  -- a filename or a file-like object to write to

           Optional argument:
            | ``file_unit``  --  the unit of the values written to file
                                 [default=angstrom]
        """
        xyz_writer = XYZWriter(f, self.symbols, file_unit=file_unit)
        for title, coordinates in zip(self.titles, self.geometries):
            xyz_writer.dump(title, coordinates)
