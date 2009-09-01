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


from molmod.io.common import slice_match, FileFormatError
from molmod.data.periodic import periodic
from molmod.molecules import Molecule
from molmod.units import angstrom

import numpy
from itertools import izip


__all__ = ["XYZReader", "XYZWriter", "XYZFile"]


class XYZReader(object):
    """A reader for XYZ trajectory files

       Use this reader as an iterator:
       >>> xr = XYZReader("somefile.xyz")
       >>> for title, coordinates in xr:
               print title
    """

    def __init__(self, f, sub=slice(None), file_unit=angstrom):
        """Initialize an XYZ reader

           Arguments:
             f  --  a filename or a file-like object
             sub  --  a slice indicating which frames to read/skip
             file_unit  --  the conversion constant to convert data into atomic
                            units (default=angstrom)

           After initialization, the following attributes are defined:
             self.symbols  --  The atom symbols
             self.numbers  --  The atom numbers
        """
        if isinstance(f, file):
            self._auto_close = False
            self._file = f
        else:
            self._auto_close = True
            self._file = file(f)
        self._sub = sub
        self._counter = 0
        self.file_unit = file_unit

        try:
            tmp = self._read()
            self.symbols = tmp[0]
            self._first = tmp[1:]
            self.numbers = numpy.zeros(len(self.symbols), int)
            for index, symbol in enumerate(self.symbols):
                atom_info = periodic[symbol]
                if atom_info is None:
                    self.numbers[index] = 0
                else:
                    self.numbers[index] = atom_info.number
        except StopIteration:
            raise FileFormatError("Could not read first frame from XYZ file. Incorrect file format.")

    def __del__(self):
        if self._auto_close:
            self._file.close()

    def _read(self):
        def read_size():
            try:
                return int(self._file.readline().strip())
            except ValueError:
                raise StopIteration

        while not slice_match(self._sub, self._counter):
            size = read_size()
            for i in xrange(size+1):
                line = self._file.readline()
                if len(line) == 0:
                    raise StopIteration
            self._counter += 1

        size = read_size()
        title = self._file.readline()[:-1]
        symbols = []
        coordinates = numpy.zeros((size, 3), float)
        for counter in xrange(size):
            line = self._file.readline()
            if len(line) == 0:
                raise StopIteration
            words = line.split()
            if len(words) < 4:
                raise StopIteration
            symbols.append(words[0])
            try:
                coordinates[counter,0] = float(words[1])
                coordinates[counter,1] = float(words[2])
                coordinates[counter,2] = float(words[3])
            except ValueError:
                raise StopIteration
        coordinates *= self.file_unit
        self._counter += 1
        return symbols, title, coordinates

    def __iter__(self):
        return self

    def next(self):
        """Read the next frame from the XYZ trajectory file

           This method is part of the iterator protocol
        """
        if self._first is not None:
            tmp = self._first
            self._first = None
            result = tmp
        else:
            result = self._read()[1:]
        if result is None:
            raise StopIteration
        return result

    def get_first_molecule(self):
        """Get the first molecule from the trajectory

           This can be useful to configure your program before handeling the
           actual trajectory.
        """
        if self._first is None:
            raise RuntimeError("get_first_molecule must be called before the first iteration.")
        else:
            title, coordinates = self._first
            molecule = Molecule(self.numbers, coordinates, title)
            return molecule


class XYZWriter(object):
    """A XYZ trajectory file writer

       This writer is designed to be used in conjunction with an iterator that
       generates coordinates in one way or the other. Example:

       xr = XYZReader("somefile.xyz")
       xw = XYZWriter("otherfile.xyz", xr.symbols[5:10])
       for title, coordinates in xr:
           xw.dump(title, -coordinates[5:10])
    """
    def __init__(self, f, symbols, file_unit=angstrom):
        """Initialize a XYZWriter object

           Arguments:
             f  -- a filename or a file-like object to write to
             symbols  --  the atom symbols
             file_unit  --  the unit of the values written to file
                            (default=angstrom)
        """
        if isinstance(f, file):
            self._auto_close = False
            self._file = f
        else:
            self._auto_close = True
            self._file = file(f, 'w')
        self.symbols = symbols
        self.file_unit = file_unit

    def __del__(self):
        if self._auto_close:
            self._file.close()

    def dump(self, title, coordinates):
        """Dump a frame to the trajectory file

           Arguments:
             title  --  the title of the frame
             coordinates  --  a numpy array with coordinates in atomic units
        """
        print >> self._file, "% 8i" % len(self.symbols)
        print >> self._file, str(title)
        for symbol, coordinate in zip(self.symbols, coordinates):
            print >> self._file, "% 2s % 12.9f % 12.9f % 12.9f" % ((symbol,) + tuple(coordinate/self.file_unit))


class XYZFile(object):
    """Data structure representing an XYZ File

       This implementation is extra layer on top the XYZReader and XYZWriter,
       and is somewhat easier to use. Example:

       xyz_file = XYZFile("some.xyz")
       mol = xyz_file.get_molecule(3)
       print mol.title
       xyz_file.geometries[0,4,2] = 5.0 # frame 0, atom 4, Z-coordinate
       xyz_file.write_to_file("other.xyz")
    """
    def __init__(self, f, sub=slice(None), file_unit=angstrom):
        """Initialize an XYZFile object

           Argument:
             f  --  a filename or a file-like object to read from
             sub  --  a slice indicating which frames to read/skip
             file_unit  --  the conversion constant to convert data into atomic
                            units (default=angstrom)

           XYZFile instances always have to following attriubtes:
             numbers  --  The atom numbers (of one frame)
             symbols  --  The atom symbols (of one frame)
             titles   --  The titles of all the frames
             geometries  --  A MxNx3 array with all the atom coordinates
        """
        xyz_reader = XYZReader(f, file_unit=file_unit)
        self.file_unit = file_unit

        self.numbers = xyz_reader.numbers
        self.symbols = xyz_reader.symbols
        self.titles = []
        if sub.stop is not None:
            start = sub.start
            if start is None: start = 0
            stride = sub.stride
            if stride is None: stride = 1
            count = (sub.stop - start)/stride
            self.geometries = numpy.zeros((count, len(self.numbers), 3), float)

            for counter, (title, coordinates) in xyz_reader:
                self.titles.append(title)
                self.geometries[counter] = coordinates
        else:
            geometries = []
            for title, coordinates in xyz_reader:
                self.titles.append(title)
                geometries.append(coordinates)
            self.geometries = numpy.array(geometries, float)

    def get_molecule(self, index=0):
        """Get a molecule from the trajectory

           Argument:
             index  --  The frame index (default=0)
        """
        return Molecule(self.numbers, self.geometries[index], self.titles[index])

    def write_to_file(self, f, file_unit=angstrom):
        """Write the trajectory to a file

           Arguments:
             f  -- a filename or a file-like object to write to
             file_unit  --  the unit of the values written to file
                            (default=angstrom)
        """
        xyz_writer = XYZWriter(f, self.symbols, file_unit=file_unit)
        for title, coordinates in izip(self.titles, self.geometries):
            xyz_writer.dump(title, coordinates)

