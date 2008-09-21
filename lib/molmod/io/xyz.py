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


from molmod.io.common import slice_match
from molmod.data.periodic import periodic
from molmod.molecules import Molecule
from molmod.units import angstrom

import numpy


__all__ = ["Error", "XYZReader", "XYZWriter", "XYZFile"]


class Error(Exception):
    pass


class XYZReader(object):
    def __init__(self, f, sub=slice(None), file_unit=angstrom):
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
            raise Error("Could not read first frame from XYZ file. Incorrect file format.")

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
        if self._first is None:
            raise Error("get_first_molecule must be called before the first iteration.")
        else:
            molecule = Molecule()
            molecule.numbers = self.numbers.copy()
            molecule.title, molecule.coordinates = self._first
            return molecule



class XYZWriter(object):
    def __init__(self, f, symbols, file_unit=angstrom):
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
        print >> self._file, "% 8i" % len(self.symbols)
        print >> self._file, str(title)
        for symbol, coordinate in zip(self.symbols, coordinates):
            print >> self._file, "% 2s % 12.9f % 12.9f % 12.9f" % ((symbol,) + tuple(coordinate/self.file_unit))


class XYZFile(object):
    def __init__(self, f, sub=slice(None), file_unit=angstrom):
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
        result = Molecule()
        result.numbers = self.numbers
        result.coordinates = self.geometries[index]
        result.title = self.titles[index]
        return result

    def write_to_file(self, f, file_unit=angstrom):
        symbols = []
        for number in self.numbers:
            atom_info = periodic[number]
            if atom_info is None:
                symbols.append("X")
            else:
                symbols.append(atom_info.symbol)
        xyz_writer = XYZWriter(f, symbols, file_unit=file_unit)
        for index, coordinates in enumerate(self.geometries):
            xyz_writer.dump("Step %i" % index, coordinates)

