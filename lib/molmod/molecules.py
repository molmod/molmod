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

from molmod.data.periodic import periodic
from molmod.units import angstrom
from molmod.graphs import cached

from StringIO import StringIO

import numpy


__all__ = ["Molecule"]


class Molecule(object):
    def __init__(self, numbers=None, coordinates=None, title=None):
        # a lot of fuzz to make the molecule a read only object.
        self._numbers = numpy.array(numbers, int)
        self._numbers.setflags(write=False)
        self._coordinates = numpy.array(coordinates, float)
        self._coordinates.setflags(write=False)
        self.title = title

    size = property(lambda self: self._numbers.shape[0])
    numbers = property(lambda self: self._numbers)
    coordinates = property(lambda self: self._coordinates)

    def dump_atoms(self, stream):
        for number, coordinate in zip(self.numbers, self.coordinates/angstrom):
            atom_info = periodic[number]
            if atom_info is None:
                symbol = "X"
            else:
                symbol = atom_info.symbol
            print >> stream, "% 2s % 12.6f % 12.6f % 12.6f" % (
                symbol,
                coordinate[0],
                coordinate[1],
                coordinate[2]
            )

    def dumps_atoms(self):
        sio = StringIO()
        self.dump_atoms(sio)
        result = sio.getvalue()
        sio.close()
        return result

    def dump(self, stream):
        print >> stream, "%5i" % len(self.numbers)
        print >> stream, str(self.title)
        self.dump_atoms(stream)

    def write_to_file(self, filename):
        f = file(filename, 'w')
        self.dump(f)
        f.close()

    @cached
    def distance_matrix(self):
        from molmod.ext import molecules_distance_matrix
        return molecules_distance_matrix(self.coordinates)


