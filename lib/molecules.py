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
from molmod.utils import cached, ReadOnly

from StringIO import StringIO

import numpy


__all__ = ["Molecule"]


class Molecule(ReadOnly):
    """Extensible datastructure for molecular geometries

       Most attributes of the molecule object are treated as constants. If you
       want to modify the molecular geometry, just create a modified molecule
       object. This facilitates the caching of derived quantities such as the
       distance matrices, while it imposes a cleaner coding style without
       a signifacant computational overhead.
    """

    def __init__(self, numbers, coordinates=None, title=None, masses=None, graph=None):
        """Initialize a Molecule object

           Mandatory arguments:
             numbers  --  numpy array (1D, N elements) with the atom numbers

           Optional keyword arguments:
             coordinates  --  numpy array (2D, Nx3 elements) Cartesian coordinates
             title  --  a string with the name of the molecule
             massess  --  a numpy array with atomic masses in atomic units
             graph  --  a MolecularGraph instance
        """
        ReadOnly.__init__(self)
        mandatory = {"numbers": numpy.array(numbers)}
        if coordinates is not None:
            coordinates = numpy.array(coordinates, float)
        if masses is not None:
            masses = numpy.array(masses, float)
        optional = {
            "coordinates": coordinates,
            "title": title,
            "masses": masses,
            "graph": graph,
        }
        self._init_attributes(mandatory, optional)

    size = property(lambda self: self.numbers.shape[0])

    @cached
    def distance_matrix(self):
        """The matrix with all atom pair distances"""
        from molmod.ext import molecules_distance_matrix
        return molecules_distance_matrix(self.coordinates)

    def dump_atoms(self, f):
        """Dump the Cartesian coordinates of the atoms to a file

           The format is the XYZ format without title and header line.
           Argument:
             f  --  a file-like object
        """
        for number, coordinate in zip(self.numbers, self.coordinates/angstrom):
            atom_info = periodic[number]
            if atom_info is None:
                symbol = "X"
            else:
                symbol = atom_info.symbol
            print >> f, "% 2s % 12.6f % 12.6f % 12.6f" % (
                symbol,
                coordinate[0],
                coordinate[1],
                coordinate[2]
            )

    def dumps_atoms(self):
        """Returns the Cartesian coordinates of the atoms as a string

           The format is the XYZ format without title and header line.
        """
        sio = StringIO()
        self.dump_atoms(sio)
        result = sio.getvalue()
        sio.close()
        return result

    def dump(self, f):
        """Dump the molecular geometry to an XYZ file

           Argument:
             f  --  a file-like object
        """
        print >> f, "%5i" % len(self.numbers)
        if hasattr(self, "title"):
            print >> f, str(self.title)
        else:
            print >> f
        self.dump_atoms(f)

    def write_to_file(self, filename):
        """Write the molecular geometry to an XYZ file

           Argument:
             filename  --  a filename
        """
        f = file(filename, 'w')
        self.dump(f)
        f.close()


