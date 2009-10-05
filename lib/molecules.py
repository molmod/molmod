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
"""Data structure for molecular geometries"""


from molmod.periodic import periodic
from molmod.units import angstrom
from molmod.utils import cached, ReadOnly
from molmod.molecular_graphs import MolecularGraph

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
        mandatory = {"numbers": numpy.array(numbers, int)}
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

    @staticmethod
    def from_file(filename):
        """Construct a molecule object based read from the given file

           The file format is inferred from the extensions. Currently supported
           formats are: *.cml, *.fchk, *.pdb, *.sdf, *.xyz

           If a file contains more than one molecule, only the first one is
           read.

           Argument:
             filename  --  the name of the file containing the molecule
        """
        # TODO: many different API's to load files. brrr...
        if filename.endswith(".cml"):
            from molmod.io import load_cml
            return load_cml(filename)[0]
        elif filename.endswith(".fchk"):
            from molmod.io import FCHKFile
            fchk = FCHKFile(filename, field_labels=[])
            return fchk.molecule
        elif filename.endswith(".pdb"):
            from molmod.io import load_pdb
            return load_pdb(filename)
        elif filename.endswith(".sdf"):
            from molmod.io import SDFReader
            return SDFReader(filename).next()
        elif filename.endswith(".xyz"):
            from molmod.io import XYZReader
            xyz_reader = XYZReader(filename)
            title, coordinates = xyz_reader.next()
            return Molecule(xyz_reader.numbers, coordinates, title)
        else:
            raise ValueError("Could not determine file format for %s." % filename)

    size = property(lambda self: self.numbers.shape[0])

    @cached
    def distance_matrix(self):
        """The matrix with all atom pair distances"""
        from molmod.ext import molecules_distance_matrix
        return molecules_distance_matrix(self.coordinates)

    @cached
    def mass(self):
        """The total mass of the molecule"""
        return self.masses.sum()

    @cached
    def com(self):
        """The center of mass of the molecule"""
        return (self.coordinates*self.masses.reshape((-1,1))).sum(axis=0)/self.mass

    @cached
    def inertia_tensor(self):
        """The intertia tensor of the molecule"""
        result = numpy.zeros((3,3), float)
        for i in xrange(self.size):
            r = self.coordinates[i] - self.com
            # the diagonal term
            result.ravel()[::4] += self.masses[i]*(r**2).sum()
            # the outer product term
            result -= self.masses[i]*numpy.outer(r,r)
        return result

    @cached
    def chemical_formula(self):
        """The chemical formula of the molecule"""
        counts = {}
        for number in self.numbers:
            counts[number] = counts.get(number, 0)+1
        items = []
        for number, count in sorted(counts.iteritems(), reverse=True):
            if count == 1:
                items.append(periodic[number].symbol)
            else:
                items.append("%s%i" % (periodic[number].symbol, count))
        return "".join(items)

    def set_default_masses(self):
        """Set self.masses based on self.numbers"""
        self.masses = numpy.array([periodic[n].mass for n in self.numbers])

    def set_default_graph(self):
        """Set self.graph to the default graph, see MolecularGraph.from_geometry"""
        self.graph = MolecularGraph.from_geometry(self)

    def write_to_file(self, filename):
        """Write the molecule geometry to a file

           The file format is inferred from the extensions. Currently supported
           formats are: *.xyz, *.cml

           Argument:
             filename  --  a filename
        """
        # TODO: give all file format writers the same API
        if filename.endswith('.cml'):
            from molmod.io import dump_cml
            dump_cml(filename, [self])
        elif filename.endswith('.xyz'):
            from molmod.io import XYZWriter
            xyz_writer = XYZWriter(filename, [periodic[n].symbol for n in self.numbers])
            if hasattr(self, "title"):
                title = self.title
            else:
                title = "Sorry, no titles today..."
            xyz_writer.dump(title, self.coordinates)
            del xyz_writer
        else:
            raise ValueError("Could not determine file format for %s." % filename)


