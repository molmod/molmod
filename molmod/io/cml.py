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
"""Basic support for the Chemical Markup Language

   Not all features of the CML standard are supported in this module, only the
   basic aspects that are relevant for computational chemistry. For more info
   on the CML format, visit: http://cml.sourceforge.net/

   In this module, only atoms, their 3D coordinates, atom numbers, bonds and
   bond orders are supported.
"""


from __future__ import division

from builtins import range
import numpy as np

from xml.sax import make_parser
from xml.sax.handler import feature_namespaces, ContentHandler

from molmod.units import angstrom
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.periodic import periodic


__all__ = ["load_cml", "dump_cml"]


class CMLMoleculeLoader(ContentHandler):
    """A ContentHandler that reads the essentials out of a CML file"""
    def __init__(self):
        self.molecules = []
        self.current_title = None # current molecule
        ContentHandler.__init__(self)

    atom_exclude = frozenset(['id', 'elementType', 'x3', 'y3', 'z3', 'x2', 'y2'])
    bond_exclude = frozenset(['id', 'atomRefs2', 'order'])
    molecule_exclude = frozenset(['id', 'xmlns'])

    def _get_extra(self, attrs, exclude):
        """Read the extra properties, taking into account an exclude list"""
        result = {}
        for key in attrs.getNames():
            if key not in exclude:
                result[str(key)] = str(attrs[key])
        return result

    def startElement(self, name, attrs):
        #print "START", name
        # If it's not a comic element, ignore it
        if name == 'molecule':
            self.current_title = str(attrs.get('id', 'No Title'))
            self.current_numbers = []
            self.current_coordinates = []
            self.current_atom_names = []
            self.current_bonds = []
            self.current_extra = self._get_extra(attrs, self.molecule_exclude)
            self.current_atoms_extra = {}
        elif self.current_title is not None:
            if name == 'atom':
                atom_name = attrs.get('id', None)
                if atom_name is None:
                    return
                symbol = attrs.get('elementType', None)
                if symbol is None:
                    return
                try:
                    x = float(attrs.get('x3', None))
                    y = float(attrs.get('y3', None))
                    z = float(attrs.get('z3', None))
                except ValueError:
                    return
                atom_record = periodic[str(symbol)]
                if atom_record is None:
                    return
                self.current_atom_names.append(atom_name)
                self.current_numbers.append(atom_record.number)
                self.current_coordinates.append([x, y, z])
                # find potential extra attributes
                extra = self._get_extra(attrs, self.atom_exclude)
                if len(extra) > 0:
                    self.current_atoms_extra[len(self.current_numbers)-1] = extra
            elif name == 'bond':
                refs = str(attrs.get('atomRefs2', ''))
                if refs.count(' ') != 1:
                    return
                name1, name2 = refs.split(' ')
                extra = self._get_extra(attrs, self.bond_exclude)
                self.current_bonds.append((name1, name2, extra))

    def endElement(self, name):
        #print "END", name
        if name == 'molecule':
            if len(self.current_numbers) > 0:
                self.current_coordinates = np.array(self.current_coordinates)*angstrom
                molecule = Molecule(self.current_numbers, self.current_coordinates, self.current_title)
                molecule.extra = self.current_extra
                molecule.atoms_extra = self.current_atoms_extra

                name_to_index = {}
                for counter, name in enumerate(self.current_atom_names):
                    name_to_index[name] = counter

                edges = set()
                current_bonds_extra = {}
                for name1, name2, extra in self.current_bonds:
                    i1 = name_to_index.get(name1)
                    i2 = name_to_index.get(name2)
                    if i1 is not None and i2 is not None:
                        edge = frozenset([i1, i2])
                        if len(extra) > 0:
                            current_bonds_extra[edge] = extra
                        edges.add(edge)

                molecule.bonds_extra = current_bonds_extra
                if len(edges) == 0:
                    molecule.graph = None
                else:
                    molecule.graph = MolecularGraph(edges, self.current_numbers)
                del self.current_atom_names
                del self.current_bonds

                self.molecules.append(molecule)
            self.current_title = None


def load_cml(cml_filename):
    """Load the molecules from a CML file

       Argument:
        | ``cml_filename``  --  The filename of a CML file.

       Returns a list of molecule objects with optional molecular graph
       attribute and extra attributes.
    """
    parser = make_parser()
    parser.setFeature(feature_namespaces, 0)
    dh = CMLMoleculeLoader()
    parser.setContentHandler(dh)
    parser.parse(cml_filename)
    return dh.molecules


def _dump_cml_molecule(f, molecule):
    """Dump a single molecule to a CML file

       Arguments:
        | ``f``  --  a file-like object
        | ``molecule``  --  a Molecule instance
    """
    extra = getattr(molecule, "extra", {})
    attr_str = " ".join("%s='%s'" % (key, value) for key, value in extra.items())
    f.write(" <molecule id='%s' %s>\n" % (molecule.title, attr_str))
    f.write("  <atomArray>\n")
    atoms_extra = getattr(molecule, "atoms_extra", {})
    for counter, number, coordinate in zip(range(molecule.size), molecule.numbers, molecule.coordinates/angstrom):
        atom_extra = atoms_extra.get(counter, {})
        attr_str = " ".join("%s='%s'" % (key, value) for key, value in atom_extra.items())
        f.write("   <atom id='a%i' elementType='%s' x3='%s' y3='%s' z3='%s' %s />\n" % (
            counter, periodic[number].symbol, coordinate[0],  coordinate[1],
            coordinate[2], attr_str,
        ))
    f.write("  </atomArray>\n")
    if molecule.graph is not None:
        bonds_extra = getattr(molecule, "bonds_extra", {})
        f.write("  <bondArray>\n")
        for edge in molecule.graph.edges:
            bond_extra = bonds_extra.get(edge, {})
            attr_str = " ".join("%s='%s'" % (key, value) for key, value in bond_extra.items())
            i1, i2 = edge
            f.write("   <bond atomRefs2='a%i a%i' %s />\n" % (i1, i2, attr_str))
        f.write("  </bondArray>\n")
    f.write(" </molecule>\n")


def dump_cml(f, molecules):
    """Write a list of molecules to a CML file

       Arguments:
        | ``f``  --  a filename of a CML file or a file-like object
        | ``molecules``  --  a list of molecule objects.
    """
    if isinstance(f, str):
        f = open(f, "w")
        close = True
    else:
        close = False
    f.write("<?xml version='1.0'?>\n")
    f.write("<list xmlns='http://www.xml-cml.org/schema'>\n")
    for molecule in molecules:
        _dump_cml_molecule(f, molecule)
    f.write("</list>\n")
    if close:
        f.close()
