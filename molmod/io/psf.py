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
"""Tools for reading and writing PSF (Protein Structure File) files

   This format is orignally developed in conjunction with the CHARMM program,
   but is also used by CP2K as a generic format to define the molecular bond
   graph and other topological aspects of a molecular model. This module just
   creates files that can be read with CP2K.

   Due to the lack of public definition of a PSF file format, several variations
   exist. Therefore we do not make any attempt to follow the original format
   strictly, nor one of these variations. It would be more interesting to
   define a new public standard for molecular topologies, which is also suitable
   for non-protein systems.
"""


from __future__ import print_function, division

import numpy as np

from molmod.periodic import periodic
from molmod.units import unified
from molmod.graphs import CriteriaSet, GraphSearch
from molmod.molecular_graphs import MolecularGraph, BondPattern, \
    BendingAnglePattern, DihedralAnglePattern, OutOfPlanePattern, \
    HasNumNeighbors
from molmod.graphs import Graph
from molmod.io.common import FileFormatError


__all__ = ["PSFFile"]


class PSFFile(object):
    "A very simplistic and limited implementation of the PSF file format"

    def __init__(self, filename=None):
        """
           Argument:
            | ``filename``  --  When not given, an empty data structure is
                                created, otherwise the file is loaded from disk
        """
        if filename is None:
            self.clear()
        else:
            self.read_from_file(filename)

    def clear(self):
        """Clear the contents of the data structure"""
        self.title = None
        self.numbers = np.zeros(0, int)
        self.atom_types = [] # the atom_types in the second column, used to associate ff parameters
        self.charges = [] # ff charges
        self.names = [] # a name that is unique for the molecule composition and connectivity
        self.molecules = np.zeros(0, int) # a counter for each molecule
        self.bonds = np.zeros((0, 2), int)
        self.bends = np.zeros((0, 3), int)
        self.dihedrals = np.zeros((0, 4), int)
        self.impropers = np.zeros((0, 4), int)

        self.name_cache = {}

    def read_from_file(self, filename):
        """Load a PSF file"""
        self.clear()
        with open(filename) as f:
            # A) check the first line
            line = next(f)
            if not line.startswith("PSF"):
                raise FileFormatError("Error while reading: A PSF file must start with a line 'PSF'.")
            # B) read in all the sections, without interpreting them
            current_section = None
            sections = {}
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                elif "!N" in line:
                    words = line.split()
                    current_section = []
                    section_name = words[1][2:]
                    if section_name.endswith(":"):
                        section_name = section_name[:-1]
                    sections[section_name] = current_section
                else:
                    current_section.append(line)
        # C) interpret the supported sections
        # C.1) The title
        self.title = sections['TITLE'][0]
        molecules = []
        numbers = []
        # C.2) The atoms and molecules
        for line in sections['ATOM']:
            words = line.split()
            self.atom_types.append(words[5])
            self.charges.append(float(words[6]))
            self.names.append(words[3])
            molecules.append(int(words[2]))
            atom = periodic[words[4]]
            if atom is None:
                numbers.append(0)
            else:
                numbers.append(periodic[words[4]].number)
        self.molecules = np.array(molecules)-1
        self.numbers = np.array(numbers)
        self.charges = np.array(self.charges)
        # C.3) The bonds section
        tmp = []
        for line in sections['BOND']:
            tmp.extend(int(word) for word in line.split())
        self.bonds = np.reshape(np.array(tmp), (-1, 2))-1
        # C.4) The bends section
        tmp = []
        for line in sections['THETA']:
            tmp.extend(int(word) for word in line.split())
        self.bends = np.reshape(np.array(tmp), (-1, 3))-1
        # C.5) The dihedral section
        tmp = []
        for line in sections['PHI']:
            tmp.extend(int(word) for word in line.split())
        self.dihedrals = np.reshape(np.array(tmp), (-1, 4))-1
        # C.6) The improper section
        tmp = []
        for line in sections['IMPHI']:
            tmp.extend(int(word) for word in line.split())
        self.impropers = np.reshape(np.array(tmp), (-1, 4))-1

    def _get_name(self, graph, group=None):
        """Convert a molecular graph into a unique name

           This method is not sensitive to the order of the atoms in the graph.
        """
        if group is not None:
            graph = graph.get_subgraph(group, normalize=True)

        fingerprint = graph.fingerprint.tobytes()
        name = self.name_cache.get(fingerprint)
        if name is None:
            name = "NM%02i" % len(self.name_cache)
            self.name_cache[fingerprint] = name
        return name

    def write_to_file(self, filename):
        """Write the data structure to a file"""
        with open(filename, 'w') as f:
            self.dump(f)

    def dump(self, f):
        """Dump the data structure to a file-like object"""
        # header
        print("PSF", file=f)
        print(file=f)

        # title
        print("      1 !NTITLE", file=f)
        print(self.title, file=f)
        print(file=f)

        # atoms
        print("% 7i !NATOM" % len(self.numbers), file=f)
        if len(self.numbers) > 0:
            for index, (number, atom_type, charge, name, molecule) in enumerate(zip(self.numbers, self.atom_types, self.charges, self.names, self.molecules)):
                atom = periodic[number]
                print("% 7i % 4s % 4i NAME % 6s % 6s % 8.4f % 12.6f 0" % (
                    index + 1,
                    name,
                    molecule + 1,
                    atom.symbol,
                    atom_type,
                    charge,
                    atom.mass/unified,
                ), file=f)
        print(file=f)

        # bonds
        print("% 7i !NBOND" % len(self.bonds), file=f)
        if len(self.bonds) > 0:
            tmp = []
            for bond in self.bonds:
                tmp.extend(bond+1)
                if len(tmp) >= 8:
                    print(" ".join("% 7i" % v for v in tmp[:8]), file=f)
                    tmp = tmp[8:]
            if len(tmp) > 0:
                print(" ".join("% 7i" % v for v in tmp), file=f)
        print(file=f)

        # bends
        print("% 7i !NTHETA" % len(self.bends), file=f)
        if len(self.bends) > 0:
            tmp = []
            for bend in self.bends:
                tmp.extend(bend+1)
                if len(tmp) >= 9:
                    print(" " + (" ".join("% 6i" % v for v in tmp[:9])), file=f)
                    tmp = tmp[9:]
            if len(tmp) > 0:
                print(" " + (" ".join("% 6i" % v for v in tmp)), file=f)
        print(file=f)

        # dihedrals
        print("% 7i !NPHI" % len(self.dihedrals), file=f)
        if len(self.dihedrals) > 0:
            tmp = []
            for dihedral in self.dihedrals:
                tmp.extend(dihedral+1)
                if len(tmp) >= 8:
                    print(" " + (" ".join("% 6i" % v for v in tmp[:8])), file=f)
                    tmp = tmp[8:]
            if len(tmp) > 0:
                print(" " + (" ".join("% 6i" % v for v in tmp)), file=f)
        print(file=f)

        # impropers
        print("% 7i !NIMPHI" % len(self.impropers), file=f)
        if len(self.impropers) > 0:
            tmp = []
            for improper in self.impropers:
                tmp.extend(improper+1)
                if len(tmp) >= 8:
                    print(" " + (" ".join("% 6i" % v for v in tmp[:8])), file=f)
                    tmp = tmp[8:]
            if len(tmp) > 0:
                print(" " + (" ".join("% 6i" % v for v in tmp)), file=f)
        print(file=f)

        # not implemented fields
        print("      0 !NDON", file=f)
        print(file=f)
        print("      0 !NNB", file=f)
        print(file=f)
        print("      0 !NGRP", file=f)
        print(file=f)

    def add_molecule(self, molecule, atom_types=None, charges=None, split=True):
        """Add the graph of the molecule to the data structure

           The molecular graph is estimated from the molecular geometry based on
           interatomic distances.

           Argument:
            | ``molecule``  --  a Molecule instance

           Optional arguments:
            | ``atom_types``  --  a list with atom type strings
            | ``charges``  --  The net atom charges
            | ``split``  --  When True, the molecule is split into disconnected
                             molecules [default=True]
        """
        molecular_graph = MolecularGraph.from_geometry(molecule)
        self.add_molecular_graph(molecular_graph, atom_types, charges, split, molecule)

    def add_molecular_graph(self, molecular_graph, atom_types=None, charges=None, split=True, molecule=None):
        """Add the molecular graph to the data structure

           Argument:
            | ``molecular_graph``  --  a MolecularGraph instance

           Optional arguments:
            | ``atom_types``  --  a list with atom type strings
            | ``charges``  --  The net atom charges
            | ``split``  --  When True, the molecule is split into disconnected
                             molecules [default=True]
        """
        # add atom numbers and molecule indices
        new = len(molecular_graph.numbers)
        if new == 0: return
        prev = len(self.numbers)
        offset = prev
        self.numbers.resize(prev + new)
        self.numbers[-new:] = molecular_graph.numbers
        if atom_types is None:
            atom_types = [periodic[number].symbol for number in molecular_graph.numbers]
        self.atom_types.extend(atom_types)
        if charges is None:
            charges = [0.0]*len(molecular_graph.numbers)
        self.charges.extend(charges)
        self.molecules.resize(prev + new)

        # add names (autogenerated)
        if split:
            groups = molecular_graph.independent_vertices
            names = [self._get_name(molecular_graph, group) for group in groups]
            group_indices = np.zeros(new, int)
            for group_index, group in enumerate(groups):
                for index in group:
                    group_indices[index] = group_index
            self.names.extend([names[group_index] for group_index in group_indices])
            if prev == 0:
                self.molecules[:] = group_indices
            else:
                self.molecules[-new:] = self.molecules[-new]+group_indices+1
        else:
            if prev == 0:
                self.molecules[-new:] = 0
            else:
                self.molecules[-new:] = self.molecules[-new]+1
            name = self._get_name(molecular_graph)
            self.names.extend([name]*new)

        self._add_graph_bonds(molecular_graph, offset, atom_types, molecule)
        self._add_graph_bends(molecular_graph, offset, atom_types, molecule)
        self._add_graph_dihedrals(molecular_graph, offset, atom_types, molecule)
        self._add_graph_impropers(molecular_graph, offset, atom_types, molecule)

    def _add_graph_bonds(self, molecular_graph, offset, atom_types, molecule):
        # add bonds
        match_generator = GraphSearch(BondPattern([CriteriaSet()]))
        tmp = sorted([(
            match.get_destination(0),
            match.get_destination(1),
        ) for match in match_generator(molecular_graph)])
        new = len(tmp)
        if new > 0:
            prev = len(self.bonds)
            self.bonds.resize((prev + len(tmp), 2))
            self.bonds[-len(tmp):] = tmp
            self.bonds[-len(tmp):] += offset

    def _add_graph_bends(self, molecular_graph, offset, atom_types, molecule):
        # add bends
        match_generator = GraphSearch(BendingAnglePattern([CriteriaSet()]))
        tmp = sorted([(
            match.get_destination(0),
            match.get_destination(1),
            match.get_destination(2),
        ) for match in match_generator(molecular_graph)])
        new = len(tmp)
        if new > 0:
            prev = len(self.bends)
            self.bends.resize((prev + len(tmp), 3))
            self.bends[-len(tmp):] = tmp
            self.bends[-len(tmp):] += offset

    def _add_graph_dihedrals(self, molecular_graph, offset, atom_types, molecule):
        # add dihedrals
        match_generator = GraphSearch(DihedralAnglePattern([CriteriaSet()]))
        tmp = sorted([(
            match.get_destination(0),
            match.get_destination(1),
            match.get_destination(2),
            match.get_destination(3),
        ) for match in match_generator(molecular_graph)])
        new = len(tmp)
        if new > 0:
            prev = len(self.dihedrals)
            self.dihedrals.resize((prev + len(tmp), 4))
            self.dihedrals[-len(tmp):] = tmp
            self.dihedrals[-len(tmp):] += offset

    def _add_graph_impropers(self, molecular_graph, offset, atom_types, molecule):
        # add improper dihedrals, only when center has three bonds
        match_generator = GraphSearch(OutOfPlanePattern([CriteriaSet(
            vertex_criteria={0: HasNumNeighbors(3)},
        )], vertex_tags={1:1}))
        tmp = sorted([(
            match.get_destination(0),
            match.get_destination(1),
            match.get_destination(2),
            match.get_destination(3),
        ) for match in match_generator(molecular_graph)])
        new = len(tmp)
        if new > 0:
            prev = len(self.impropers)
            self.impropers.resize((prev + len(tmp), 4))
            self.impropers[-len(tmp):] = tmp
            self.impropers[-len(tmp):] += offset

    def get_graph(self):
        """Return the bond graph represented by the data structure"""
        return Graph(self.bonds)

    def get_molecular_graph(self):
        """Return the molecular graph represented by the data structure"""
        return MolecularGraph(self.bonds, self.numbers)

    def get_groups(self):
        """Return a list of groups of atom indexes

           Each atom in a group belongs to the same molecule or residue.
        """
        groups = []
        for a_index, m_index in enumerate(self.molecules):
            if m_index >= len(groups):
                groups.append([a_index])
            else:
                groups[m_index].append(a_index)
        return groups
