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


import numpy
from molmod.data.periodic import periodic
from molmod.units import unified
from molmod.graphs import CriteriaSet, MatchGenerator
from molmod.molecular_graphs import MolecularGraph, generate_molecular_graph, \
    BondMatchDefinition, BendingAngleMatchDefinition, \
    DihedralAngleMatchDefinition, MolecularExactMatchDefinition, HasAtomNumber
from molmod.graphs import Graph


__all__ = ["Error", "PSFFile"]

class Error(Exception):
    pass


class PSFFile(object):
    "A very simplistic and limited implementation of the PSF file format."

    def __init__(self, filename=None):
        if filename is None:
            self.clear()
        else:
            self.read_from_file(filename)

    def clear(self):
        self.title = None
        self.numbers = numpy.zeros(0,int)
        self.atom_types = [] # the atom_types in the second column, used to associate ff parameters
        self.charges = [] # ff charges
        self.names = [] # a name that is unique for the molecule composition and connectivity
        self.molecules = numpy.zeros(0,int) # a counter for each molecule
        self.bonds = numpy.zeros((0,2),int)
        self.bends = numpy.zeros((0,3),int)
        self.dihedrals = numpy.zeros((0,4),int)

        self.name_cache = {}

    def read_from_file(self, filename):
        self.clear()
        f = file(filename)
        # A) check the first line
        line = f.next()
        if not line.startswith("PSF"):
            raise Error("Error while reading: A PSF file must start with a line 'PSF'.")
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
        f.close()
        # C) interpret the supported sections
        # C.1) The title
        self.title = sections['TITLE'][0]
        self.molecules = []
        self.numbers = []
        # C.2) The atoms and molecules
        for line in sections['ATOM']:
            words = line.split()
            self.atom_types.append(words[5])
            self.charges.append(float(words[6]))
            self.names.append(words[3])
            self.molecules.append(int(words[2]))
            atom = periodic[words[4]]
            if atom is None:
                self.numbers.append(0)
            else:
                self.numbers.append(periodic[words[4]].number)
        self.molecules = numpy.array(self.molecules)-1
        self.numbers = numpy.array(self.numbers)
        self.charges = numpy.array(self.charges)
        # C.3) The bonds section
        tmp = []
        for line in sections['BOND']:
            tmp.extend(int(word) for word in line.split())
        self.bonds = numpy.reshape(numpy.array(tmp), (-1,2))-1
        # C.4) The bends section
        tmp = []
        for line in sections['THETA']:
            tmp.extend(int(word) for word in line.split())
        self.bends = numpy.reshape(numpy.array(tmp), (-1,3))-1
        # C.5) The dihedral section
        tmp = []
        for line in sections['PHI']:
            tmp.extend(int(word) for word in line.split())
        self.dihedrals = numpy.reshape(numpy.array(tmp), (-1,4))-1

    def _get_name(self, graph, group=None):
        def compare(g1, g2):
            if len(g1.pairs) != len(g2.pairs): return False
            if len(g1.nodes) != len(g2.nodes): return False
            atom_criteria = dict((index, HasAtomNumber(number)) for index, number in zip(g1.nodes, g1.numbers))
            md = MolecularExactMatchDefinition(g1, [CriteriaSet(atom_criteria)])
            try:
                match = MatchGenerator(md,debug=False)(g2,one_match=True).next()
                return True
            except StopIteration:
                return False

        #print "numbers", graph.numbers
        if group is not None:
            #print "selecting", group
            graph = graph.subgraph(group)
        #print "numbers", graph.numbers
        for name, ref_graph in self.name_cache.iteritems():
            #print "trying", name, ref_graph, "==", graph
            if compare(graph, ref_graph):
                #print "YES", name
                return name
            #print "NO"
        name = "NM%02i" % len(self.name_cache)
        #print "setting", name, graph
        self.name_cache[name] = graph
        #print "-"*20
        return name

    def write_to_file(self, filename):
        f = file(filename, 'w')
        self.dump(f)
        f.close()

    def dump(self, f):
        # header
        print >> f, "PSF"
        print >> f

        # title
        print >> f, "      1 !NTITLE"
        print >> f, self.title
        print >> f

        # atoms
        print >> f, "% 7i !NATOM" % len(self.numbers)
        if len(self.numbers) > 0:
            for index, (number, atom_type, charge, name, molecule) in enumerate(zip(self.numbers,self.atom_types,self.charges,self.names,self.molecules)):
                atom = periodic[number]
                print >> f, "% 7i % 4s % 4i NAME % 6s % 6s % 8.4f % 12.6f 0" % (
                    index + 1,
                    name,
                    molecule + 1,
                    atom.symbol,
                    atom_type,
                    charge,
                    atom.mass/unified,
                )
        print >> f

        # bonds
        print >> f, "% 7i !NBOND" % len(self.bonds)
        if len(self.bonds) > 0:
            tmp = []
            for bond in self.bonds:
                tmp.extend(bond+1)
                if len(tmp) >= 8:
                    print >> f, " ".join("% 7i" % v for v in tmp[:8])
                    tmp = tmp[8:]
            if len(tmp) > 0:
                print >> f, " ".join("% 7i" % v for v in tmp)
        print >> f

        # bends
        print >> f, "% 7i !NTHETA" % len(self.bends)
        if len(self.bends) > 0:
            tmp = []
            for bend in self.bends:
                tmp.extend(bend+1)
                if len(tmp) >= 9:
                    print >> f, " " + (" ".join("% 6i" % v for v in tmp[:9]))
                    tmp = tmp[9:]
            if len(tmp) > 0:
                print >> f, " " + (" ".join("% 6i" % v for v in tmp))
        print >> f

        # dihedrals
        print >> f, "% 7i !NPHI" % len(self.dihedrals)
        if len(self.dihedrals) > 0:
            tmp = []
            for dihedral in self.dihedrals:
                tmp.extend(dihedral+1)
                if len(tmp) >= 8:
                    print >> f, " " + (" ".join("% 6i" % v for v in tmp[:8]))
                    tmp = tmp[8:]
            if len(tmp) > 0:
                print >> f, " " + (" ".join("% 6i" % v for v in tmp))
        print >> f

        # not implemented fields
        print >> f, "      0 !NIMPHI"
        print >> f
        print >> f, "      0 !NDON"
        print >> f
        print >> f, "      0 !NNB"
        print >> f
        print >> f, "      0 !NGRP"
        print >> f

    def add_molecule(self, molecule, atom_types=None, charges=None, split=True):
        molecular_graph = generate_molecular_graph(molecule)
        self.add_molecular_graph(molecular_graph, atom_types, charges, split)

    def add_molecular_graph(self, molecular_graph, atom_types=None, charges=None, split=True):
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
            groups = molecular_graph.get_indexes_per_independent_graph()
            names = [self._get_name(molecular_graph, group) for group in groups]
            group_indices = numpy.zeros(new, int)
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


        # add bonds
        match_generator = MatchGenerator(BondMatchDefinition([CriteriaSet()]))
        tmp = [(
            match.get_destination(0),
            match.get_destination(1),
        ) for match in match_generator(molecular_graph)]
        tmp.sort()
        new = len(tmp)
        if new > 0:
            prev = len(self.bonds)
            self.bonds.resize((prev + len(tmp), 2))
            self.bonds[-len(tmp):] = tmp
            self.bonds[-len(tmp):] += offset

        # add bends
        match_generator = MatchGenerator(BendingAngleMatchDefinition([CriteriaSet()]))
        tmp = [(
            match.get_destination(0),
            match.get_destination(1),
            match.get_destination(2),
        ) for match in match_generator(molecular_graph)]
        tmp.sort()
        new = len(tmp)
        if new > 0:
            prev = len(self.bends)
            self.bends.resize((prev + len(tmp), 3))
            self.bends[-len(tmp):] = tmp
            self.bends[-len(tmp):] += offset

        # add dihedrals
        match_generator = MatchGenerator(DihedralAngleMatchDefinition([CriteriaSet()]))
        tmp = [(
            match.get_destination(0),
            match.get_destination(1),
            match.get_destination(2),
            match.get_destination(3),
        ) for match in match_generator(molecular_graph)]
        tmp.sort()
        new = len(tmp)
        if new > 0:
            prev = len(self.dihedrals)
            self.dihedrals.resize((prev + len(tmp), 4))
            self.dihedrals[-len(tmp):] = tmp
            self.dihedrals[-len(tmp):] += offset

    def get_graph(self):
        return Graph(pairs=set(frozenset(bond) for bond in self.bonds))

    def get_molecular_graph(self, labels=None):
        return MolecularGraph(set(frozenset(bond) for bond in self.bonds), self.numbers, labels)

    def get_groups(self):
        groups = []
        for a_index, m_index in enumerate(self.molecules):
            if m_index >= len(groups):
                groups.append([a_index])
            else:
                groups[m_index].append(a_index)
        return groups



