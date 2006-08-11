# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --


from molmod.graphs import Graph, SymmetricGraph, MatchFilterParameterized, Criterium
from molmod.binning import IntraAnalyseNeighboringObjects, PositionedObject, SparseBinnedObjects
from molmod.data import bonds

import math, numpy

__all__ = [
    "MolecularCriterium",
    "BinaryMolecularCriterium", "MolecularOr", "MolecularAnd",
    "AtomNumberCriterium", "AtomNumberRequire", "AtomNumberRefuse",
    "BondTypeCriterium", "BondTypeRequire", "BondTypeRefuse",
    "NumNeighboursCriterium", "NumNeighboursRequire", "NumNeighboursRefuse",
    "CriteriaSets", "CriteriaSet", "BondSets", "BendSets", "DihedralSets",
    "OutOfPlaneAngleSets", "OutOfPlaneDistanceSets", "LongRangePairSets",
    "MolecularGraph",
]

# Elementary criteria for MatchFilters

class MolecularCriterium(Criterium):
    def __init__(self, *parameters):
        self.molecular_graph = None
        Criterium.__init__(self, *parameters)

    def set_molecular_graph(self, molecular_graph):
        self.molecular_graph = molecular_graph


class BinaryMolecularCriterium(MolecularCriterium):
    def __init__(self, criterium1, criterium2):
        self.criterium1 = criterium1
        self.criterium2 = criterium2
        MolecularCriterium.__init__(self, None)

    def get_tag(self):
        return (self.__class__, self.criterium1.get_tag(), self.criterium2.get_tag())

    def set_molecular_graph(self, molecular_graph):
        self.molecular_graph = molecular_graph
        self.criterium1.set_molecular_graph(molecular_graph)
        self.criterium2.set_molecular_graph(molecular_graph)


class MolecularOr(BinaryMolecularCriterium):
    def __call__(self, argument):
        return self.criterium1(argument) or self.criterium2(argument)


class MolecularAnd(BinaryMolecularCriterium):
    def __call__(self, argument):
        return self.criterium1(argument) and self.criterium2(argument)


class AtomNumberCriterium(MolecularCriterium):
    def __init__(self, number):
        self.number = number
        MolecularCriterium.__init__(self, number)


class AtomNumberRequire(AtomNumberCriterium):
    def __call__(self, index):
        return self.molecular_graph.molecule.numbers[index] == self.number


class AtomNumberRefuse(AtomNumberCriterium):
    def __call__(self, index):
        return self.molecular_graph.molecule.numbers[index] != self.number


class BondTypeCriterium(MolecularCriterium):
    def __init__(self, bond_type):
        self.bond_type = bond_type
        MolecularCriterium.__init__(self, bond_type)


class BondTypeRequire(BondTypeCriterium):
    def __call__(self, pair):
        return self.molecular_graph.bonds[pair] == self.bond_type


class BondTypeRefuse(BondTypeCriterium):
    def __call__(self, pair):
        return self.molecular_graph.bonds[pair] != self.bond_type


class NumNeighboursCriterium(MolecularCriterium):
    def __init__(self, num_neighbors):
        self.num_neighbors = num_neighbors
        MolecularCriterium.__init__(self, num_neighbors)


class NumNeighboursRequire(NumNeighboursCriterium):
    def __call__(self, index):
        return len(self.molecular_graph.neighbors[index]) == self.num_neighbors


class NumNeighboursRefuse(NumNeighboursCriterium):
    def __call__(self, index):
        return len(self.molecular_graph.neighbors[index]) != self.num_neighbors


# Predefined sets of criteria: bonds, angles, dihedrals

class CriteriaSets(object):
    def __init__(self, subpairs, initiator, calculation_tags, sets):
        self.subpairs = subpairs   # pairs of bonded atoms
        self.initiator = initiator # the central ellement, the one that is transformed onto itself by most of the symmetries
        self.calculation_tags = calculation_tags # tags that indicate which atoms are similar due to symmetrie in the topology, not atom labels taken into account
        self.sets = sets

    def yield_criteria(self):
        for item in self.sets:
            yield item.tag, item.atom_criteria, item.bond_criteria, item.filter_tags

    def yield_tags(self):
        for item in self.sets:
            yield item.tag


class CriteriaSet(object):
    def __init__(self, tag, short=None, extensive=None, filter_tags=True):
        self.tag = tag
        self.atom_criteria = {}
        self.bond_criteria = {}
        self.filter_tags = filter_tags
        if short != None:
            atoms, bonds = short
            if atoms != None:
                for index, number in enumerate(atoms):
                    if number > 0:
                        self.atom_criteria[index] = AtomNumberRequire(number)
                    elif number < 0:
                        self.atom_criteria[index] = AtomNumberRefuse(number)
            if bonds != None:
                for index, bond_type in bonds.iteritems:
                    if bond_type > 0:
                        self.bond_criteria[frozenset(index)] = BondTypeRequire(bond_type)
                    elif bond_type < 0:
                        self.bond_criteria[frozenset(index)] = BondTypeRefuse(bond_type)
        if extensive != None:
            atoms, bonds = extensive
            if atoms != None:
                self.atom_criteria.update(atoms)
            if bonds != None:
                self.bond_criteria.update(bonds)


class BondSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(0, 1)], 0, {0: 0, 1: 0}, sets)


class BendSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(0, 1), (1, 2)], 1, {0: 1, 1: 0, 2: 1}, sets)


class DihedralSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(0, 1), (1, 2), (2, 3)], 1, {0: 0, 1: 1, 2: 1, 3: 0}, sets)


class OutOfPlaneAngleSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(1, 0), (1, 2), (1, 3)], 1, {0: 0, 1: 1, 2: 2, 3: 2}, sets)


class OutOfPlaneDistanceSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(0, 1), (0, 2), (0, 3)], 0, {0: 0, 1: 1, 2: 1, 3: 1}, sets)


class LongRangePairSets(CriteriaSets):
    def __init__(self, sets):
        CriteriaSets.__init__(self, [(0, 1)], 0, {0: 0, 1: 0}, sets)

    def bond_excludes(self):
        return BondSets(self.sets)

    def bend_excludes(self):
        new_sets = [
            CriteriaSet(
                set.tag,
                extensive=({0: set.atom_criteria[0], 2: set.atom_criteria[1]}, None),
                filter_tags=set.filter_tags
            )
            for set
            in self.sets
        ]
        return BendSets(new_sets)

    def dihedral_excludes(self):
        new_sets = [
            CriteriaSet(
                set.tag,
                extensive=({0: set.atom_criteria[0], 3: set.atom_criteria[1]}, None),
                filter_tags=set.filter_tags
            )
            for set
            in self.sets
        ]
        return DihedralSets(new_sets)


#

class MolecularGraph(Graph):
    def __init__(self, molecule):
        self.molecule = molecule

        def yield_positioned_atoms():
            for index in xrange(len(self.molecule.numbers)):
                yield PositionedObject(index, self.molecule.coordinates[index])

        binned_atoms = SparseBinnedObjects(yield_positioned_atoms(), bonds.max_length*bonds.bond_tolerance)

        def compare_function(positioned1, positioned2):
            delta = positioned2.coordinate - positioned1.coordinate
            distance = math.sqrt(numpy.dot(delta, delta))
            if distance < binned_atoms.gridsize:
                bond_order = bonds.bonded(self.molecule.numbers[positioned1.id], self.molecule.numbers[positioned2.id], distance)
                if bond_order != None:
                    return bond_order, distance

        bond_data = list(
            (frozenset([positioned.id for positioned in key]), data)
            for key, data
            in IntraAnalyseNeighboringObjects(binned_atoms, compare_function)()
        )
        pairs = set(key for key, data in bond_data)
        self.bond_orders = dict([(key, data[0]) for key, data in bond_data])
        self.bond_lengths = dict([(key, data[1]) for key, data in bond_data])
        Graph.__init__(self, pairs)

    def yield_subgraphs(self, criteria_sets):
        subgraph = SymmetricGraph(criteria_sets.subpairs, criteria_sets.initiator)
        for tag, atom_criteria, bond_criteria, filter_tags in criteria_sets.yield_criteria():
            for atom_criterium in atom_criteria.itervalues():
                atom_criterium.set_molecular_graph(self)
            for bond_criterium in bond_criteria.itervalues():
                bond_criterium.set_molecular_graph(self)
            graph_filter = MatchFilterParameterized(
                subgraph,
                criteria_sets.calculation_tags,
                atom_criteria,
                bond_criteria,
                filter_tags
            )

            for match in subgraph.yield_matching_subgraphs(self):
                for parsed in graph_filter.parse(match):
                    yield tag, parsed
