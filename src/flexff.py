# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

from molmod.molecular_graphs import BondMatchDefinition, BendingAngleMatchDefinition, DihedralAngleMatchDefinition
from molmod.graphs import MatchGenerator, CriteriaSet


__all__ = [
    "EnergyTerm",
    "ValenceTerm",
    "BondStretchTerm",
    "ThreeBodyTerm",
    "BendingCosineTerm",
    "UreyBradleyTerm",
    "FourBodyTerm",
    "DihedralCosineTerm",
    "OneFourTerm",
    "NonbondTerm",
    "ForceField",
]

class EnergyTerm(object):
    def __init__(self, label, calculate_eg):
        self.label = label
        self.calculate_eg = calculate_eg

    def __call__(self, coordinates, unit_cell):
        raise NotImplementedError

    def init_graph(self, graph):
        raise NotImplementedError


class ValenceTerm(EnergyTerm):
    def __init__(self, label, calculate_eg, match_definition):
        self.match_definition = match_definition
        EnergyTerm.__init__(self, label, calculate_eg)

    def calculate_ic(self, coordinates):
        raise NotImplementedError

    def do_chain_rule(self, coordinates, ic, gscalar):
        raise NotImplementedError

    def init_graph(self, graph):
        # the internal coordinates are enumerated
        indices = []
        for match in MatchGenerator(self.match_definition)(graph):
            indices.append([val for key,val in sorted(match.forward.iteritems())])
        self.indices = numpy.array(indices, int)

    def __call__(self, coordinates, gradient_sum, unit_cell):
        energy_sum = 0.0
        for row in self.indices:
            sub_coordinates = numpy.array([
                unit_cell.shortest_vector(coordinates[index] - coordinates[row[0]])
                for index in row
            ])
            ic, jacobian = self.calculate_ic(sub_coordinates)
            energy, gscalar = self.calculate_eg(ic)
            energy_sum += energy
            for index, gradient_atom in zip(row, gscalar*jacobian):
                gradient_sum[index] += gradient_atom
        return energy_sum


class BondStretchTerm(ValenceTerm):
    def __init__(self, label, calculate_eg, atom_criteria):
        ValenceTerm.__init__(self, label, calculate_eg, BondMatchDefinition(
            [CriteriaSet("foo", dict((index, criterion) for index, criterion in enumerate(atom_criteria)))]
        ))

    def calculate_ic(self, coordinates):
        delta = coordinates[0] - coordinates[1]
        distance = numpy.linalg.norm(delta)
        jacobian = numpy.zeros(coordinates.shape, float)
        jacobian[0] = delta/distance
        jacobian[1] = -jacobian[0]
        return distance, jacobian



class ThreeBodyTerm(ValenceTerm):
    def __init__(self, label, calculate_eg, atom_criteria):
        ValenceTerm.__init__(self, label, calculate_eg, BendingAngleMatchDefinition(
            [CriteriaSet("foo", dict((index, criterion) for index, criterion in enumerate(atom_criteria)))]
        ))


class BendingCosineTerm(ThreeBodyTerm):
    def calculate_ic(self, coordinates):
        jacobian = numpy.zeros(coordinates.shape, float)
        a = coordinates[0] - coordinates[1]
        b = coordinates[2] - coordinates[1]
        dot_ab = numpy.dot(a, b)
        dot_aa = numpy.dot(a, a)
        dot_bb = numpy.dot(b, b)
        root = numpy.sqrt(dot_aa*dot_bb)
        jacobian[0] = (root*b - dot_ab*dot_bb*a/root)/root**2
        jacobian[2] = (root*a - dot_ab*dot_aa*b/root)/root**2
        jacobian[1] = -jacobian[0]-jacobian[2]
        return dot_ab/root, jacobian


class UreyBradleyTerm(ThreeBodyTerm):
    def calculate_ic(self, coordinates):
        delta = coordinates[0] - coordinates[2]
        distance = numpy.linalg.norm(delta)
        jacobian = numpy.zeros(coordinates.shape, float)
        jacobian[0] = delta/distance
        jacobian[2] = -jacobian[0]
        return distance, jacobian


class FourBodyTerm(ValenceTerm):
    def __init__(self, label, calculate_eg, atom_criteria):
        ValenceTerm.__init__(self, label, calculate_eg, DihedralAngleMatchDefinition(
            [CriteriaSet("foo", dict((index, criterion) for index, criterion in enumerate(atom_criteria)))]
        ))


class DihedralCosineTerm(FourBodyTerm):
    def calculate_ic(self, coordinates):
        jacobian = numpy.zeros(coordinates.shape, float)
        a = coordinates[0] - coordinates[1]
        b = coordinates[2] - coordinates[1]
        c = coordinates[3] - coordinates[2]
        dot_ab = numpy.dot(a, b)
        dot_cb = numpy.dot(c, b)
        dot_bb = numpy.dot(b, b)
        d = a - b*dot_ab/dot_bb
        e = c - b*dot_cb/dot_bb
        dot_de = numpy.dot(d, e)
        dot_dd = numpy.dot(d, d)
        dot_ee = numpy.dot(e, e)
        root = numpy.sqrt(dot_dd*dot_ee)
        gd = (root*e - dot_de*dot_ee*d/root)/root**2
        ge = (root*d - dot_de*dot_dd*e/root)/root**2
        jda = numpy.identity(3) - numpy.outer(b,b)/numpy.sqrt(dot_bb)
        jdb = -(
            (numpy.outer(b, a)*dot_bb - 2*numpy.outer(b, b)*dot_ab)/dot_bb/dot_bb +
            numpy.identity(3)*dot_ab/dot_bb
        )
        jec = numpy.identity(3) - numpy.outer(b,b)/numpy.sqrt(dot_bb)
        jeb = -(
            (numpy.outer(b, c)*dot_bb - 2*numpy.outer(b, b)*dot_cb)/dot_bb/dot_bb +
            numpy.identity(3)*dot_cb/dot_bb
        )
        jacobian[0] = numpy.dot(jda, gd)
        jacobian[3] = numpy.dot(jec, ge)
        jacobian[2] = numpy.dot(jdb, gd) + numpy.dot(jeb, ge)
        jacobian[1] = -jacobian[2]
        jacobian[1] -= jacobian[0]
        jacobian[2] -= jacobian[3]
        return dot_de/root, jacobian



class OneFourTerm(FourBodyTerm):
    def calculate_ic(self, coordinates):
        delta = coordinates[0] - coordinates[3]
        distance = numpy.linalg.norm(delta)
        jacobian = numpy.zeros(coordinates.shape, float)
        jacobian[0] = delta/distance
        jacobian[3] = -jacobian[0]
        return distance, jacobian


class NonbondTerm(EnergyTerm):
    def __init__(self, label, calculate_eg, atom_criteria, nonbond_filter, cutoff):
        self.atom0_criterion, self.atom1_criterion = atom_criteria
        self.nonbond_filter = nonbond_filter
        self.cutoff = cutoff
        EnergyTerm.__init__(self, label, calculate_eg)

    def yield_pairs(self, graph):
        self.atom0_criterion.init_graph(graph)
        self.atom1_criterion.init_graph(graph)
        n = len(graph.molecule.numbers)
        for atom0 in xrange(n):
            for atom1 in xrange(n):
                if ((self.atom0_criterion(atom0) and self.atom1_criterion(atom1)) or
                    (self.atom0_criterion(atom1) and self.atom1_criterion(atom0))):
                    distance = graph.get_distance(atom0, atom1)
                    if self.nonbond_filter(distance):
                        #print atom0, atom1
                        yield atom0, atom1

    def init_graph(self, graph):
        graph.init_index()
        graph.init_distances()
        # The nonbonding pairs are enumerated here
        self.indices = numpy.array(list(self.yield_pairs(graph)), int)
        #for row in self.indices:
        #    print [graph.molecule.numbers[index] for index in row]

    def _get_neighbor_cells(self, unit_cell):
        # enumerate all periodic images that might be relevant for the nonbonding
        # interactions in the cutoff radius
        signs = numpy.array([
            [ 1,  1,  1],
            [ 1,  1, -1],
            [ 1, -1,  1],
            [ 1, -1, -1],
            [-1,  1,  1],
            [-1,  1, -1],
            [-1, -1,  1],
            [-1, -1, -1],
        ], int)
        #neighbor_cells = [[0, 0, 0]]
        neighbor_cells = []

        def yield_shell(p):
            for a in xrange(0, p+1):
                for b in xrange(0, p-a+1):
                    c = p-a-b
                    tmp = numpy.array([a,b,c], float)
                    for sign in signs:
                        if ((sign < 0) * (tmp == 0)).any():
                            continue
                        yield tmp*sign


        def pair(la, lb):
            for a in la:
                for b in lb:
                    yield a, b



        found_new = True
        p = 1
        while found_new:
            #print "p", p, len(neighbor_cells)
            found_new = False
            for indices in yield_shell(p):
                #print "indices", indices
                corners_center = numpy.dot(0.5*signs, unit_cell.cell.transpose())
                corners_other = numpy.dot(0.5*signs + indices, unit_cell.cell.transpose())
                for corner_center, corner_other in pair(corners_center, corners_other):
                    distance = numpy.linalg.norm(corner_other - corner_center)
                    if distance < self.cutoff:
                        neighbor_cells.append(indices)
                        found_new = True
                        break
            p += 1
        return numpy.array(neighbor_cells, int)

    def __call__(self, coordinates, gradient_sum, unit_cell):
        energy_sum = 0.0
        energy_cutoff, foo = self.calculate_eg(self.cutoff)
        del foo
        if not unit_cell.cell_active.all():
            raise NotImplementedError("Only 3D periodic cells are supported.")
        neighbor_cells = self._get_neighbor_cells(unit_cell)
        #print neighbor_cells
        for atom0, atom1 in self.indices:
            for neighbor_cell in neighbor_cells:
                delta = coordinates[atom0] - coordinates[atom1] + numpy.dot(unit_cell.cell, neighbor_cell)
                distance = numpy.linalg.norm(delta)
                #print distance, self.cutoff
                if distance < self.cutoff:
                    #print "Neighbor cell", neighbor_cell, atom0, atom1
                    energy, gscalar = self.calculate_eg(distance)
                    direction = delta/distance
                    energy_sum += 0.5*(energy - energy_cutoff)
                    gradient_sum[atom0] += 0.5*gscalar*direction
                    gradient_sum[atom1] -= 0.5*gscalar*direction
        for atom0, atom1 in self.indices:
            if atom0 > atom1:
                delta = coordinates[atom0] - coordinates[atom1]
                distance = numpy.linalg.norm(delta)
                if distance < self.cutoff:
                    #print "Same cell", atom0, atom1
                    energy, gscalar = self.calculate_eg(distance)
                    energy_sum += energy - energy_cutoff
                    direction = delta/distance
                    gradient_sum[atom0] += gscalar*direction
                    gradient_sum[atom1] -= gscalar*direction
        return energy_sum


class ForceField(object):
    def __init__(self, graph, unit_cell, terms):
        self.graph = graph
        self.unit_cell = unit_cell
        self.terms = terms
        for term in terms:
            term.init_graph(graph)

    def __call__(self, coordinates):
        energy_sum = 0.0
        gradient_sum = numpy.zeros(coordinates.shape, float)
        for term in self.terms:
            energy_sum += term(coordinates, gradient_sum, self.unit_cell)
        return energy_sum, gradient_sum

