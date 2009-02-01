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


from molmod.molecular_graphs import MolecularGraph
from molmod.molecules import Molecule
from molmod.minimizer import Minimizer, GoldenLineSearch, NewtonGLineSearch
from molmod.units import angstrom, deg
import molmod.ic as ic

import numpy


__all__ = [
    "Patch", "PatchError", "Restraint", "MethylPatch", "AmonoAcidPatch"
]


class Patch(object):
    def execute(self, mol, graph, indexes, **parameters):
        if mol.size != len(graph.numbers):
            raise PatchError("Molecule and graph do not match.")
        if (mol.numbers != graph.numbers).all():
            raise PatchError("Molecule and graph do not match.")

        self.check(graph, indexes)
        selection = self.select(graph, indexes)
        remove_indexes = self.get_remove(graph, selection)
        add_pairs, add_numbers = self.get_add(graph, selection)

        # map of old indexes to new indexes
        old_to_new = {}
        j = 0
        last = 0
        for i in xrange(mol.size):
            if i == remove_indexes[j]:
                j += 1
            else:
                old_to_new[i] = last
                last += 1
        # create a new graph
        new_pairs = []
        for i,j in graph.pairs:
            i = old_to_new.get(i)
            if i is None: continue
            j = old_to_new.get(j)
            if j is None: continue
            new_pairs.append((i,j))
        for i,j in add_pairs:
            if i < 0:
                i = len(old_to_new)-i-1
            else:
                i = old_to_new[i]
            if j < 0:
                j = len(old_to_new)-j-1
            else:
                j = old_to_new[j]
            new_pairs.append((i,j))
        new_numbers = list(graph.numbers)
        for remove_index in remove_indexes:
            del new_numbers[remove_index]
        new_numbers.extend(add_numbers)
        new_graph = MolecularGraph(new_pairs, new_numbers)
        # create the fixed atoms
        x_fixed = {}
        for i,j in selection.iteritems():
            new_i = old_to_new.get(i)
            if new_i is not None:
                x_fixed[j] = mol.coordinates[i]
        x_fixed = list(cor for i,cor in sorted(x_fixed.iteritems()))
        num_fixed = len(x_fixed)

        # prepare the optimization
        self.update_restraints(graph, selection, parameters)
        x_init = numpy.concatenate([
            numpy.random.uniform(-0.5,0.5,3)
            for add_number in add_numbers
        ])

        # the function to minimize
        def cost(x, do_gradient=False):
            x = x.reshape((-1,3))
            result = 0.0
            if do_gradient:
                gradient = numpy.zeros(x.shape, float)
            for restraint in self.restraints:
                args = []
                for index in restraint.indexes:
                    if index < num_fixed:
                        args.append(x_fixed[index])
                    else:
                        args.append(x[index-num_fixed])
                if do_gradient:
                    args.append(1)
                    q, gq = restraint.fn(*args)
                    f = 2*restraint.k*(q-restraint.q0)
                    for j, index in enumerate(restraint.indexes):
                        if index >= num_fixed:
                            gradient[index-num_fixed] += f*gq[j]
                else:
                    q, = restraint.fn(*args)

                result += restraint.k*(q-restraint.q0)**2
            if do_gradient:
                return result, gradient.ravel()
            else:
                return result

        # do the actual minimization
        minimizer = Minimizer(
            x_init, cost, NewtonGLineSearch, 1e-6, 1e-10, 1e1, 400, 50,
            do_gradient=True, absftol=True, verbose=False)
        x_opt = minimizer.x

        # create a new molecule
        new_coordinates = list(mol.coordinates)
        for remove_index in remove_indexes:
            del new_coordinates[remove_index]
        for i in xrange(len(add_numbers)):
            new_coordinates.append(x_opt[i*3:(i+1)*3])
        new_mol = Molecule(new_numbers, new_coordinates)
        new_mol.remove_indexes = remove_indexes
        new_mol.add_numbers = add_numbers
        new_mol.add_pairs = new_graph.pairs[-len(add_pairs):]
        new_mol.add_coordinates = new_mol.coordinates[-len(add_numbers):]
        new_mol.old_to_new = old_to_new
        return new_mol, new_graph

    def check(self, graph, indexes):
        """Make sure that the indexes in the mol/graph are compatible with the fragment."""
        raise NotImplementedError

    def select(self, graph, indexes):
        """Assign indexes to the relevant atoms in the molecule"""
        raise NotImplementedError

    def get_remove(self, graph, selection):
        """Return the atom indexes to be removed"""
        raise NotImplementedError

    def get_add(self, graph, selection):
        """Return the added pairs and atom numbers"""
        raise NotImplementedError

    def update_restraints(self, graph, selection, parameters):
        """Update some rest values in the restraints"""
        raise NotImplementedError


class PatchError(Exception):
    pass


class Restraint(object):
    def __init__(self, fn, indexes, k, q0):
        self.fn = fn
        self.indexes = indexes
        self.k = k
        self.q0 = q0


class MethylPatch(Patch):
    def __init__(self):
        self.restraints = [
            Restraint(ic.bond_length, (1,2), 10.0, 1.52*angstrom),
            Restraint(ic.bond_length, (2,3), 10.0, 1.10*angstrom),
            Restraint(ic.bond_length, (2,4), 10.0, 1.10*angstrom),
            Restraint(ic.bond_length, (2,5), 10.0, 1.10*angstrom),
            Restraint(ic.bend_cos, (0,1,2), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (1,2,3), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (1,2,4), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (1,2,5), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (3,2,4), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (4,2,5), 2.0, -1.0/3),
            Restraint(ic.bend_cos, (5,2,3), 2.0, -1.0/3),
            Restraint(ic.dihed_cos, (0,1,2,3), 1.0, -1.0),
        ]

    def check(self, graph, indexes):
        l = len(indexes)
        if l != 1:
            raise PatchError("Expecting only one anchor atom, %i given." % l)
        index = indexes[0]
        l = len(graph.neighbors[index])
        if l != 1:
            raise PatchError("The anchor atom can only have one neighbor, but has %i." % l)

    def select(self, graph, indexes):
        index = indexes[0]
        n, = graph.neighbors[index]
        nns = graph.neighbors[n]
        nn = None
        for nnc in nns:
            if nn is None or graph.numbers[nn] < graph.numbers[nnc]:
                nn = nnc
        return { # make sure that the removed atom is last
            0: nn, 1: n, 2: index
        }

    def get_remove(self, graph, selection):
        return [selection[2]]

    def get_add(self, graph, selection):
        numbers = [6,1,1,1]
        pairs = [(selection[1], -1), (-1, -2), (-1, -3), (-1, -4)]
        return pairs, numbers

    def update_restraints(self, graph, selection, parameters):
        angle = parameters.get("theta")
        if angle is None:
            angle = 109.45*deg
        self.restraints[4].q0 = numpy.cos(angle)
        bond = parameters.get("bond")
        if bond is None:
            bond = 1.52*angstrom
        self.restraints[0].q0 = bond


class AmonoAcidPatch(Patch):
    pass


