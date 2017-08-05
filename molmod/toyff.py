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
"""Tools to guess initial geometries quickly based on the molecular graph

   The ToyFF is not meant for accurate geometries, but rather to generate
   molecular geometries from scratch. The routines below start with just a
   random set of coordinates and turn that into a rough molecular geometry.
   Post-processing with a more reliable force field is mandatory.
"""


import numpy as np
import pkg_resources

from molmod.molecules import Molecule
from molmod.periodic import periodic
from molmod.ext import ff_dm_quad, ff_dm_reci, ff_bond_quad, ff_bond_hyper


__all__ = ["guess_geometry", "tune_geometry", "ToyFF", "SpecialAngles"]


def guess_geometry(graph, unit_cell=None, verbose=False):
    """Construct a molecular geometry based on a molecular graph.

       This routine does not require initial coordinates and will give a very
       rough picture of the initial geometry. Do not expect all details to be
       in perfect condition. A subsequent optimization with a more accurate
       level of theory is at least advisable.

       Argument:
        | ``graph``  --  The molecular graph of the system, see
                         :class:molmod.molecular_graphs.MolecularGraph

       Optional argument:
        | ``unit_cell``  --  periodic boundry conditions, see
                             :class:`molmod.unit_cells.UnitCell`
        | ``verbose``  --  Show optimizer progress when True
    """

    N = len(graph.numbers)
    from molmod.minimizer import Minimizer, ConjugateGradient, \
        NewtonLineSearch, ConvergenceCondition, StopLossCondition

    search_direction = ConjugateGradient()
    line_search = NewtonLineSearch()
    convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6)
    stop_loss = StopLossCondition(max_iter=500, fun_margin=0.1)

    ff = ToyFF(graph, unit_cell)
    x_init = np.random.normal(0, 1, N*3)

    #  level 1 geometry optimization: graph based
    ff.dm_quad = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    #  level 2 geometry optimization: graph based + pauli repulsion
    ff.dm_quad = 1.0
    ff.dm_reci = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    # Add a little noise to avoid saddle points
    x_init += np.random.uniform(-0.01, 0.01, len(x_init))

    #  level 3 geometry optimization: bond lengths + pauli
    ff.dm_quad = 0.0
    ff.dm_reci = 0.2
    ff.bond_quad = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    #  level 4 geometry optimization: bond lengths + bending angles + pauli
    ff.bond_quad = 0.0
    ff.bond_hyper = 1.0
    ff.span_quad = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    x_opt = x_init

    mol = Molecule(graph.numbers, x_opt.reshape((N, 3)))
    return mol


def tune_geometry(graph, mol, unit_cell=None, verbose=False):
    """Fine tune a molecular geometry, starting from a (very) poor guess of
       the initial geometry.

       Do not expect all details to be in perfect condition. A subsequent
       optimization with a more accurate level of theory is at least advisable.

       Arguments:
        | ``graph``  --  The molecular graph of the system, see
                         :class:molmod.molecular_graphs.MolecularGraph
        | ``mol``  --  A :class:molmod.molecules.Molecule class with the initial
                       guess of the coordinates

       Optional argument:
        | ``unit_cell``  --  periodic boundry conditions, see
                             :class:`molmod.unit_cells.UnitCell`
        | ``verbose``  --  Show optimizer progress when True
    """

    N = len(graph.numbers)
    from molmod.minimizer import Minimizer, ConjugateGradient, \
        NewtonLineSearch, ConvergenceCondition, StopLossCondition

    search_direction = ConjugateGradient()
    line_search = NewtonLineSearch()
    convergence = ConvergenceCondition(grad_rms=1e-6, step_rms=1e-6)
    stop_loss = StopLossCondition(max_iter=500, fun_margin=1.0)

    ff = ToyFF(graph, unit_cell)
    x_init = mol.coordinates.ravel()

    #  level 3 geometry optimization: bond lengths + pauli
    ff.dm_reci = 0.2
    ff.bond_quad = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    #  level 4 geometry optimization: bond lengths + bending angles + pauli
    ff.bond_quad = 0.0
    ff.bond_hyper = 1.0
    ff.span_quad = 1.0
    minimizer = Minimizer(x_init, ff, search_direction, line_search, convergence, stop_loss, anagrad=True, verbose=verbose)
    x_init = minimizer.x

    x_opt = x_init

    mol = Molecule(graph.numbers, x_opt.reshape((N, 3)))
    return mol


class ToyFF(object):
    """A force field implementation for generating geometries.

       See :func:guess_geomtry and :func:tune_geomtry for two practical use
       cases.
    """

    def __init__(self, graph, unit_cell=None):
        """
           Argument:
            | ``graph``  --  the molecular graph from which the force field terms
                             are extracted. See
                             :class:molmod.molecular_graphs.MolecularGraph

           Optional argument:
            | ``unit_cell``  --  periodic boundry conditions, see
                                 :class:`molmod.unit_cells.UnitCell`
        """
        from molmod.bonds import bonds

        if unit_cell is None:
            self.matrix = None
            self.reciprocal = None
        else:
            self.matrix = unit_cell.matrix
            self.reciprocal = unit_cell.reciprocal

        self.dm = graph.distances.astype(int)
        dm = self.dm.astype(float)
        self.dm0 = dm**2
        self.dmk = (dm+0.1)**(-3)
        self.vdw_radii = np.array([periodic[number].vdw_radius for number in graph.numbers], dtype=float)
        self.covalent_radii = np.array([periodic[number].covalent_radius for number in graph.numbers], dtype=float)

        bond_edges = []
        bond_lengths = []
        for i, j in graph.edges:
            bond_edges.append((i, j))
            bond_lengths.append(bonds.get_length(graph.numbers[i], graph.numbers[j]))
        self.bond_edges = np.array(bond_edges, int)
        self.bond_lengths = np.array(bond_lengths, float)

        special_angles = SpecialAngles()

        span_edges = []
        span_lengths = []
        for i, neighbors in graph.neighbors.items():
            number_i = graph.numbers[i]
            if (number_i >= 5 and number_i <=8):
                valence = len(neighbors) + abs(number_i-6)
            elif number_i >= 13 and number_i <= 16:
                valence = len(neighbors) + abs(number_i-14)
            else:
                valence = -1
            if valence < 2 or valence > 6:
                default_angle = np.pi/180.0*115.0
            elif valence == 2:
                default_angle = np.pi
            elif valence == 3:
                default_angle = np.pi/180.0*125.0
            elif valence == 4:
                default_angle = np.pi/180.0*109.0
            elif valence == 5:
                default_angle = np.pi/180.0*100.0
            elif valence == 6:
                default_angle = np.pi/180.0*90.0
            for j in neighbors:
                number_j = graph.numbers[j]
                for k in neighbors:
                    if j < k and not frozenset([j, k]) in graph.edges:
                        number_k = graph.numbers[k]

                        triplet = (
                            number_j, len(graph.neighbors[j]),
                            number_i, len(graph.neighbors[i]),
                            number_k, len(graph.neighbors[k]),
                        )

                        angle = special_angles.get_angle(triplet)
                        if angle is None:
                            angle = default_angle

                        dj = bonds.get_length(number_i, number_j)
                        dk = bonds.get_length(number_i, number_k)
                        d = np.sqrt(dj**2+dk**2-2*dj*dk*np.cos(angle))
                        span_edges.append((j, k))
                        span_lengths.append(d)
        self.span_edges = np.array(span_edges, int)
        self.span_lengths = np.array(span_lengths, float)

        self.dm_quad = 0.0
        self.dm_reci = 0.0
        self.bond_quad = 0.0
        self.span_quad = 0.0
        self.bond_hyper = 0.0
        self.bond_hyper_scale = 5.0

    def __call__(self, x, do_gradient=False):
        """Compute the energy (and gradient) for a set of Cartesian coordinates

           Argument:
            | ``x``  --  the Cartesian coordinates
            | ``do_gradient``  --  when set to True, the gradient is also
                                   computed and returned. [default=False]
        """
        x = x.reshape((-1, 3))
        result = 0.0

        gradient = np.zeros(x.shape, float)
        if self.dm_quad > 0.0:
            result += ff_dm_quad(x, self.dm0, self.dmk, self.dm_quad,
                                 gradient, self.matrix, self.reciprocal)
        if self.dm_reci:
            result += ff_dm_reci(x, self.vdw_radii, self.dm, self.dm_reci,
                                 gradient, self.matrix, self.reciprocal)
        if self.bond_quad:
            result += ff_bond_quad(x, self.bond_edges, self.bond_lengths,
                                   self.bond_quad, gradient, self.matrix,
                                   self.reciprocal)
        if self.span_quad:
            result += ff_bond_quad(x, self.span_edges, self.span_lengths,
                                   self.span_quad, gradient, self.matrix,
                                   self.reciprocal)
        if self.bond_hyper:
            result += ff_bond_hyper(x, self.bond_edges, self.bond_lengths,
                                    self.bond_hyper_scale, self.bond_hyper,
                                    gradient, self.matrix, self.reciprocal)

        if do_gradient:
            return result, gradient.ravel()
        else:
            return result


class SpecialAngles(object):
    """A database with precomputed valence angles from small molecules"""
    def __init__(self):
        self._angle_dict = {}
        with pkg_resources.resource_stream(__name__, 'data/toyff_angles.txt') as f:
            for line in f:
                if line[0] != '#':
                    key = tuple(int(word) for word in line[0:line.index(b':')].split(b","))
                    value = np.pi/180.0*float(line[line.index(b':')+1:-1])
                    self._angle_dict[key] = value

    def get_angle(self, triplet):
        """Get a rest angle for a given triplet

           A triplet consists of a tuple with six elements: (n0, v0, n1, v1, n2, v2)
           The indexes refer to consecutive atoms forming a valence angle. n0,
           n1 and n2 are the atom numbers of the angle and v0, v1 and v2 are the
           valences of the corresponding atoms. n1 and v1 are the values for the
           central atom in the angle.
        """
        return self._angle_dict.get(triplet)
