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


from molmod.molecules import Molecule
from molmod.data.periodic import periodic

from molmod.ext import ff_dm_quad, ff_dm_reci, ff_bond_quad, ff_bond_hyper

import numpy


__all__ = ["guess_geometry", "tune_geometry", "ToyFF"]


def guess_geometry(graph):
    """Construct a molecular geometry based on a molecular graph.

       This routine does not require initial coordinates and will give a very
       rough picture of the initial geometry. Do not expect all details to be
       in perfect condition. A subsequent optimization with a more accurate
       level of theory is at least advisable.

       Arguments:
         graph  --  The molecular graph of the system
    """

    N = len(graph.numbers)
    from molmod.minimizer import Minimizer, NewtonGLineSearch

    ff = ToyFF(graph)
    x_init = numpy.random.normal(0,1,N*3)

    #  level 1 geometry optimization: graph based
    ff.dm_quad = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-10, 1e-8, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    #  level 2 geometry optimization: graph based + pauli repulsion
    ff.dm_quad = 1.0
    ff.dm_reci = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-10, 1e-8, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    # Add a little noise to avoid saddle points
    x_init += numpy.random.uniform(-0.01, 0.01, len(x_init))

    #  level 3 geometry optimization: bond lengths + pauli
    ff.dm_quad = 0.0
    ff.dm_reci = 0.2
    ff.bond_quad = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-3, 1e-3, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    #  level 4 geometry optimization: bond lengths + bending angles + pauli
    ff.bond_quad = 0.0
    ff.bond_hyper = 1.0
    ff.span_quad = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-6, 1e-6, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    x_opt = x_init

    mol = Molecule(graph.numbers, x_opt.reshape((N,3)))
    return mol


def tune_geometry(graph, mol):
    """Fine tune a molecular geometry, starting from a (very) poor guess of
       the initial geometry.

       Do not expect all details to be in perfect condition. A subsequent
       optimization with a more accurate level of theory is at least advisable.

       Arguments:
         graph  --  The molecular graph of the system
         mol  --  The initial guess of the coordinates
    """

    N = len(graph.numbers)
    from molmod.minimizer import Minimizer, NewtonGLineSearch

    ff = ToyFF(graph)
    x_init = mol.coordinates.ravel()

    #  level 3 geometry optimization: bond lengths + pauli
    ff.dm_reci = 0.2
    ff.bond_quad = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-3, 1e-3, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    #  level 4 geometry optimization: bond lengths + bending angles + pauli
    ff.bond_quad = 0.0
    ff.bond_hyper = 1.0
    ff.span_quad = 1.0
    minimizer = Minimizer(x_init, ff, NewtonGLineSearch, 1e-6, 1e-6, 2*N, 500, 50, do_gradient=True, verbose=False)
    x_init = minimizer.x

    x_opt = x_init

    mol = Molecule(graph.numbers, x_opt.reshape((N,3)))
    return mol


class ToyFF(object):
    """A force field implementation for generating geometries.

       See guess_geomtry and tune_geomtry for two practical use cases.
    """

    def __init__(self, graph):
        from molmod.data.bonds import bonds

        self.dm = graph.distances.astype(numpy.int32)
        dm = self.dm.astype(float)
        self.dm0 = dm**2
        self.dmk = dm**(-3)
        self.vdw_radii = numpy.array([periodic[number].vdw_radius for number in graph.numbers], dtype=float)
        self.covalent_radii = numpy.array([periodic[number].covalent_radius for number in graph.numbers], dtype=float)

        num_bonds = len(graph.pairs)
        bond_pairs = []
        bond_lengths = []
        for counter, (i,j) in enumerate(graph.pairs):
            bond_pairs.append((i,j))
            bond_lengths.append(bonds.get_length(graph.numbers[i],graph.numbers[j]))
        self.bond_pairs = numpy.array(bond_pairs, numpy.int32)
        self.bond_lengths = numpy.array(bond_lengths, float)

        span_pairs = []
        span_lengths = []
        for i, neighbors in graph.neighbors.iteritems():
            number_i = graph.numbers[i]
            if (number_i >= 5 and number_i <=8):
                order = len(neighbors) + abs(number_i-6)
            elif number_i >= 13 and number_i <= 16:
                order = len(neighbors) + abs(number_i-14)
            else:
                order = -1
            if order < 2 or order > 6:
                angle = numpy.pi/180.0*115.0
            elif order == 2:
                angle = numpy.pi
            elif order == 3:
                angle = numpy.pi/180.0*125.0
            elif order == 4:
                angle = numpy.pi/180.0*109.0
            elif order == 5:
                angle = numpy.pi/180.0*100.0
            elif order == 6:
                angle = numpy.pi/180.0*90.0
            for j in neighbors:
                number_j = graph.numbers[j]
                for k in neighbors:
                    if j<k and not frozenset([j,k]) in graph.pairs:
                        number_k = graph.numbers[k]
                        dj = bonds.get_length(number_i,number_j)
                        dk = bonds.get_length(number_i,number_k)
                        d = numpy.sqrt(dj**2+dk**2-2*dj*dk*numpy.cos(angle))
                        span_pairs.append((j,k))
                        span_lengths.append(d)
        self.span_pairs = numpy.array(span_pairs, numpy.int32)
        self.span_lengths = numpy.array(span_lengths, float)

        self.dm_quad = 0.0
        self.dm_reci = 0.0
        self.bond_quad = 0.0
        self.span_quad = 0.0
        self.bond_hyper = 0.0

    def __call__(self, x, do_gradient=False):
        x = x.reshape((-1,3))
        result = 0.0

        gradient = numpy.zeros(x.shape, float)
        if self.dm_quad > 0.0:
            result += ff_dm_quad(x, self.dm0, self.dmk, self.dm_quad, gradient)
        if self.dm_reci:
            result += ff_dm_reci(1.0*self.vdw_radii, x, self.dm, self.dm_reci, gradient)
        if self.bond_quad:
            result += ff_bond_quad(x, self.bond_pairs, self.bond_lengths, self.bond_quad, gradient)
        if self.span_quad:
            result += ff_bond_quad(x, self.span_pairs, self.span_lengths, self.span_quad, gradient)
        if self.bond_hyper:
            result += ff_bond_hyper(x, self.bond_pairs, self.bond_lengths, 5.0, self.bond_hyper, gradient)

        if do_gradient:
            return result, gradient.ravel()
        else:
            return result

