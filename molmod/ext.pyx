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
"""C extensions for molmod."""


from __future__ import division


import numpy as np
cimport numpy as np

cimport ff
cimport graphs
cimport molecules
cimport similarity
cimport unit_cells


#
#  ff.c
#

def ff_dm_quad(double[:, ::1] cor not None, double[:, ::1] dm0 not None,
               double[:, ::1] dmk not None, double amp, double[:, ::1] gradient not None,
               double[:, ::1] matrix=None, double[:, ::1] reciprocal=None):
    cdef size_t natom = cor.shape[0]
    if cor.shape[1] != 3:
        raise TypeError('cor argument must have three columns.')
    if dm0.shape[0] != natom or dm0.shape[1] != natom:
        raise TypeError('dm0 must have shape (natom, natom).')
    if dmk.shape[0] != natom or dmk.shape[1] != natom:
        raise TypeError('dmk must have shape (natom, natom).')
    if gradient.shape[0] != cor.shape[0] or gradient.shape[1] != cor.shape[1]:
        raise TypeError('gradient must have same shape as cor.')
    if (matrix is None) ^ (reciprocal is None):
        raise TypeError('Either both matrix and reciprocal or given, or both are not given.')
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    if matrix is None:
        return ff.ff_dm_quad(
            natom, long(matrix is not None), &cor[0, 0], &dm0[0, 0], &dmk[0, 0], amp,
            &gradient[0, 0], NULL, NULL)
    else:
        return ff.ff_dm_quad(
            natom, long(matrix is not None), &cor[0, 0], &dm0[0, 0], &dmk[0, 0], amp,
            &gradient[0, 0], &matrix[0, 0], &reciprocal[0, 0])


def ff_dm_reci(double[:, ::1] cor not None, double[::1] radii not None,
               long[:, ::1] dm0 not None, double amp, double[:, ::1] gradient not None,
               double[:, ::1] matrix=None, double[:, ::1] reciprocal=None):
    cdef size_t natom = cor.shape[0]
    if cor.shape[1] != 3:
        raise TypeError('cor argument must have three columns.')
    if radii.shape[0] != natom:
        raise TypeError('radii must have shape (natom,).')
    if dm0.shape[0] != natom or dm0.shape[1] != natom:
        raise TypeError('dm0 must have shape (natom, natom).')
    if gradient.shape[0] != cor.shape[0] or gradient.shape[1] != cor.shape[1]:
        raise TypeError('gradient must have same shape as cor.')
    if (matrix is None) ^ (reciprocal is None):
        raise TypeError('Either both matrix and reciprocal or given, or both are not given.')
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    if matrix is None:
        return ff.ff_dm_reci(
            natom, long(matrix is not None), &cor[0, 0], &radii[0], &dm0[0, 0], amp,
            &gradient[0, 0], NULL, NULL)
    else:
        return ff.ff_dm_reci(
            natom, long(matrix is not None), &cor[0, 0], &radii[0], &dm0[0, 0], amp,
            &gradient[0, 0], &matrix[0, 0], &reciprocal[0, 0])


def ff_bond_quad(double[:, ::1] cor not None, long[:, ::1] pairs not None,
                 double[::1] lengths not None, double amp,
                 double[:, ::1] gradient not None,
                 double[:, ::1] matrix=None, double[:, ::1] reciprocal=None):
    cdef size_t npair = pairs.shape[0]
    if cor.shape[1] != 3:
        raise TypeError('cor argument must have three columns.')
    if pairs.shape[1] != 2:
        raise TypeError('pairs argument must have two columns.')
    if np.asarray(pairs).min() < 0 or np.asarray(pairs).max() >= cor.shape[0]:
        raise ValueError('The pairs array contains atom indexes that are out of bounds.')
    if lengths.shape[0] != npair:
        raise TypeError('lengths must have shape (npair,).')
    if gradient.shape[0] != cor.shape[0] or gradient.shape[1] != cor.shape[1]:
        raise TypeError('gradient must have same shape as cor.')
    if (matrix is None) ^ (reciprocal is None):
        raise TypeError('Either both matrix and reciprocal or given, or both are not given.')
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    if matrix is None:
        return ff.ff_bond_quad(
            npair, long(matrix is not None), &cor[0, 0], &pairs[0, 0], &lengths[0], amp,
            &gradient[0, 0], NULL, NULL)
    else:
        return ff.ff_bond_quad(
            npair, long(matrix is not None), &cor[0, 0], &pairs[0, 0], &lengths[0], amp,
            &gradient[0, 0], &matrix[0, 0], &reciprocal[0, 0])


def ff_bond_hyper(double[:, ::1] cor not None, long[:, ::1] pairs not None,
                  double[::1] lengths not None, double scale, double amp,
                  double[:, ::1] gradient not None,
                  double[:, ::1] matrix=None, double[:, ::1] reciprocal=None):
    cdef size_t npair = pairs.shape[0]
    if cor.shape[1] != 3:
        raise TypeError('cor argument must have three columns.')
    if pairs.shape[1] != 2:
        raise TypeError('pairs argument must have two columns.')
    #if np.asarray(pairs).min() < 0 or np.asarray(pairs).max() >= cor.shape[0]:
    #    raise ValueError('The pairs array contains atom indexes that are out of bounds.')
    if lengths.shape[0] != npair:
        raise TypeError('lengths must have shape (npair,).')
    if gradient.shape[0] != cor.shape[0] or gradient.shape[1] != cor.shape[1]:
        raise TypeError('gradient must have same shape as cor.')
    if (matrix is None) ^ (reciprocal is None):
        raise TypeError('Either both matrix and reciprocal or given, or both are not given.')
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    if matrix is None:
        return ff.ff_bond_hyper(
            npair, long(matrix is not None), &cor[0, 0], &pairs[0, 0], &lengths[0], scale,
            amp, &gradient[0, 0], NULL, NULL)
    else:
        return ff.ff_bond_hyper(
            npair, long(matrix is not None), &cor[0, 0], &pairs[0, 0], &lengths[0], scale,
            amp, &gradient[0, 0], &matrix[0, 0], &reciprocal[0, 0])


#
# graphs.c
#


def graphs_floyd_warshall(long[:, ::1] dm):
    cdef size_t nvertex = dm.shape[0]
    if dm.shape[1] != nvertex:
        raise TypeError('dm must be a square matrix.')
    graphs.graphs_floyd_warshall(nvertex, &dm[0, 0])


#
# molecules.c
#

def molecules_distance_matrix(double[:, ::1] cor not None, double[:, ::1] matrix=None,
                              double[:, ::1] reciprocal=None):
    cdef size_t natom = cor.shape[0]
    if cor.shape[1] != 3:
        raise TypeError('cor argument must have three columns.')
    if (matrix is None) ^ (reciprocal is None):
        raise TypeError('Either both matrix and reciprocal or given, or both are not given.')
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    cdef np.ndarray[double, ndim=2] dm = np.zeros((natom, natom), float)
    if matrix is None:
        molecules.molecules_distance_matrix(
            natom, &cor[0, 0], long(matrix is not None), NULL, NULL, &dm[0, 0])
    else:
        molecules.molecules_distance_matrix(
            natom, &cor[0, 0], long(matrix is not None), &matrix[0, 0], &reciprocal[0, 0],
            &dm[0, 0])
    return dm


#
# similarity.c
#


def similarity_table_labels(const long[::1] labels not None):
    cdef size_t n = labels.shape[0]
    cdef np.ndarray[long, ndim=2] labels_table = np.zeros(((n*(n-1))//2, 2), dtype=int)
    similarity.similarity_table_labels(n, &labels[0], &labels_table[0, 0])
    return labels_table


def similarity_table_distances(const double[:, ::1] distance_matrix not None):
    cdef size_t n = distance_matrix.shape[0]
    if distance_matrix.shape[1] != n:
        raise TypeError('distance_matrix must be sqaure')
    cdef np.ndarray[double] distances_table = np.zeros(((n*(n-1))//2), dtype=float)
    similarity.similarity_table_distances(n, &distance_matrix[0, 0], &distances_table[0])
    return distances_table


def similarity_measure(long[:, ::1] labels1 not None, double[::1] distances1 not None,
                       long[:, ::1] labels2 not None, double[::1] distances2 not None,
                       double margin, double cutoff):
    cdef size_t n1 = labels1.shape[0]
    cdef size_t n2 = labels2.shape[0]
    if labels1.shape[1] != 2:
        raise TypeError('labels1 must have two columns.')
    if distances1.shape[0] != n1:
        raise TypeError('labels1 and distances1 must have the same numbers of rows.')
    if labels2.shape[1] != 2:
        raise TypeError('labels1 must have two columns.')
    if distances2.shape[0] != n2:
        raise TypeError('labels2 and distances2 must have the same numbers of rows.')
    return similarity.similarity_measure(n1, &labels1[0, 0], &distances1[0],
                                         n2, &labels2[0, 0], &distances2[0],
                                         margin, cutoff)


#
# unit_cell.c
#


def unit_cell_get_radius_indexes(double[:, ::1] matrix not None,
                                 double[:, ::1] reciprocal not None,
                                 double radius, long[::1] max_ranges not None,
                                 long[:, ::1] indexes not None):
    if matrix is not None and matrix.shape[0] != 3 and matrix.shape[1] != 3:
        raise TypeError('matrix must be an array with shape (3, 3)')
    if reciprocal is not None and reciprocal.shape[0] != 3 and reciprocal.shape[1] != 3:
        raise TypeError('reciprocal must be an array with shape (3, 3)')
    if radius < 0:
        raise ValueError('Radius must be positive.')
    if max_ranges.shape[0] != 3:
        raise TypeError('max_ranges must have shape (3, ).')
    if indexes.shape[1] != 3:
        raise TypeError('indexes must have three columns.')
    return unit_cells.unit_cell_get_radius_indexes(&matrix[0, 0], &reciprocal[0, 0], radius,
                                                   &max_ranges[0], &indexes[0, 0])
