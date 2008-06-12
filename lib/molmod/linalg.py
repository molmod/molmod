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


__all__ = ["safe_solve", "extended_solve"]


def safe_solve(A, b, cutoff=0.0):
    U, W, Vt = numpy.linalg.svd(A, full_matrices=False)
    rank = sum(abs(W)>(max(abs(W))*cutoff))
    W = W[:rank]
    Ut = U[:,:rank].transpose()
    V = Vt[:rank].transpose()
    return numpy.dot(V, (numpy.dot(Ut, b)/W))


def extended_solve(A, b, cutoff=0.0):
    U, W, Vt = numpy.linalg.svd(A, full_matrices=True)
    rank = sum(abs(W)>(max(abs(W))*cutoff))
    W = W[:rank]
    Ut = U[:,:rank].transpose()
    nullspace = Vt[rank:].transpose()
    V = Vt[:rank].transpose()
    solution = numpy.dot(V, (numpy.dot(Ut, b)/W))
    return solution, nullspace




