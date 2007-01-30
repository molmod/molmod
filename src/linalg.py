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


import numpy


def generic_solve(A, b, cutoff=0.0):
    V, S, Wt = numpy.linalg.svd(A, True)
    print S
    rank = sum(abs(S)>(max(abs(S))*cutoff))
    S = S[:rank]
    V = V[:,:rank]
    nullspace = Wt[rank:].transpose()
    Wt = Wt[:rank]
    particular = numpy.dot(Wt.transpose(), (numpy.dot(V.transpose(), b)/S))
    error = ((numpy.dot(A, particular) - b)**2).sum()

    return particular, nullspace, error
