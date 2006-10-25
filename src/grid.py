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

__all__ = ["yield_3d_indices", "yield_cusp_grid", "hierarchical_grid"]


def yield_3d_indices(num):
    for i in xrange(num):
        for j in xrange(num):
            for k in xrange(num):
                yield numpy.array([i, j, k], int)


def yield_cusp_grid(center, corner, spacing, size, div, subdiv, level):
    spacing /= div
    subspacing = spacing/subdiv
    center_index = numpy.floor((center - corner)/subspacing).astype(int)
    #print "center_index", center_index

    for index in yield_3d_indices(size*div*subdiv):
        if level == 1 or (abs(index - center_index) > (size/2)*subdiv).any():
            #print index, point, spacing**3
            yield (index+0.5) * subspacing + corner, subspacing**3

    if level > 1:
        new_corner = (center_index - (size/2)*subdiv)*subspacing + corner
        for data in yield_cusp_grid(center, new_corner, spacing, size, div, subdiv, level-1):
            yield data


def hierarchical_grid(coordinates, callback, div=2, min_spacing=1e-4, max_spacing=1e-0, margin=8, alpha1=0.1, alpha2=0.3, threshold=1):
    from molmod.helpers import gridf
    center = sum(coordinates)/len(coordinates)
    size = (coordinates.max(axis=0) - coordinates.min(axis=0)).max() + 2*margin
    gridf.recursive_grid(center, size, div, min_spacing, max_spacing, alpha1, alpha2, threshold, coordinates, callback)
