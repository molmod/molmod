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

__all__ = ["hierarchical_grid"]


def hierarchical_grid(coordinates, callback, div=2, min_spacing=1e-4, max_spacing=1e-0, margin=8, alpha1=0.1, alpha2=0.3, threshold=1):
    from molmod.helpers import gridf
    center = sum(coordinates)/len(coordinates)
    size = (coordinates.max(axis=0) - coordinates.min(axis=0)).max() + 2*margin
    gridf.recursive_grid(center, size, div, min_spacing, max_spacing, alpha1, alpha2, threshold, coordinates, callback)


class Grid(object):
    def __init__(self, points, volumes=None, data={}):
        self.__dict__.update(data)
        self.points = points
        self.volumes = volumes
