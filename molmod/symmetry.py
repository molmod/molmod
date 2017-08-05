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
"""Tools to analyze the symmetry of molecules"""


from molmod.units import angstrom
from molmod.transformations import fit_rmsd


__all__ = ["compute_rotsym"]


def compute_rotsym(molecule, graph, threshold=1e-3*angstrom):
    """Compute the rotational symmetry number

       Arguments:
        | ``molecule``  --  The molecule
        | ``graph``  --  The corresponding bond graph

       Optional argument:
        | ``threshold``  --  only when a rotation results in an rmsd below the
                             given threshold, the rotation is considered to
                             transform the molecule onto itself.
    """
    result = 0
    for match in graph.symmetries:
        permutation = list(j for i,j in sorted(match.forward.items()))
        new_coordinates = molecule.coordinates[permutation]
        rmsd = fit_rmsd(molecule.coordinates, new_coordinates)[2]
        if rmsd < threshold:
            result += 1
    return result
