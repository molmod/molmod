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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


import math, numpy

from molmod.molecular_graphs import MolecularGraph


def lone_pair_2(bond1, bond2, angle):
    """
    Returns two vectors of unit length that point in the direction of the
    lone pairs. (The atom under study is supposed to be of type O.)

    Arguments:
    bond1 -- the relative vector from the central atom to the first
             bonded atom
    bond2 -- ...
    angle -- the angle between the two lone pairs.

    Returns:
    lone1 -- the relative vector with unit length pointing along the
             first lone pair
    lone2 -- ...
    """

    in_plane = bond1/math.sqrt(numpy.dot(bond1, bond1)) + bond2/math.sqrt(numpy.dot(bond2, bond2))
    length = math.sqrt(numpy.dot(in_plane, in_plane))
    assert length > 0, "The two bonds opposite."
    in_plane /= -length

    ortho = numpy.cross(bond1, bond2)
    length = math.sqrt(numpy.dot(ortho, ortho))
    assert length > 0, "The two bonds are parallel."
    ortho /= length

    c = math.cos(0.5*angle)
    s = math.sin(0.5*angle)
    lone1 = c*in_plane + s*ortho
    lone2 = c*in_plane - s*ortho

    return lone1, lone2


def lone_pair_1(bond1, bond2, bond3):
    """
    Returns a vector pointing in the direction of the lone pair with
    unit length. (The atom under study is supposed to be of type N.)

    Arguments:
    bond1 -- the relative vector from the central atom to the first
             bonded atom
    bond2 -- ...
    bond3 -- ...

    Returns:
    lone -- the relative vector with unit length pointing along the
            lone pair
    """

    lone = -(
        bond1/math.sqrt(numpy.dot(bond1, bond1)) +
        bond2/math.sqrt(numpy.dot(bond2, bond2)) +
        bond3/math.sqrt(numpy.dot(bond3, bond3))
    )
    length = math.sqrt(numpy.dot(lone, lone))
    assert length > 0, "The three bonds sum to zero."
    lone /= length
    return lone


def all_lone_pairs(molecule, singles=[7], doubles=[8], angle=1.910):
    """
    Returns a list with pairs (index, lone), where index indicates the
    atom and lone is a relative vector with unit length pointing along
    a lone pair on that atom.

    Arguments:
    molecule -- A molecule for which the lone pairs should be calculated.
    singles -- A list of atom number which should be treated like nitrogen.
    doubles -- A list of atom numbers which should be treated like oxygen.
    angle -- The angle between two lone pairs on the same atom. (in rad)

    Returns:
    lone_pairs -- a list with pairs (index, lone)
    """

    result = []
    mgraph = MolecularGraph(molecule)

    for index, (number, coordinate) in enumerate(zip(molecule.numbers, molecule.coordinates)):
        if number in singles:
            neighbors = mgraph.neighbors[index]
            result.append((
                index,
                lone_pair_1(
                    molecule.coordinates[neighbors[0]] - coordinate,
                    molecule.coordinates[neighbors[1]] - coordinate,
                    molecule.coordinates[neighbors[2]] - coordinate
                )
            ))
        elif number in doubles:
            neighbors = mgraph.neighbors[index]
            lone1, lone2 = lone_pair_2(
                molecule.coordinates[neighbors[0]] - coordinate,
                molecule.coordinates[neighbors[1]] - coordinate,
                angle
            )
            result.append((index, lone1))
            result.append((index, lone2))

    return result



