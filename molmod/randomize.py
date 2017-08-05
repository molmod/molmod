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
"""Tools for randomizing molecular geometries with reasonable distortions"""


from __future__ import print_function

from builtins import range
from random import shuffle, sample
import copy

import numpy as np

from molmod.molecules import Molecule
from molmod.graphs import GraphError
from molmod.transformations import Translation, Complete
from molmod.vectors import random_orthonormal, random_unit


__all__ = [
    "MolecularDistortion",
    "RandomManipulation", "RandomStretch", "RandomTorsion",
    "RandomBend", "RandomDoubleStretch",
    "iter_halfs_bond", "iter_halfs_bend", "iter_halfs_double",
    "generate_manipulations", "check_nonbond",
    "randomize_molecule", "randomize_molecule_low",
    "single_random_manipulation", "single_random_manipulation_low",
    "random_dimer",
]


class MolecularDistortion(object):
    """A geometeric manipulation (rotation + translation) of a part of molecule

       The data structure also comes with a straight forward human readable
       file format, which makes it easy to save the distortions for later
       reference.
    """
    @classmethod
    def read_from_file(cls, filename):
        """Construct a MolecularDistortion object from a file"""
        with open(filename) as f:
            lines = list(line for line in f if line[0] != '#')
        r = []
        t = []
        for line in lines[:3]:
            values = list(float(word) for word in line.split())
            r.append(values[:3])
            t.append(values[3])
        transformation = Complete(r, t)
        affected_atoms = set(int(word) for word in lines[3].split())
        return cls(affected_atoms, transformation)

    def __init__(self, affected_atoms, transformation):
        """Initialize a new MolecularDistortion object

           Arguments:
             affected_atoms  --  a list of atoms that undergo the transformation
             transformation  --  a transformation object
        """
        self.affected_atoms = affected_atoms
        self.transformation = Complete.cast(transformation)

    def apply(self, coordinates):
        """Apply this distortion to Cartesian coordinates"""
        for i in self.affected_atoms:
            coordinates[i] = self.transformation*coordinates[i]

    def write_to_file(self, filename):
        """Write the object to a file"""
        r = self.transformation.r
        t = self.transformation.t
        with open(filename, "w") as f:
            print("# A (random) transformation of a part of a molecule:", file=f)
            print("# The translation vector is in atomic units.", file=f)
            print("#     Rx             Ry             Rz              T", file=f)
            print("% 15.9e % 15.9e % 15.9e % 15.9e" % (r[0, 0], r[0, 1], r[0, 2], t[0]), file=f)
            print("% 15.9e % 15.9e % 15.9e % 15.9e" % (r[1, 0], r[1, 1], r[1, 2], t[1]), file=f)
            print("% 15.9e % 15.9e % 15.9e % 15.9e" % (r[2, 0], r[2, 1], r[2, 2], t[2]), file=f)
            print("# The indexes of the affected atoms:", file=f)
            print(" ".join(str(i) for i in self.affected_atoms), file=f)



class RandomManipulation(object):
    """Base class for a random manipulation"""
    num_hinge_atoms = None

    def __init__(self, affected_atoms, max_amplitude, hinge_atoms):
        """Initialize a RandomManipulation object

           Arguments:
             affected_atoms  --  a list of atoms that undergo the transformation
             max_amplitude  --  the maximum displacement (unit depends on
                                actual implementation)
             hinge_atoms  --  atoms that are invariant under the transformation
        """
        if len(hinge_atoms) != self.num_hinge_atoms:
            raise ValueError("The number of hinge atoms must be %i, got %i." % (
                self.num_hinge_atoms,
                len(hinge_atoms)
            ))
        self.affected_atoms = affected_atoms
        self.max_amplitude = max_amplitude
        self.hinge_atoms = hinge_atoms

    def apply(self, coordinates):
        """Generate, apply and return a random manipulation"""
        transform = self.get_transformation(coordinates)
        result = MolecularDistortion(self.affected_atoms, transform)
        result.apply(coordinates)
        return result

    def get_transformation(self, coordinates):
        """Construct a transformation object"""
        raise NotImplementedError


class RandomStretch(RandomManipulation):
    """A random variation in a bond length by displacing a part of a molecule"""
    num_hinge_atoms = 2

    def get_transformation(self, coordinates):
        """Construct a transformation object"""
        atom1, atom2 = self.hinge_atoms
        direction = coordinates[atom1] - coordinates[atom2]
        direction /= np.linalg.norm(direction)
        direction *= np.random.uniform(-self.max_amplitude, self.max_amplitude)
        result = Translation(direction)
        return result


class RandomTorsion(RandomManipulation):
    """A random bond torsion by rotation a part of a molecule"""
    num_hinge_atoms = 2

    def get_transformation(self, coordinates):
        """Construct a transformation object"""
        atom1, atom2 = self.hinge_atoms
        center = coordinates[atom1]
        axis = coordinates[atom1] - coordinates[atom2]
        axis /= np.linalg.norm(axis)
        angle = np.random.uniform(-self.max_amplitude, self.max_amplitude)
        return Complete.about_axis(center, angle, axis)


class RandomBend(RandomManipulation):
    """A random bend by rotation a part of a molecule"""
    num_hinge_atoms = 3

    def get_transformation(self, coordinates):
        """Construct a transformation object"""
        atom1, atom2, atom3 = self.hinge_atoms
        center = coordinates[atom2]
        a = coordinates[atom1] - coordinates[atom2]
        b = coordinates[atom3] - coordinates[atom2]
        axis = np.cross(a, b)
        norm = np.linalg.norm(axis)
        if norm < 1e-5:
            # We suppose that atom3 is part of the affected atoms
            axis = random_orthonormal(a)
        else:
            axis /= np.linalg.norm(axis)
        angle = np.random.uniform(-self.max_amplitude, self.max_amplitude)
        return Complete.about_axis(center, angle, axis)


class RandomDoubleStretch(RandomManipulation):
    """A random bend by rotation a part of a molecule, works also on single rings"""
    num_hinge_atoms = 4

    def get_transformation(self, coordinates):
        """Construct a transformation object"""
        atom1, atom2, atom3, atom4 = self.hinge_atoms
        a = coordinates[atom1] - coordinates[atom2]
        a /= np.linalg.norm(a)
        b = coordinates[atom3] - coordinates[atom4]
        b /= np.linalg.norm(b)
        direction = 0.5*(a+b)
        direction *= np.random.uniform(-self.max_amplitude, self.max_amplitude)
        result = Translation(direction)
        return result


def iter_halfs_bond(graph):
    """Select a random bond (pair of atoms) that divides the molecule in two"""
    for atom1, atom2 in graph.edges:
        try:
            affected_atoms1, affected_atoms2 = graph.get_halfs(atom1, atom2)
            yield affected_atoms1, affected_atoms2, (atom1, atom2)
        except GraphError:
            # just try again
            continue


def iter_halfs_bend(graph):
    """Select randomly two consecutive bonds that divide the molecule in two"""
    for atom2 in range(graph.num_vertices):
        neighbors = list(graph.neighbors[atom2])
        for index1, atom1 in enumerate(neighbors):
            for atom3 in neighbors[index1+1:]:
                try:
                    affected_atoms = graph.get_halfs(atom2, atom1)[0]
                    # the affected atoms never contain atom1!
                    yield affected_atoms, (atom1, atom2, atom3)
                    continue
                except GraphError:
                    pass
                try:
                    affected_atoms = graph.get_halfs(atom2, atom3)[0]
                    # the affected atoms never contain atom3!
                    yield affected_atoms, (atom3, atom2, atom1)
                except GraphError:
                    pass


def iter_halfs_double(graph):
    """Select two random non-consecutive bonds that divide the molecule in two"""
    edges = graph.edges
    for index1, (atom_a1, atom_b1) in enumerate(edges):
        for atom_a2, atom_b2 in edges[:index1]:
            try:
                affected_atoms1, affected_atoms2, hinge_atoms = graph.get_halfs_double(atom_a1, atom_b1, atom_a2, atom_b2)
                yield affected_atoms1, affected_atoms2, hinge_atoms
            except GraphError:
                pass



def generate_manipulations(
    molecule, bond_stretch_factor=0.15, torsion_amplitude=np.pi,
    bending_amplitude=0.30
):
    """Generate a (complete) set of manipulations

       The result can be used as input for the functions 'randomize_molecule'
       and 'single_random_manipulation'

       Arguments:
         molecule  --  a reference geometry of the molecule, with graph
                       attribute
         bond_stretch_factor  --  the maximum relative change in bond length by
                                  one bond stretch manipulatio
         torsion_amplitude  --  the maximum change a dihdral angle
         bending_amplitude  --  the maximum change in a bending angle

       The return value is a list of RandomManipulation objects. They can be
       used to generate different random distortions of the original molecule.
    """
    do_stretch = (bond_stretch_factor > 0)
    do_double_stretch = (bond_stretch_factor > 0)
    do_bend = (bending_amplitude > 0)
    do_double_bend = (bending_amplitude > 0)
    do_torsion = (torsion_amplitude > 0)

    results = []
    # A) all manipulations that require one bond that cuts the molecule in half
    if do_stretch or do_torsion:
        for affected_atoms1, affected_atoms2, hinge_atoms in iter_halfs_bond(molecule.graph):
            if do_stretch:
                length = np.linalg.norm(
                    molecule.coordinates[hinge_atoms[0]] -
                    molecule.coordinates[hinge_atoms[1]]
                )
                results.append(RandomStretch(
                    affected_atoms1, length*bond_stretch_factor, hinge_atoms
                ))
            if do_torsion and len(affected_atoms1) > 1 and len(affected_atoms2) > 1:
                results.append(RandomTorsion(
                    affected_atoms1, torsion_amplitude, hinge_atoms
                ))
    # B) all manipulations that require a bending angle that cuts the molecule
    #    in two parts
    if do_bend:
        for affected_atoms, hinge_atoms in iter_halfs_bend(molecule.graph):
            results.append(RandomBend(
                affected_atoms, bending_amplitude, hinge_atoms
            ))
    # C) all manipulations that require two bonds that separate two halfs
    if do_double_stretch or do_double_bend:
        for affected_atoms1, affected_atoms2, hinge_atoms in iter_halfs_double(molecule.graph):
            if do_double_stretch:
                length1 = np.linalg.norm(
                    molecule.coordinates[hinge_atoms[0]] -
                    molecule.coordinates[hinge_atoms[1]]
                )
                length2 = np.linalg.norm(
                    molecule.coordinates[hinge_atoms[2]] -
                    molecule.coordinates[hinge_atoms[3]]
                )
                results.append(RandomDoubleStretch(
                    affected_atoms1, 0.5*(length1+length2)*bond_stretch_factor, hinge_atoms
                ))
            if do_double_bend and len(affected_atoms1) > 2 and len(affected_atoms2) > 2:
                if hinge_atoms[0] != hinge_atoms[2]:
                    results.append(RandomTorsion(
                        affected_atoms1, bending_amplitude, (hinge_atoms[0], hinge_atoms[2])
                    ))
                if hinge_atoms[1] != hinge_atoms[3]:
                    results.append(RandomTorsion(
                        affected_atoms2, bending_amplitude, (hinge_atoms[1], hinge_atoms[3])
                    ))
    # Neglect situations where three or more cuts are required.
    return results


def check_nonbond(molecule, thresholds):
    """Check whether all nonbonded atoms are well separated.

       If a nonbond atom pair is found that has an interatomic distance below
       the given thresholds. The thresholds dictionary has the following format:
       {frozenset([atom_number1, atom_number2]): distance}

       When random geometries are generated for sampling the conformational
       space of a molecule without strong repulsive nonbonding interactions, try
       to underestimate the thresholds at first instance and exclude bond
       stretches and bending motions for the random manipuulations. Then compute
       the forces projected on the nonbonding distance gradients. The distance
       for which the absolute value of these gradients drops below 100 kJ/mol is
       a coarse guess of a proper threshold value.
    """

    # check that no atoms overlap
    for atom1 in range(molecule.graph.num_vertices):
        for atom2 in range(atom1):
            if molecule.graph.distances[atom1, atom2] > 2:
                distance = np.linalg.norm(molecule.coordinates[atom1] - molecule.coordinates[atom2])
                if distance < thresholds[frozenset([molecule.numbers[atom1], molecule.numbers[atom2]])]:
                    return False
    return True


def randomize_molecule(molecule, manipulations, nonbond_thresholds, max_tries=1000):
    """Return a randomized copy of the molecule.

       If no randomized molecule can be generated that survives the nonbond
       check after max_tries repetitions, None is returned. In case of success,
       the randomized molecule is returned. The original molecule is not
       altered.
    """
    for m in range(max_tries):
        random_molecule = randomize_molecule_low(molecule, manipulations)
        if check_nonbond(random_molecule, nonbond_thresholds):
            return random_molecule


def randomize_molecule_low(molecule, manipulations):
    """Return a randomized copy of the molecule, without the nonbond check."""

    manipulations = copy.copy(manipulations)
    shuffle(manipulations)
    coordinates = molecule.coordinates.copy()
    for manipulation in manipulations:
        manipulation.apply(coordinates)
    return molecule.copy_with(coordinates=coordinates)


def single_random_manipulation(molecule, manipulations, nonbond_thresholds, max_tries=1000):
    """Apply a single random manipulation.

       If no randomized molecule can be generated that survives the nonbond
       check after max_tries repetitions, None is returned. In case of success,
       the randomized molecule and the corresponding transformation is returned.
       The original molecule is not altered.
    """
    for m in range(max_tries):
        random_molecule, transformation = single_random_manipulation_low(molecule, manipulations)
        if check_nonbond(random_molecule, nonbond_thresholds):
            return random_molecule, transformation
    return None


def single_random_manipulation_low(molecule, manipulations):
    """Return a randomized copy of the molecule, without the nonbond check."""

    manipulation = sample(manipulations, 1)[0]
    coordinates = molecule.coordinates.copy()
    transformation = manipulation.apply(coordinates)
    return molecule.copy_with(coordinates=coordinates), transformation


def random_dimer(molecule0, molecule1, thresholds, shoot_max):
    """Create a random dimer.

       molecule0 and molecule1 are placed in one reference frame at random
       relative positions. Interatomic distances are above the thresholds.
       Initially a dimer is created where one interatomic distance approximates
       the threshold value. Then the molecules are given an additional
       separation in the range [0, shoot_max].

       thresholds has the following format:
       {frozenset([atom_number1, atom_number2]): distance}
    """

    # apply a random rotation to molecule1
    center = np.zeros(3, float)
    angle = np.random.uniform(0, 2*np.pi)
    axis = random_unit()
    rotation = Complete.about_axis(center, angle, axis)
    cor1 = np.dot(molecule1.coordinates, rotation.r)

    # select a random atom in each molecule
    atom0 = np.random.randint(len(molecule0.numbers))
    atom1 = np.random.randint(len(molecule1.numbers))

    # define a translation of molecule1 that brings both atoms in overlap
    delta = molecule0.coordinates[atom0] - cor1[atom1]
    cor1 += delta

    # define a random direction
    direction = random_unit()
    cor1 += 1*direction

    # move molecule1 along this direction until all intermolecular atomic
    # distances are above the threshold values
    threshold_mat = np.zeros((len(molecule0.numbers), len(molecule1.numbers)), float)
    distance_mat = np.zeros((len(molecule0.numbers), len(molecule1.numbers)), float)
    for i1, n1 in enumerate(molecule0.numbers):
        for i2, n2 in enumerate(molecule1.numbers):
            threshold = thresholds.get(frozenset([n1, n2]))
            threshold_mat[i1, i2] = threshold**2
    while True:
        cor1 += 0.1*direction
        distance_mat[:] = 0
        for i in 0, 1, 2:
            distance_mat += np.subtract.outer(molecule0.coordinates[:, i], cor1[:, i])**2
        if (distance_mat > threshold_mat).all():
            break

    # translate over a random distance [0, shoot] along the same direction
    # (if necessary repeat until no overlap is found)
    while True:
        cor1 += direction*np.random.uniform(0, shoot_max)
        distance_mat[:] = 0
        for i in 0, 1, 2:
            distance_mat += np.subtract.outer(molecule0.coordinates[:, i], cor1[:, i])**2
        if (distance_mat > threshold_mat).all():
            break

    # done
    dimer = Molecule(
        np.concatenate([molecule0.numbers, molecule1.numbers]),
        np.concatenate([molecule0.coordinates, cor1])
    )
    dimer.direction = direction
    dimer.atom0 = atom0
    dimer.atom1 = atom1
    return dimer
