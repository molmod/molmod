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
# --


from molmod import context
from molmod.units import u

import os


__all__ = ["ame2003", "nubtab03"]


class Ame2003(object):
    """An interface to a subset of the data from Ame2003.

    If you use this interface, also refer to:

    The AME2003 atomic mass evaluation (I). Evaluation of input data, adjustment
    procedures. A.H. Wapstra, G. Audi, and C. Thibault. Nuclear Physics A729,
    129 (2003).

    The AME2003 atomic mass evaluation (II). Tables, graphs, and references. G.
    Audi, A.H. Wapstra, and C. Thibault. Nuclear Physics A729, 337 (2003).
    """
    def __init__(self, filename):
        self.masses = {}

        def add_mass(N, Z, mass):
            n_masses = self.masses.setdefault(Z, {})
            n_masses[Z+N] = mass

        f = file(filename)
        for i in xrange(39):
            f.next()

        for line in f:
            N = int(line[ 5:10])
            Z = int(line[10:15])
            mass = float(line[96:114].replace(" ", "").replace("#", ""))*1e-6*u
            add_mass(N, Z, mass)

        f.close()


ame2003 = Ame2003(os.path.join(context.share_path, "mass.mas03"))


class NubTab03(object):
    """An interface to a subset of the data that from NubTab03.

    If you use this interface, also refer to:

    The NUBASE evaluation of nuclear and decay properties. G. Audi, O. Bersillon,
    J. Blachot and A.H. Wapstra, Nuclear Physics A729, 3-128 (2003)
    """

    def __init__(self, filename):
        self.abundances = {}

        def add_abundance(A, Z, abundance):
            n_abundances = self.abundances.setdefault(Z, {})
            n_abundances[A] = abundance

        f = file(filename)
        for line in f:
            Z = int(line[0:3])
            A = int(line[4:7])
            properties = dict(word.split("=") for word in line[106:].split(";") if word.count("=")==1)
            abundance = properties.get('IS')
            if abundance is None: continue
            abundance = float(abundance.split()[0])
            add_abundance(A, Z, abundance)

        f.close()


nubtab03 = NubTab03(os.path.join(context.share_path, "nubtab03.asc"))




