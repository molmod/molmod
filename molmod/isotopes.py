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
"""Databases with isotope information

   This module defines two database interfaces: ame2003 and nubtab03. The
   AME2003 [1]_, [2]_ database contains all the isotope masses and can be used as
   follows:

   >>> from molmod.isotopes import ame2003
   >>> print ame2003.masses[7][15]

   This would print the mass of nitrogen 15. The abundances are provided by the
   NUBTAB03 database [3]_. Use it as follows:

   >>> from molmod.isotopes import nubtab03
   >>> print nubtab.abundances[6][12]

   This would print the abundance of the carbon 12 isotope.

   References:

   .. [1] The AME2003 atomic mass evaluation (I). Evaluation of input data, adjustment
          procedures. A.H. Wapstra, G. Audi, and C. Thibault. Nuclear Physics A729,
          129 (2003).

   .. [2] The AME2003 atomic mass evaluation (II). Tables, graphs, and references. G.
          Audi, A.H. Wapstra, and C. Thibault. Nuclear Physics A729, 337 (2003).

   .. [3] The NUBASE evaluation of nuclear and decay properties. G. Audi, O. Bersillon,
          J. Blachot and A.H. Wapstra, Nuclear Physics A729, 3-128 (2003)
"""


from builtins import range
import pkg_resources

from molmod.units import amu


__all__ = ["Ame2003", "ame2003", "NubTab03", "nubtab03"]


class Ame2003(object):
    """An interface to a subset of the data from Ame2003.

       This object contains an attribute masses. This is a dictionary whose keys
       are the proton numbers (Z) and values are the corresponding values are
       again dictionaries. The latter dictionaries have the mass number (A) as
       keys and the corresponding isotope masses in atomic units as values. E.g.
       self.masses[6][12] is the mass of carbon 12.

       If you use this interface, cite the following references:

       The AME2003 atomic mass evaluation (I). Evaluation of input data, adjustment
       procedures. A.H. Wapstra, G. Audi, and C. Thibault. Nuclear Physics A729,
       129 (2003).

       The AME2003 atomic mass evaluation (II). Tables, graphs, and references. G.
       Audi, A.H. Wapstra, and C. Thibault. Nuclear Physics A729, 337 (2003).
    """
    def __init__(self):
        """
           An object of this type is created in this module, so there is not
           need to construct it manually.
        """
        self.masses = {}

        def add_mass(N, Z, mass):
            """Put a new mass into the dictionary"""
            n_masses = self.masses.setdefault(Z, {})
            n_masses[Z+N] = mass

        with pkg_resources.resource_stream(__name__, 'data/mass.mas03') as f:
            for i in range(39):
                next(f)
            for line in f:
                N = int(line[ 5:10])
                Z = int(line[10:15])
                mass = float(line[96:114].replace(b" ", b"").replace(b"#", b""))*1e-6*amu
                add_mass(N, Z, mass)


ame2003 = Ame2003()


class NubTab03(object):
    """An interface to a subset of the data that from NubTab03.

       This object contains an attribute abundances. This is a dictionary whose
       keys are the proton numbers (Z) and values are the corresponding values
       are again dictionaries. The latter dictionaries have the mass number (A)
       as keys and the corresponding isotope abundances as values. E.g.
       self.masses[6][12] is the abundance of carbon 12.

       If you use this interface, cite the following reference:

       The NUBASE evaluation of nuclear and decay properties. G. Audi, O. Bersillon,
       J. Blachot and A.H. Wapstra, Nuclear Physics A729, 3-128 (2003)
    """

    def __init__(self):
        """
           Argument:
            | ``filename``  --  the nubtab03 data file

           An object of this type is created in this module, so there is not
           need to construct it externally.
        """
        self.abundances = {}

        def add_abundance(A, Z, abundance):
            """Put a new abundance into the dictionary"""
            n_abundances = self.abundances.setdefault(Z, {})
            n_abundances[A] = abundance

        with pkg_resources.resource_stream(__name__, 'data/nubtab03.asc') as f:
            for line in f:
                Z = int(line[0:3])
                A = int(line[4:7])
                properties = dict(word.split(b"=") for word in line[106:].split(b";") if word.count(b"=")==1)
                abundance = properties.get(b'IS')
                if abundance is None: continue
                abundance = float(abundance.split()[0])
                add_abundance(A, Z, abundance)


nubtab03 = NubTab03()
