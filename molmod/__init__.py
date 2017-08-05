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
"""Collection of utilities to automate molecular modeling research

   The scope of MolMod is twofold. (A) On the one hand does it provide the
   necessary building blocks for user-made batch scripts that facilitate
   molecular modeling research. Often it is just a matter of taking away the
   burden of manual repetitive work when processing computational data or when
   preparing simulations. By providing a package with various common building
   blocks, one also avoids the repetitive work in the implementation of
   user-made scripts. (B) The second purpose of MolMod is to serve as a toolkit
   for more serious and long-term molecular modeling projects such as Zeobuilder
   or MD-Tracks. Such programs also benefit from the general-purpose nature of
   MolMod. This imposes a certain level of modularity and hence facilitates
   maintainability.

   Most of the submodules are loaded directly into the molmod namespace. There
   are some exceptions. Due to long loading times the submodules with large
   amounts of data must be imported explicitly: molmod.bonds molmod.isotopes,
   molmod.periodic. Example:

   from molmod.periodic import periodic

   One must also load the molmod.io subpackage explicitely. It is not loaded
   by default to avoid namespace flooding.
"""


from molmod.version import __version__


import numpy as np
np.seterr(divide='raise', invalid='raise')

from molmod.binning import *
from molmod.clusters import *
from molmod.constants import *
from molmod.graphs import *
from molmod.ic import *
from molmod.log import *
from molmod.minimizer import *
from molmod.molecules import *
from molmod.molecular_graphs import *
from molmod.pairff import *
from molmod.quaternions import *
from molmod.randomize import *
from molmod.similarity import *
from molmod.symmetry import *
from molmod.toyff import *
from molmod.transformations import *
from molmod.unit_cells import *
from molmod.units import *
from molmod.utils import *
from molmod.vectors import *
from molmod.zmatrix import *
