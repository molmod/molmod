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

import sys, os


__all__ = ["context"]


class Context(object):
    """Global variable to find and use share directory"""
    def __init__(self):
        """Initialize the Context object

           This is done once when importing a molmod module. There is no need
           to do this externally a second time.
        """
        # find the data files
        fn_datadir = os.path.join(os.path.dirname(__file__), "datadir.txt")
        if os.path.isfile(fn_datadir):
            f = file(fn_datadir)
            datadir = f.readline().strip()
            f.close()
            self.share_dir = os.path.join(datadir, "share", "molmod")
        else:
            # When running from the build directory for the tests.
            self.share_dir = "../share"
        if not os.path.isdir(self.share_dir):
            raise RuntimeError("Share dir '%s' does not exist." % self.share_dir)

    def get_share_filename(self, filename):
        """Retrieve the full path for a given filename in the share folder"""
        result = os.path.join(self.share_dir, filename)
        if not os.path.isfile(result):
            raise ValueError("Data file '%s' not found." % result)
        return result


context = Context()


from molmod.binning import *
from molmod.clusters import *
from molmod.graphs import *
from molmod.ic import *
from molmod.minimizer import *
from molmod.molecules import *
from molmod.molecular_graphs import *
from molmod.quaternions import *
from molmod.randomize import *
from molmod.similarity import *
from molmod.toyff import *
from molmod.transformations import *
from molmod.unit_cells import *
from molmod.units import *
from molmod.vectors import *
from molmod.volume import *
from molmod.zmatrix import *


