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
"""Tools for parsing CPMD trajectory files"""


from builtins import range
import numpy as np

from molmod.io.common import SlicedReader


__all__ = ["CPMDTrajectoryReader"]


class CPMDTrajectoryReader(SlicedReader):
    """A Reader for CPMD trajectory files

       Use this reader as an iterator::

         >>> ldr = CPMDTrajectoryReader("TRAJECTORY")
         >>> for pos, vel in ldr:
         ...     print pos[4,2]
    """
    def __init__(self, f, sub=slice(None)):
        """
           Arguments:
            | ``f``  --  a filename or a file-like object

           Optional argument:
            | ``sub``  --  a slice object indicating which time frames to skip/read
        """
        SlicedReader.__init__(self, f, sub)
        # first determine the number of atoms
        s = None
        self.num_atoms = 0
        for line in self._f:
            if s is None:
                s = line[:7]
            elif s != line[:7]:
                break
            self.num_atoms += 1
        self._f.seek(0) # go back to the beginning of the file

    def _read_frame(self):
        """Read and return the next time frame"""
        pos = np.zeros((self.num_atoms, 3), float)
        vel = np.zeros((self.num_atoms, 3), float)
        for i in range(self.num_atoms):
            line = next(self._f)
            words = line.split()
            pos[i, 0] = float(words[1])
            pos[i, 1] = float(words[2])
            pos[i, 2] = float(words[3])
            vel[i, 0] = float(words[4])
            vel[i, 1] = float(words[5])
            vel[i, 2] = float(words[6])
        return pos, vel

    def _skip_frame(self):
        """Skip the next time frame"""
        for i in range(self.num_atoms):
            line = next(self._f)
