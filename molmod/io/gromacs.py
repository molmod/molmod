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
"""Readers for Gromacs related file formats"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.units import picosecond, nanometer
from molmod.io.common import SlicedReader


__all__ = ["GroReader"]


class GroReader(SlicedReader):
    """A reader from .gro trajectory files from gromacs

       Use this reader as an iterator::

         >>> gr = GroReader("somefile.gro")
         >>> for time, pos, vel, cell in gr:
         ...     print pos
    """
    def __init__(self, f, sub=slice(None)):
        """
           Argument:
            | ``f``  --  a filename or a file-like object

           Optional argument:
            | ``sub``  --  a slice object to indicate the frames to be read/skipped.
        """
        SlicedReader.__init__(self, f, sub)
        self.num_atoms = None
        pos = self._read_frame()[1]
        self.num_atoms = len(pos)
        self._f.seek(0)

    def _get_line(self):
        """Get a line or raise StopIteration"""
        line = self._f.readline()
        if len(line) == 0:
            raise StopIteration
        return line

    def _read_frame(self):
        """Read one frame"""
        # Read the first line, ignore the title and try to get the time. The
        # time field is optional.
        line = self._get_line()
        pos = line.rfind("t=")
        if pos >= 0:
            time = float(line[pos+2:])*picosecond
        else:
            time = 0.0
        # Read the second line, the number of atoms must match with the first
        # frame.
        num_atoms = int(self._get_line())
        if self.num_atoms is not None and self.num_atoms != num_atoms:
            raise ValueError("The number of atoms must be the same over the entire file.")
        # Read the atom lines
        pos = np.zeros((num_atoms, 3), np.float32)
        vel = np.zeros((num_atoms, 3), np.float32)
        for i in range(num_atoms):
            words = self._get_line()[22:].split()
            pos[i, 0] = float(words[0])
            pos[i, 1] = float(words[1])
            pos[i, 2] = float(words[2])
            vel[i, 0] = float(words[3])
            vel[i, 1] = float(words[4])
            vel[i, 2] = float(words[5])
        pos *= nanometer
        vel *= nanometer/picosecond
        # Read the cell line
        cell = np.zeros((3, 3), np.float32)
        words = self._get_line().split()
        if len(words) >= 3:
            cell[0, 0] = float(words[0])
            cell[1, 1] = float(words[1])
            cell[2, 2] = float(words[2])
        if len(words) == 9:
            cell[1, 0] = float(words[3])
            cell[2, 0] = float(words[4])
            cell[0, 1] = float(words[5])
            cell[2, 1] = float(words[6])
            cell[0, 2] = float(words[7])
            cell[1, 2] = float(words[8])
        cell *= nanometer
        return time, pos, vel, cell

    def _skip_frame(self):
        """Skip one frame"""
        self._get_line()
        num_atoms = int(self._get_line())
        if self.num_atoms is not None and self.num_atoms != num_atoms:
            raise ValueError("The number of atoms must be the same over the entire file.")
        for i in range(num_atoms+1):
            self._get_line()
