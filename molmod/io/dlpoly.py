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
"""Readers for DLPoly file formats"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.units import picosecond, amu, angstrom, atm, deg
from molmod.io.common import SlicedReader, FileFormatError


__all__ = ["DLPolyHistoryReader", "DLPolyOutputReader"]


class DLPolyHistoryReader(SlicedReader):
    """A Reader for the DLPoly history file format.

       Use this object as an iterator::

         >>> hr = HistoryReader("somefile.txt")
         >>> for frame in hr:
         ...     print frame["cell"]

    """
    def __init__(self, f, sub=slice(None), pos_unit=angstrom,
        vel_unit=angstrom/picosecond, frc_unit=amu*angstrom/picosecond**2,
        time_unit=picosecond, mass_unit=amu, restart=False,
    ):
        """
           Arguments:
             | ``f``  --  a filename or a file-like object

           Optional arguments:
             | ``sub``  --  a slice indicating the frames to be skipped/selected
             | ``pos_unit``, ``vel_unit``, ``frc_unit``, ``time_unit``,
               ``mass_unit``  --  The conversion factors for the unit conversion
                from the units in the data file to atomic units. The defaults of
                these optional arguments correspond to the defaults of dlpoly.

           When the file starts with a line that satisfies the following
           conditions, it is assumed that this is a history restart file:

           * line consists of 6 words
           * first word equals 'timestep'
           * the following for words are integers
           * the last word is a float
        """
        SlicedReader.__init__(self, f, sub)
        self._counter = 1 # make our counter compatible with dlpoly
        self.pos_unit = pos_unit
        self.vel_unit = vel_unit
        self.frc_unit = frc_unit
        self.time_unit = time_unit
        self.mass_unit = mass_unit
        restart = self._detect_restart()
        if restart is None:
            try:
                self.header = next(self._f)[:-1]
                integers = tuple(int(word) for word in next(self._f).split())
                if len(integers) != 3:
                    raise FileFormatError("Second line must contain three integers.")
                self.keytrj, self.imcon, self.num_atoms = integers
            except StopIteration:
                raise FileFormatError("File is too short. Could not read header.")
            except ValueError:
                raise FileFormatError("Second line must contain three integers.")
        else:
            self.header = ''
            self.num_atoms, self.keytrj, self.imcon = restart
        self._frame_size = 4 + self.num_atoms*(self.keytrj+2)

    def _detect_restart(self):
        words = next(self._f).split()
        self._f.seek(0)
        if len(words) != 6:
            return
        if words[0] != 'timestep':
            return
        for i in 1, 2, 3, 4:
            if not words[i].isdigit():
                return
        try:
            float(words[5])
        except ValueError:
            return
        return int(words[2]), int(words[3]), int(words[4])


    def _read_frame(self):
        """Read a single frame from the trajectory"""
        # auxiliary read function
        def read_three(msg):
            """Read three words as floating point numbers"""
            line = next(self._f)
            try:
                return [float(line[:12]), float(line[12:24]), float(line[24:])]
            except ValueError:
                raise FileFormatError(msg)

        frame = {}
        # read the frame header line
        words = next(self._f).split()
        if len(words) != 6:
            raise FileFormatError("The first line of each time frame must contain 6 words. (%i'th frame)" % self._counter)
        if words[0] != "timestep":
            raise FileFormatError("The first word of the first line of each time frame must be 'timestep'. (%i'th frame)" % self._counter)
        try:
            step = int(words[1])
            frame["step"] = step
            if int(words[2]) != self.num_atoms:
                raise FileFormatError("The number of atoms has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            if int(words[3]) != self.keytrj:
                raise FileFormatError("keytrj has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            if int(words[4]) != self.imcon:
                raise FileFormatError("imcon has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            frame["timestep"] = float(words[5])*self.time_unit
            frame["time"] = frame["timestep"]*step # this is ugly, or wait ... dlpoly is a bit ugly. we are not to blame!
        except ValueError:
            raise FileFormatError("Could not convert all numbers on the first line of the current time frame. (%i'th frame)" % self._counter)
        # the three cell lines
        cell = np.zeros((3, 3), float)
        frame["cell"] = cell
        cell_msg = "The cell lines must consist of three floating point values. (%i'th frame, %i'th step)" % (self._counter, step)
        for i in range(3):
            cell[:, i] = read_three(cell_msg)
        cell *= self.pos_unit
        # the atoms
        symbols = []
        frame["symbols"] = symbols
        masses = np.zeros(self.num_atoms, float)
        frame["masses"] = masses
        charges = np.zeros(self.num_atoms, float)
        frame["charges"] = charges
        pos = np.zeros((self.num_atoms, 3), float)
        frame["pos"] = pos
        if self.keytrj > 0:
            vel = np.zeros((self.num_atoms, 3), float)
            frame["vel"] = vel
        if self.keytrj > 1:
            frc = np.zeros((self.num_atoms, 3), float)
            frame["frc"] = frc
        for i in range(self.num_atoms):
            # the atom header line
            words = next(self._f).split()
            if len(words) != 4:
                raise FileFormatError("The atom header line must contain 4 words. (%i'th frame, %i'th step, %i'th atom)" % (self._counter, step, i+1))
            symbols.append(words[0])
            try:
                masses[i] = float(words[2])*self.mass_unit
                charges[i] = float(words[3])
            except ValueError:
                raise FileFormatError("The numbers in the atom header line could not be interpreted.")
            # the pos line
            pos_msg = "The position lines must consist of three floating point values. (%i'th frame, %i'th step, %i'th atom)" % (self._counter, step, i+1)
            pos[i] = read_three(pos_msg)
            if self.keytrj > 0:
                vel_msg = "The velocity lines must consist of three floating point values. (%i'th frame, %i'th step, %i'th atom)" % (self._counter, step, i+1)
                vel[i] = read_three(vel_msg)
            if self.keytrj > 1:
                frc_msg = "The force lines must consist of three floating point values. (%i'th frame, %i'th step, %i'th atom)" % (self._counter, step, i+1)
                frc[i] = read_three(frc_msg)
        pos *= self.pos_unit # convert to au
        if self.keytrj > 0:
            vel *= self.vel_unit # convert to au
        if self.keytrj > 1:
            frc *= self.frc_unit # convert to au
        return frame

    def _skip_frame(self):
        """Skip a single frame from the trajectory"""
        for i in range(self._frame_size):
            next(self._f)


class DLPolyOutputReader(SlicedReader):
    """A Reader for DLPoly output files.

       Use this object as an iterator::

         >>> outr = OutputReader("somefile.txt")
         >>> for row in outr:
         ...     print row[5]

       The variable row in the example above is a concatenation of all the
       values that belong to one time frame. (line after line)
    """

    _marker = " " + "-"*130

    def __init__(self, f, sub=slice(None), skip_equi_period=True,
        pos_unit=angstrom, time_unit=picosecond, angle_unit=deg,
        e_unit=amu/(angstrom/picosecond)**2
    ):
        """
           Arguments:
            | ``f``  --  a filename or a file-like object

           Optional arguments:
            | ``sub``  --  a slice indicating the frames to be skipped/selected
            | ``skip_equi_period``  -- When True, the equilibration period is
                                       not read [default=True]
            | ``pos_unit``, ``time_unit``, ``angle_unit``, ``e_unit`` --  The
               conversion factors for the unit conversion from the units in the
               data file to atomic units. The defaults of these optional
               arguments correspond to the defaults of dlpoly.

        """
        SlicedReader.__init__(self, f, sub)
        self._counter = 1 # make our counter compatible with dlpoly
        self.skip_equi_period = skip_equi_period

        self._conv = [
            1,         e_unit,      1, e_unit, e_unit, e_unit,     e_unit,     e_unit,     e_unit, e_unit,
            time_unit, e_unit,      1, e_unit, e_unit, e_unit,     e_unit,     e_unit,     e_unit, e_unit,
            1,         pos_unit**3, 1, e_unit, e_unit, angle_unit, angle_unit, angle_unit, e_unit, 1000*atm,
        ]

        # find the line that gives the number of equilibration steps:
        try:
            while True:
                line = next(self._f)
                if line.startswith(" equilibration period"):
                    self.equi_period = int(line[30:])
                    break
        except StopIteration:
            raise FileFormatError("DL_POLY OUTPUT file is too short. Could not find line with the number of equilibration steps.")
        except ValueError:
            raise FileFormatError("Could not read the number of equilibration steps. (expecting an integer)")

    def goto_next_frame(self):
        """Continue reading until the next frame is reached"""
        marked = False
        while True:
            line = next(self._f)[:-1]
            if marked and len(line) > 0 and not line.startswith(" --------"):
                try:
                    step = int(line[:10])
                    return step, line
                except ValueError:
                    pass
            marked = (len(line) == 131 and line == self._marker)

    def _read_frame(self):
        """Read a single frame from the trajectory"""
        # optionally skip the equilibration
        if self.skip_equi_period:
            while True:
                step, line = self.goto_next_frame()
                self._counter += 1
                if step >= self.equi_period:
                    break
            self.skip_equi_period = False
        else:
            step, line = self.goto_next_frame()

        # read the three lines
        try:
            row = [step]
            for i in range(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
            line = next(self._f)[:-1]
            row.append(float(line[:10]))
            for i in range(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
            line = next(self._f)[:-1]
            row.append(float(line[:10]))
            for i in range(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
        except ValueError:
            raise FileFormatError("Some numbers in the output file could not be read. (expecting floating point numbers)")

        # convert all the numbers to atomic units
        for i in range(30):
            row[i] *= self._conv[i]

        # done
        return row

    def _skip_frame(self):
        """Skip a single frame from the trajectory"""
        self.goto_next_frame()
