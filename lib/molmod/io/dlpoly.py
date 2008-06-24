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


from molmod.units import ps, amu, A, atm, deg
from molmod.io.common import slice_match

import numpy


__all__ = ["Error", "HistoryReader"]


class Error(Exception):
    pass


class HistoryReader(object):
    def __init__(self, filename, sub=slice(None), pos_unit=A, vel_unit=A/ps, frc_unit=amu*A/ps**2, time_unit=ps, mass_unit=amu):
        self._f = file(filename)
        self._sub = sub
        self.pos_unit = pos_unit
        self.vel_unit = vel_unit
        self.frc_unit = frc_unit
        self.time_unit = time_unit
        self.mass_unit = mass_unit
        try:
            self.header = self._f.next()[:-1]
            integers = tuple(int(word) for word in self._f.next().split())
            if len(integers) != 3:
                raise Error("Second line must contain three integers.")
            self.keytrj, self.imcon, self.num_atoms = integers
        except StopIteration:
            raise Error("File is too short. Could not read header.")
        except ValueError:
            raise Error("Second line must contain three integers.")
        self._counter = 1
        self._frame_size = 4 + self.num_atoms*(self.keytrj+2)

    def __del__(self):
        self._f.close()

    def __iter__(self):
        return self

    def next(self):
        # auxiliary read function
        def read_three(msg):
            # read three words as floating point numbers
            line = self._f.next()
            try:
                return [float(line[:12]), float(line[12:24]), float(line[24:])]
            except ValueError:
                raise Error(msg)

        # skip frames as requested
        while not slice_match(self._sub, self._counter):
            for i in xrange(self._frame_size):
                self._f.next()
            self._counter += 1

        frame = {}
        # read the frame header line
        words = self._f.next().split()
        if len(words) != 6:
            raise Error("The first line of each time frame must contain 6 words. (%i'th frame)" % self._counter)
        if words[0] != "timestep":
            raise Error("The first word of the first line of each time frame must be 'timestep'. (%i'th frame)" % self._counter)
        try:
            step = int(words[1])
            frame["step"] = step
            if int(words[2]) != self.num_atoms:
                raise Error("The number of atoms has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            if int(words[3]) != self.keytrj:
                raise Error("keytrj has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            if int(words[4]) != self.imcon:
                raise Error("imcon has changed. (%i'th frame, %i'th step)" % (self._counter, step))
            frame["timestep"] = float(words[5])*self.time_unit
            frame["time"] = frame["timestep"]*step # this is ugly, or wait ... dlpoly is a bit ugly. we are not to blame!
        except ValueError:
            raise Error("Could not convert all numbers on the first line of the current time frame. (%i'th frame)" % self._counter)
        # the three cell lines
        cell = numpy.zeros((3,3), float)
        frame["cell"] = cell
        cell_msg = "The cell lines must consist of three floating point values. (%i'th frame, %i'th step)" % (self._counter, step)
        for i in xrange(3):
            cell[:,i] = read_three(cell_msg)
        cell *= self.pos_unit
        # the atoms
        symbols = []
        frame["symbols"] = symbols
        masses = numpy.zeros(self.num_atoms, float)
        frame["masses"] = masses
        charges = numpy.zeros(self.num_atoms, float)
        frame["charges"] = charges
        pos = numpy.zeros((self.num_atoms,3), float)
        frame["pos"] = pos
        if self.keytrj > 0:
            vel = numpy.zeros((self.num_atoms,3), float)
            frame["vel"] = vel
        if self.keytrj > 1:
            frc = numpy.zeros((self.num_atoms,3), float)
            frame["frc"] = frc
        for i in xrange(self.num_atoms):
            # the atom header line
            words = self._f.next().split()
            if len(words) != 4:
                raise Error("The atom header line must contain 4 words. (%i'th frame, %i'th step, %i'th atom)" % (self._counter, step, i+1))
            symbols.append(words[0])
            try:
                masses[i] = float(words[2])*self.mass_unit
                charges[i] = float(words[3])
            except ValueError:
                raise Error("The numbers in the atom header line could not be interpreted.")
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
        # done
        self._counter += 1
        return frame


class OutputReader(object):
    _marker = " " + "-"*130

    def __init__(self, filename, sub=slice(None), skip_equi_period=True, pos_unit=A, time_unit=ps, angle_unit=deg, e_unit=amu/(A/ps)**2):
        self._f = file(filename)
        self._sub = sub
        self.skip_equi_period = skip_equi_period
        self._counter = 1

        self._conv = [
            1,         e_unit,      1, e_unit, e_unit, e_unit,     e_unit,     e_unit,     e_unit, e_unit,
            time_unit, e_unit,      1, e_unit, e_unit, e_unit,     e_unit,     e_unit,     e_unit, e_unit,
            1,         pos_unit**3, 1, e_unit, e_unit, angle_unit, angle_unit, angle_unit, e_unit, 1000*atm,
        ]
        self.last_step = None

        # find the line that gives the number of equilibration steps:
        try:
            while True:
                line = self._f.next()
                if line.startswith(" equilibration period"):
                    self.equi_period = int(line[30:])
                    break
        except StopIteration:
            raise Error("DL_POLY OUTPUT file is too short. Could not find line with the number of equilibration steps.")
        except ValueError:
            raise Error("Could not read the number of equilibration steps. (expecting an integer)")

    def __del__(self):
        self._f.close()

    def __iter__(self):
        return self

    def next(self):
        def goto_next_frame():
            marked = False
            while True:
                line = self._f.next()[:-1]
                if marked and len(line) > 0 and not line.startswith(" --------"):
                    try:
                        step = int(line[:10])
                        return step, line
                    except ValueError:
                        pass
                marked = (len(line) == 131 and line == self._marker)

        while True:
            step, line = goto_next_frame()
            if (not self.skip_equi_period or step >= self.equi_period) and \
               step != self.last_step:
                break

        # skip frames as requested
        while not slice_match(self._sub, self._counter):
            step, line = goto_next_frame()
            self._counter += 1

        # now really read these three lines
        try:
            row = [step]
            for i in xrange(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
            line = self._f.next()[:-1]
            row.append(float(line[:10]))
            for i in xrange(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
            line = self._f.next()[:-1]
            row.append(float(line[:10]))
            for i in xrange(9):
                row.append(float(line[10+i*12:10+(i+1)*12]))
        except ValueError:
            raise Error("Some numbers in the output file could not be read. (expecting floating point numbers)")

        # convert all the numbers to atomic units
        for i in xrange(30):
            row[i] *= self._conv[i]

        # done
        self.last_step = step
        return row


