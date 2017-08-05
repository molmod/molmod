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
"""Tools for reading Gamess punch files"""


from __future__ import division

from builtins import range
import numpy as np

from molmod.units import angstrom, amu


__all__ = ["PunchFile"]


class PunchFile(object):
    """Reader for GAMESS Punch files.

       After initialization, the data from the file is available as attributes
       of the created object. The following attributes are created if the
       corresponding data is present in the punch file: ``title``, ``symmetry``,
       ``symbols``, ``numbers``, ``coordinates``, ``energy``, ``gradient`` and
       ``hessian``.

    """
    def __init__(self, filename):
        """
           Argument:
            | ``filename`` --  The file to load from.
        """
        self.filename = filename
        self._read(filename)

    def _read(self, filename):
        """Internal routine that reads all data from the punch file."""
        data = {}
        parsers = [
            FirstDataParser(), CoordinateParser(), EnergyGradParser(),
            SkipApproxHessian(), HessianParser(), MassParser(),
        ]
        with open(filename) as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                # at each line, a parsers checks if it has to process a piece of
                # file. If that happens, the parser gets control over the file
                # and reads as many lines as it needs to collect data for some
                # attributes.
                for parser in parsers:
                    if parser.test(line, data):
                        parser.read(line, f, data)
                        break
        self.__dict__.update(data)


class PunchParser(object):
    """Base class for all parsers.

       A parser has two essential methods: test and read. The test routine
       checks the current line and previously collected data to see if the
       parser should be activated. If that happens, the read routine will be
       called by PunchFile._read to perform the actual extraction of data.
    """

    def test(self, line, data):
        """Returns True of this parser can process some of the following lines

           Arguments:
            | ``line``  --  The current line.
            | ``data``  --  The data dictionary constructed so far.
        """
        raise NotImplementedError

    def read(self, line, f, data):
        """Reads part of a file and adds the extracted information to data.

           Arguments:
            | ``line``  --  The current line.
            | ``f``  --  A file object to read the following lines.
            | ``data``  --  The data dictionary constructed so far.
        """
        raise NotImplementedError


class FirstDataParser(PunchParser):
    """Extracts ``title``, ``symmetry`` and ``symbols`` from the punch file."""

    def __init__(self):
        self.used = False

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return not self.used and line == " $DATA\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        self.used = True
        data["title"] = f.readline().strip()
        data["symmetry"] = f.readline().split()[0]
        if data["symmetry"] != "C1":
            raise NotImplementedError("Only C1 symmetry is supported.")
        symbols = []
        while line != " $END      \n":
            line = f.readline()
            if line[0] != " ":
                symbols.append(line.split()[0])
        data["symbols"] = symbols



class CoordinateParser(PunchParser):
    """Extracts ``numbers`` and ``coordinates`` from the punch file."""

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return "symbols" in data and line == " COORDINATES OF SYMMETRY UNIQUE ATOMS (ANGS)\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        f.readline()
        f.readline()
        N = len(data["symbols"])
        # if the data are already read before, just overwrite them
        numbers = data.get("numbers")
        if numbers is None:
            numbers = np.zeros(N, int)
            data["numbers"] = numbers
        coordinates = data.get("coordinates")
        if coordinates is None:
            coordinates = np.zeros((N,3), float)
            data["coordinates"] = coordinates
        for i in range(N):
            words = f.readline().split()
            numbers[i] = int(float(words[1]))
            coordinates[i,0] = float(words[2])*angstrom
            coordinates[i,1] = float(words[3])*angstrom
            coordinates[i,2] = float(words[4])*angstrom


class EnergyGradParser(PunchParser):
    """Extracts ``energy`` and ``gradient`` from the punch file."""

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return "symbols" in data and line == " $GRAD\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        data["energy"] = float(f.readline().split()[1])
        N = len(data["symbols"])
        # if the data are already read before, just overwrite them
        gradient = data.get("gradient")
        if gradient is None:
            gradient = np.zeros((N,3), float)
            data["gradient"] = gradient
        for i in range(N):
            words = f.readline().split()
            gradient[i,0] = float(words[2])
            gradient[i,1] = float(words[3])
            gradient[i,2] = float(words[4])


class SkipApproxHessian(PunchParser):
    """Prevents the approximate Hessian from being parsed."""

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return line == "CAUTION, APPROXIMATE HESSIAN!\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        line = f.readline()
        assert(line == " $HESS\n")
        while line != " $END\n":
            line = f.readline()


class HessianParser(PunchParser):
    """Extracts ``hessian`` from the punch file."""

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return "symbols" in data and line == " $HESS\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        assert("hessian" not in data)
        f.readline()
        N = len(data["symbols"])
        hessian = np.zeros((3*N, 3*N), float)
        tmp = hessian.ravel()
        counter = 0
        while True:
            line = f.readline()
            if line == " $END\n":
                break
            line = line[5:-1]
            for j in range(len(line)//15):
                tmp[counter] = float(line[j*15:(j+1)*15])
                counter += 1
        data["hessian"] = hessian


class MassParser(PunchParser):
    """Extracts ``masses`` from the punch file."""

    def test(self, line, data):
        """See :meth:`PunchParser.test`"""
        return "symbols" in data and line == "ATOMIC MASSES\n"

    def read(self, line, f, data):
        """See :meth:`PunchParser.read`"""
        N = len(data["symbols"])
        masses = np.zeros(N, float)
        counter = 0
        while counter < N:
            words = f.readline().split()
            for word in words:
                masses[counter] = float(word)*amu
                counter += 1
        data["masses"] = masses
