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



from molmod.molecules import Molecule

import re


__all__ = ["ReadError", "GDMA"]


class ReadError(Exception):
    pass


class GDMA(object):
    # written for gdma 2.2.01
    # http://www-stone.ch.cam.ac.uk/

    def __init__(self, filename):
        f = file(filename)
        self.load_from_file(f)
        f.close()

    def load_from_file(self, f):
        # 1) find the line "Using ___ density matrix from file ___"
        r = re.compile(r"^Using (?P<density>\S+) density matrix from file (?P<fchkfile>\S+)")
        self.density = None
        while True:
            line = f.readline()
            if len(line) == 0:
                raise ReadError("Could not find line \"Using ___ density matrix from file ___\"")
            match = r.search(line)
            if match is not None:
                self.density = match.group("density")
                break
        # 2) find the line "Multipole moments in atomic units ..."
        atomic_units = False
        while True:
            line = f.readline()
            if len(line) == 0:
                raise ReadError("Only atomic units are supported.")
            if line.startswith("Multipole moments in atomic units"):
                atomic_units = True
                break
        f.readline()

        # 3) read each multipole
        self.coordinates = []
        self.names = []
        self.ranks = []
        self.relative_radii = []
        self.multipole_norms = []
        self.multipoles = []

        line = f.readline()
        while not line.startswith("Total multipoles"):
            # read the name and the coordinates
            words = line.replace("=", "").split()
            self.names.append(words[0])
            self.coordinates.append([float(words[2]), float(words[4]), float(words[6])])
            # read the rank and radius
            words = f.readline().split()
            self.ranks.append(int(words[3]))
            self.relative_radii.append(float(words[7]))

            norm_lines, poles_lines, line = self.read_multipole(f)
            self.multipole_norms.append(norm_lines)
            self.multipoles.append(poles_lines)

        # read the origin
        words = f.readline().split()
        self.origin = [float(words[6][:-1]), float(words[9][:-1]), float(words[12])]
        # read the total multipole
        self.total_multipole_norms, self.total_multipole, foo = self.read_multipole(f)

    def read_multipole(self, f):
        # read the monopole line
        words = f.readline().split()
        last_norm = float(words[2])
        last_poles = [float(words[2])]
        # read the multipole lines
        norm_lines = []
        poles_lines = []
        while True:
            line = f.readline()
            if line == "\n": continue
            if len(line) == 0:
                raise ReadError("Unexpexted end of file.")
            words = line.replace("=", "").split()
            if not line.startswith("   "):
                norm_lines.append(last_norm)
                poles_lines.append(last_poles)
                if line.startswith("|Q"):
                    last_norm = float(words[1])
                    last_poles = [float(words[3]), float(words[5]), float(words[7])]
                else:
                    return norm_lines, poles_lines, line
            else:
                last_poles.extend(float(word) for index, word in enumerate(words) if index%2 == 1)







