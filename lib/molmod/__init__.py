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


import sys, os


__all__ = ["context"]


class Error(Exception):
    pass


class Context(object):
    def __init__(self):
        share_dirs = set([
            os.path.join(sys.prefix, "share/molmod"),
            os.path.join(sys.prefix, "local/share/molmod"),
            os.path.join(str(os.getenv("HOME")), "local/share/molmod"),
            "/usr/share/molmod",
            "/usr/local/share/molmod",
        ])
        self._share_dirs = []
        for share_dir in share_dirs:
            if os.path.isdir(share_dir):
                self._share_dirs.append(share_dir)
        if len(self._share_dirs) == 0:
            raise Error("Could not find shared files.")

    def get_share_filename(self, filename):
        for share_dir in self._share_dirs:
            result = os.path.join(share_dir, filename)
            if os.path.isfile(result):
                return result
        raise ValueError("No file '%s' found in the share directories." % filename)


context = Context()




