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
"""Common functionality used by the molmod.io package"""

__all__ = ["slice_match", "FileFormatError"]


def slice_match(sub, counter):
    """Efficiently test if counter is in xrange(*sub)

       Arguments:
         sub  --  a slice object
         counter  -- an integer

       The function returns True if the counter is in
       xrange(sub.start, sub.stop, sub.step).
    """

    if sub.start is not None and counter < sub.start:
        return False
    if sub.stop is not None and counter >= sub.stop:
        raise StopIteration
    if sub.step is not None:
        if sub.start is None:
            if counter % sub.step != 0:
                return False
        else:
            if (counter - sub.start) % sub.step != 0:
                return False
    return True


class FileFormatError(Exception):
    """Is raised when unexpected data is encountered while reading a file."""
    pass


