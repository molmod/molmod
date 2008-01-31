# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import copy, numpy


__all__ = ["FileError", "ArrayState"]


class FileError(Exception):
    pass


class ArrayState(object):
    def __init__(self):
        self._fields = {}

    def set_field(self, name, array):
        if not isinstance(array, numpy.ndarray):
            raise TypeError("The second argument must be a numpy array.")
        if array.dtype.fields is not None:
            raise ValueError("Record arrays are not supported yet.")
        if len(name) > 40:
            raise ValueError("Name can count at most 40 characters.")
        self._fields[name] = array

    def get(self):
        return copy.deepcopy(self._fields)

    def set(self, new_fields):
        for name in new_fields:
            if name not in self._fields:
                raise ValueError("new_fields contains an unknown field.")
        if len(new_fields) != len(self._fields):
            raise ValueError("new_fields contains too many fields.")
        for name, array in self._fields.iteritems():
            array[:] = new_fields[name]

    def dump(self, filename):
        f = file(filename, "w")
        for name in sorted(self._fields):
            array = self._fields[name]
            # print the header line
            print >> f, "% 40s  kind=%s  shape=%s" % (
                name,
                array.dtype.kind,
                ("%s" % (array.shape,)).replace(" ", ""),
            )
            # print the numbers
            counter = 0
            for value in array.flat:
                counter += 1
                print >> f, "% 20s" % value,
                if counter % 4 == 0:
                    print >> f
            if counter % 4 != 0:
                print >> f
        f.close()

    def load(self, filename):
        f = file(filename, "r")
        array = None
        name = None
        counter = None
        counter_limit = None
        num_names = 0

        for line in f:
            if name is None:
                # read a header line
                words = line.split()
                name = words[0]
                array = self._fields.get(name)
                if array is None:
                    raise FileError("Wrong header: unknown field %s" % name)

                if not words[1].startswith("kind="):
                    raise FileError("Malformatted array header line. (kind)")
                kind = words[1][5:]
                if kind != array.dtype.kind:
                    raise FileError("Wrong header: kind of field %s does not match. Got %s, expected %s" % (name, kind, array.dtype.kind))

                if not (words[2].startswith("shape=(") and words[2].endswith(")")):
                    raise FileError("Malformatted array header line. (shape)")
                shape = words[2][7:-1]
                if shape[-1]==',':
                    shape = shape[:-1]
                try:
                    shape = tuple(int(word) for word in shape.split(","))
                except ValueError:
                    raise FileError("Malformatted array header. (shape)")
                if shape != array.shape:
                    raise FileError("Wrong header: shape of field %s does not match. Got %s, expected %s" % (name, shape, array.shape))

                counter = 0
                counter_limit = array.size

                self._fields[name]
            else:
                words = line.split()
                for word in words:
                    if counter >= counter_limit:
                        raise FileError("Wrong array data: too many values.")
                    array.flat[counter] = float(word)
                    counter += 1
                if counter == counter_limit:
                    # time for the next array
                    num_names += 1
                    name = None
                    array = None
                    counter = None
                    counter_limit = None
        if num_names != len(self._fields):
            raise FileError("Some fields are missing in the file.")
        f.close()
