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


__all__ = [
    "FileError", "StateAttr", "ImmutableAttr", "ArrayAttr", "NumberState"
]


class FileError(Exception):
    pass


class StateAttr(object):
    def get(self, copy=False):
        raise NotImplementedError

    def get_kind(self, value):
        raise NotImplementedError

    def set(self, value):
        raise NotImplementedError


class ImmutableAttr(StateAttr):
    def __init__(self, owner, name):
        self.owner = owner
        self.name = name

    def get(self, copy=False):
        return getattr(self.owner, self.name)

    def get_kind(self, value):
        if isinstance(value, float):
            return 'f'
        elif isinstance(value, int) or isinstance(value, long):
            return 'i'
        else:
            raise ValueError("Only integer or floating point values can be stored.")

    def set(self, value):
        setattr(self.owner, self.name, value)

    def dump(self, f, name):
        # print the header line
        value = self.get()
        kind = self.get_kind(value)
        print >> f, "% 40s  kind=%s  value=%s" % (name, kind, value)


class ArrayAttr(StateAttr):
    def __init__(self, array):
        self.array = array
        if array.dtype.fields is not None:
            raise ValueError("Record arrays are not supported yet.")

    def get(self, copy=False):
        if copy:
            return self.array.copy()
        else:
            return self.array

    def get_kind(self, value):
        return value.dtype.kind

    def set(self, value):
        self.array[:] = value

    def dump(self, f, name):
        # print the header line
        print >> f, "% 40s  kind=%s  shape=%s" % (
            name,
            self.array.dtype.kind,
            ("%s" % (self.array.shape,)).replace(" ", ""),
        )
        # print the numbers
        counter = 0
        for value in self.array.flat:
            counter += 1
            print >> f, "% 20s" % value,
            if counter % 4 == 0:
                print >> f
        if counter % 4 != 0:
            print >> f

    def load(self, f, skip):
        counter = 0
        counter_limit = self.array.size
        convert = self.array.dtype.type
        while counter < counter_limit:
            line = f.readline()
            words = line.split()
            for word in words:
                if counter >= counter_limit:
                    raise FileError("Wrong array data: too many values.")
                if not skip:
                    self.array.flat[counter] = convert(word)
                counter += 1


class NumberState(object):
    def __init__(self):
        self._fields = {}

    def set_field(self, name, attr):
        if not isinstance(attr, StateAttr):
            raise TypeError("The second argument must a StateAttr instance.")
        if len(name) > 40:
            raise ValueError("Name can count at most 40 characters.")
        self._fields[name] = attr

    def get(self):
        return dict((name, attr.get(copy=True)) for name, attr in self._fields.iteritems())

    def set(self, new_fields):
        for name in new_fields:
            if name not in self._fields:
                raise ValueError("new_fields contains an unknown field.")
        if len(new_fields) != len(self._fields):
            raise ValueError("new_fields contains too many fields.")
        for name, attr in self._fields.iteritems():
            attr.set(new_fields[name])

    def dump(self, filename):
        f = file(filename, "w")
        for name in sorted(self._fields):
            self._fields[name].dump(f, name)
        f.close()

    def load(self, filename, subset=None):
        f = file(filename, "r")
        name = None
        num_names = 0

        while True:
            # read a header line
            line = f.readline()
            if len(line) == 0:
                break

            # process the header line
            words = line.split()
            name = words[0]
            attr = self._fields.get(name)
            if attr is None:
                raise FileError("Wrong header: unknown field %s" % name)

            if not words[1].startswith("kind="):
                raise FileError("Malformatted array header line. (kind)")
            kind = words[1][5:]
            if kind != attr.get_kind(attr.get()):
                raise FileError("Wrong header: kind of field %s does not match. Got %s, expected %s" % (name, kind, array.dtype.kind))

            skip = ((subset is not None) and (name not in subset))

            if (words[2].startswith("shape=(") and words[2].endswith(")")):
                if not isinstance(attr, ArrayAttr):
                    raise FileError("field '%s' is not an array." % name)
                shape = words[2][7:-1]
                if shape[-1]==',':
                    shape = shape[:-1]
                try:
                    shape = tuple(int(word) for word in shape.split(","))
                except ValueError:
                    raise FileError("Malformatted array header. (shape)")
                if shape != attr.get().shape:
                    raise FileError("Wrong header: shape of field %s does not match. Got %s, expected %s" % (name, shape, array.shape))
                attr.load(f, skip)
            elif words[2].startswith("value="):
                if not isinstance(attr, ImmutableAttr):
                    raise FileError("field '%s' is not a single value." % name)
                if not skip:
                    if kind == 'i':
                        attr.set(int(words[2][6:]))
                    else:
                        attr.set(float(words[2][6:]))
            else:
                raise FileError("Malformatted array header line. (shape/value)")

            num_names += 1

        if num_names != len(self._fields) and subset is None:
            raise FileError("Some fields are missing in the file.")
        f.close()


