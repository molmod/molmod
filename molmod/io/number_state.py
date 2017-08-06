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
"""Persistance, i.e. storage on disk, for objects with numerical attributes"""


from __future__ import print_function

import numpy as np

from molmod.io.common import FileFormatError


__all__ = ["NumberState"]


class StateAttr(object):
    """Base class for NumberState attributes"""

    def __init__(self, owner, name):
        """
           Arguments:
            | ``owner``  --  the instance to read the attribute from
            | ``name``  --  the name of the attribute
        """
        self.owner = owner
        self.name = name

    def get(self, copy=False):
        """Return the value of the attribute"""
        raise NotImplementedError

    def get_kind(self, value):
        """Return the kind (type) of the attribute"""
        raise NotImplementedError

    def set(self, value):
        """Set the value of the attribute"""
        raise NotImplementedError


class ScalarAttr(StateAttr):
    """A scalar attribute for NumberState objects"""

    def get(self, copy=False):
        """Return the value of the attribute"""
        return getattr(self.owner, self.name)

    def get_kind(self, value):
        """Return the kind (type) of the attribute"""
        if isinstance(value, float):
            return 'f'
        elif isinstance(value, int):
            return 'i'
        else:
            raise ValueError("Only integer or floating point values can be stored.")

    def set(self, value):
        """Set the value of the attribute"""
        setattr(self.owner, self.name, value)

    def dump(self, f, name):
        """Write the attribute to a file-like object"""
        # print the header line
        value = self.get()
        kind = self.get_kind(value)
        print("% 40s  kind=%s  value=%s" % (name, kind, value), file=f)


class ArrayAttr(StateAttr):
    """An array attribute for the NumberState object"""
    def __init__(self, owner, name):
        """Initialize a ArrayAttr object

           Arguments:
             ``owner``  --  the instance to read the attribute from
             ``name``  --  the name of the attribute
        """
        StateAttr.__init__(self, owner, name)
        array = self.get()
        if array.dtype.fields is not None:
            raise ValueError("Record arrays are not supported yet.")

    def get(self, copy=False):
        """Return the value of the attribute"""
        array = getattr(self.owner, self.name)
        if copy:
            return array.copy()
        else:
            return array

    def get_kind(self, value):
        """Return the kind (type) of the attribute"""
        return value.dtype.kind

    def set(self, value):
        """Set the value of the attribute"""
        array = self.get()
        array[:] = value

    def dump(self, f, name):
        """Write the attribute to a file-like object"""
        array = self.get()
        # print the header line
        print("% 40s  kind=%s  shape=(%s)" % (
            name,
            array.dtype.kind,
            ",".join([str(int(size_axis)) for size_axis in array.shape]),
        ), file=f)
        # print the numbers
        counter = 0
        for value in array.flat:
            counter += 1
            print("% 20s" % value, end=' ', file=f)
            if counter % 4 == 0:
                print(file=f)
        if counter % 4 != 0:
            print(file=f)

    def load(self, f, skip):
        """Load the array data from a file-like object"""
        array = self.get()
        counter = 0
        counter_limit = array.size
        convert = array.dtype.type
        while counter < counter_limit:
            line = f.readline()
            words = line.split()
            for word in words:
                if counter >= counter_limit:
                    raise FileFormatError("Wrong array data: too many values.")
                if not skip:
                    array.flat[counter] = convert(word)
                counter += 1


class NumberState(object):
    """Component class for data structures with human-readable persistence.

       The format used to save and load the object is similar to a formatted
       checkpoint file from the Gaussian package. Some additional info is
       stored such as the shape of an array and the exact data type of the
       array elements.

       The attributes that contain data to be read from or to be written to
       files are set up in the constructor of the owner class. This is a
       typical simple example::

         >>> class Foo(object):
         ...     def __init__(self, a, b):
         ...         self.a = a
         ...         self.b = b
         ...         self.state = NumberState(self, ["a", "b"])

       In this example a is an array and b is a single scalar. One can now
       read/write these attributes to a file as follows:

         >>> foo = Foo(a, b)
         >>> foo.state.dump("somefile.txt")
         >>> foo.state.load("somefile.txt")
    """
    def __init__(self, owner, names):
        """
           Arguments:
            | ``owner``  --  the object whose attributes are dumped and loaded
            | ``names``  --  a list of attribute names to dump and load
        """
        self._owner = owner
        self._fields = {}
        for name in names:
            value = getattr(owner, name)
            if isinstance(value, np.ndarray):
                self._register(name, ArrayAttr)
            elif isinstance(value, int) or isinstance(value, float):
                self._register(name, ScalarAttr)
            else:
                raise TypeError("Can not handle attribute %s=%s" % (name, value))

    def _register(self, name, AttrCls):
        """Register a new attribute to take care of with dump and load

           Arguments:
            | ``name``  --  the name to be used in the dump file
            | ``AttrCls``  --  an attr class describing the attribute
        """
        if not issubclass(AttrCls, StateAttr):
            raise TypeError("The second argument must a StateAttr instance.")
        if len(name) > 40:
            raise ValueError("Name can count at most 40 characters.")
        self._fields[name] = AttrCls(self._owner, name)

    def get(self, subset=None):
        """Return a dictionary object with the registered fields and their values

           Optional rgument:
            | ``subset``  --  a list of names to restrict the number of fields
                              in the result
        """
        if subset is None:
            return dict((name, attr.get(copy=True)) for name, attr in self._fields.items())
        else:
            return dict((name, attr.get(copy=True)) for name, attr in self._fields.items() if name in subset)

    def set(self, new_fields, subset=None):
        """Assign the registered fields based on a dictionary

           Argument:
            | ``new_fields``  --  the dictionary with the data to be assigned to
                                  the attributes

           Optional argument:
            | ``subset``  --  a list of names to restrict the fields that are
                              effectively overwritten
        """
        for name in new_fields:
            if name not in self._fields and (subset is None or name in subset):
                raise ValueError("new_fields contains an unknown field '%s'." % name)
        if subset is not None:
            for name in subset:
                if name not in self._fields:
                    raise ValueError("name '%s' in subset is not a known field in self._fields." % name)
                if name not in new_fields:
                    raise ValueError("name '%s' in subset is not a known field in new_fields." % name)
        if subset is None:
            if len(new_fields) != len(self._fields):
                raise ValueError("new_fields contains too many fields.")
        for name, attr in self._fields.items():
            if name in subset:
                attr.set(new_fields[name])

    def dump(self, filename):
        """Dump the registered fields to a file

           Argument:
            | ``filename``  --  the file to write to
        """
        with open(filename, "w") as f:
            for name in sorted(self._fields):
                self._fields[name].dump(f, name)

    def load(self, filename, subset=None):
        """Load data into the registered fields

           Argument:
            | ``filename``  --  the filename to read from

           Optional argument:
            | ``subset``  --  a list of field names that are read from the file.
                              If not given, all data is read from the file.
        """
        with open(filename, "r") as f:
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
                    raise FileFormatError("Wrong header: unknown field %s" % name)

                if not words[1].startswith("kind="):
                    raise FileFormatError("Malformatted array header line. (kind)")
                kind = words[1][5:]
                expected_kind = attr.get_kind(attr.get())
                if kind != expected_kind:
                    raise FileFormatError("Wrong header: kind of field %s does not match. Got %s, expected %s" % (name, kind, expected_kind))

                skip = ((subset is not None) and (name not in subset))

                print(words)
                if (words[2].startswith("shape=(") and words[2].endswith(")")):
                    if not isinstance(attr, ArrayAttr):
                        raise FileFormatError("field '%s' is not an array." % name)
                    shape = words[2][7:-1]
                    if shape[-1] == ', ':
                        shape = shape[:-1]
                    try:
                        shape = tuple(int(word) for word in shape.split(","))
                    except ValueError:
                        raise FileFormatError("Malformatted array header. (shape)")
                    expected_shape = attr.get().shape
                    if shape != expected_shape:
                        raise FileFormatError("Wrong header: shape of field %s does not match. Got %s, expected %s" % (name, shape, expected_shape))
                    attr.load(f, skip)
                elif words[2].startswith("value="):
                    if not isinstance(attr, ScalarAttr):
                        raise FileFormatError("field '%s' is not a single value." % name)
                    if not skip:
                        if kind == 'i':
                            attr.set(int(words[2][6:]))
                        else:
                            attr.set(float(words[2][6:]))
                else:
                    raise FileFormatError("Malformatted array header line. (shape/value)")

                num_names += 1

            if num_names != len(self._fields) and subset is None:
                raise FileFormatError("Some fields are missing in the file.")
