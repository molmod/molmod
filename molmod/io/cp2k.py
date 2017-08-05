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
"""Tools for generating CP2K input files and a Reader for unit cell trajectories"""


from __future__ import print_function

from molmod.io.common import FileFormatError


__all__ = ["CP2KSection", "CP2KKeyword", "CP2KInputFile"]


class CP2KSection(object):
    """Data structure representing a section in a CP2K input file

       Instances of this class are iterable and support full element access.
    """

    def __init__(self, name="", children=None, section_parameters=""):
        """
           Optional arguments:
            | ``name``  --  the name of the section
            | ``children``  --  a list of CP2KSection and CP2KKeyword objects
            | ``section_parameters``  --  the value assigned to the section
        """
        if not isinstance(name, str):
            raise TypeError("A name must be a string, got %s." % name)
        if children is not None:
            for child in children:
                if not (isinstance(child, CP2KSection) or isinstance(child, CP2KKeyword)):
                    raise TypeError("All children must be CP2KSection or CP2KKeyword objects, found: %s." % child)
        if not isinstance(section_parameters, str):
            raise TypeError("The section_parameters argument must be a string, got %s." % section_parameters)
        self.__name = name.upper()
        self.__index = {}
        self.__order = []
        if children is not None:
            for child in children:
                self.append(child)
        self.section_parameters = section_parameters

    def _consistent(self):
        """Checks the constency between self.__index and self.__order"""
        if len(self.__order) != sum(len(values) for values in self.__index.values()):
            return False
        import copy
        tmp = copy.copy(self.__order)
        for key, values in self.__index.items():
            for value in values:
                if value.name != key:
                    return False
                if value in tmp:
                    tmp.remove(value)
                else:
                    return False
                if isinstance(value, CP2KSection):
                    if not value._consistent():
                        return False
        return True

    def __iter__(self):
        return iter(self.__order)

    def __eq__(self, other):
        if not (isinstance(other, CP2KSection) and self.name == other.name):
            return False
        if len(self.__index) != len(other.__index):
            #print "len(self.__index) != len(other.__index)"
            return False
        for key, lself in self.__index.items():
            lother = other.__index.get(key)
            if lother is None:
                #print "lother==None"
                return False
            if len(lother) != len(lself):
                #print len(lother), "==", len(lself)
                return False
            for iself, iother in zip(lself, lother):
                if not iself == iother: return False
        return True

    def __len__(self):
        return len(self.__order)

    def __getitem__(self, key):
        if isinstance(key, str):
            l = self.__index[key]
            if len(l) == 1:
                return l[0]
            else:
                return l
        elif isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], str) and isinstance(key[1], int):
            return self.__index[key[0]][key[1]]
        else:
            raise TypeError("Unsupported key: %s" % key)

    def __setitem__(self, key, value):
        def delete_all(key):
            """Delete all children whose name is equal to the argument key"""
            indexes = []
            for index, item in enumerate(self.__order):
                if item.name == key:
                    indexes.append(index)
            for index in indexes:
                del self.__order[index]
            if len(indexes) > 0:
                return indexes[0]
            else:
                return -1

        if isinstance(key, str):
            if isinstance(value, list):
                for item in value:
                    if not (isinstance(item, CP2KSection) or isinstance(item, CP2KKeyword)):
                        raise TypeError("The value must be a CP2KSection or a CP2KKeyword, got: %s." % value)
                    if item.name != key:
                        raise KeyError("The name of the child must correspond to the given key, %s!=%s" % (key, item.name))
                self.__index[key] = value
                # update the ordered list
                index = delete_all(key)
                for item in value[::-1]:
                    self.__order.insert(index, item)
            elif isinstance(value, CP2KSection) or isinstance(value, CP2KKeyword):
                if value.name != key:
                    raise KeyError("The name of the child must correspond to the given key, %s!=%s" % (key, value.name))
                self.__index[key] = [value]
                index = delete_all(key)
                self.__order.insert(index, value)
        elif isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], str) and isinstance(key[1], int):
            if value.name != key[0]:
                raise KeyError("The name of the child must correspond to the given key, %s!=%s" % (key, value.name))
            index = self.__order.index(self.__index[key[0]][key[1]])
            self.__index[key[0]][key[1]] = value
            self.__order[index] = value
        else:
            raise TypeError("Unsupported key: %s" % key)

    def __delitem__(self, key):
        if isinstance(key, str):
            items = self.__index[key]
            for item in items:
                self.__order.remove(item)
            del self.__index[key]
        elif isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], str) and isinstance(key[1], int):
            self.__order.remove(self.__index[key[0]][key[1]])
            del self.__index[key[0]][key[1]]
        else:
            raise TypeError("Unsupported key: %s" % key)

    def get_name(self):
        """Return the name of this section"""
        return self.__name

    name = property(get_name)

    def append(self, child):
        """Add a child section or keyword"""
        if not (isinstance(child, CP2KSection) or isinstance(child, CP2KKeyword)):
            raise TypeError("The child must be a CP2KSection or a CP2KKeyword, got: %s." % child)
        l = self.__index.setdefault(child.name, [])
        l.append(child)
        self.__order.append(child)

    def insert(self, index, child):
        """Insert a child section or keyword"""
        if not (isinstance(child, CP2KSection) or isinstance(child, CP2KKeyword)):
            raise TypeError("The child must be a CP2KSection or a CP2KKeyword, got: %s." % child)
        l = self.__index.setdefault(child.name, [])
        l.append(child)
        self.__order.insert(index, child)

    def dump_children(self, f, indent=''):
        """Dump the children of the current section to a file-like object"""
        for child in self.__order:
            child.dump(f, indent+'  ')

    def dump(self, f, indent=''):
        """Dump this section and its children to a file-like object"""
        print(("%s&%s %s" % (indent, self.__name, self.section_parameters)).rstrip(), file=f)
        self.dump_children(f, indent)
        print("%s&END %s" % (indent, self.__name), file=f)

    def readline(self, f):
        """A helper method that only reads uncommented lines"""
        while True:
            line = f.readline()
            if len(line) == 0:
                raise EOFError
            line = line[:line.find('#')]
            line = line.strip()
            if len(line) > 0:
                return line

    def load_children(self, f):
        """Load the children of this section from a file-like object"""
        while True:
            line = self.readline(f)
            if line[0] == '&':
                if line[1:].startswith("END"):
                    check_name = line[4:].strip().upper()
                    if check_name != self.__name:
                        raise FileFormatError("CP2KSection end mismatch, pos=%s", f.tell())
                    break
                else:
                    section = CP2KSection()
                    section.load(f, line)
                    self.append(section)
            else:
                keyword = CP2KKeyword()
                keyword.load(line)
                self.append(keyword)


    def load(self, f, line=None):
        """Load this section from a file-like object"""
        if line is None:
            # in case the file contains only a fragment of an input file,
            # this is useful.
            line = f.readlin()
        words = line[1:].split()
        self.__name = words[0].upper()
        self.section_parameters = " ".join(words[1:])
        try:
            self.load_children(f)
        except EOFError:
            raise FileFormatError("Unexpected end of file, section '%s' not ended." % self.__name)


class CP2KKeyword(object):
    """Data structure representing a keyword-value pair in a CP2K input file"""

    def __init__(self, name="", value="", unit=None):
        """
           Optional arguments:
            | ``name``  --  The keyword name
            | ``value``  --  The associated value
            | ``unit``  --  The unit (user must guarantee validity of the unit)
        """
        self.__name = name.upper()
        self.__value = value
        self.__unit = unit

    def __eq__(self, other):
        #print (self.name, other.name), (self.value, other.value)
        return (
            isinstance(other, CP2KKeyword) and
            self.name == other.name and
            self.value == other.value and
            self.unit == other.unit
        )

    def dump(self, f, indent=''):
        """Dump this keyword to a file-like object"""
        if self.__unit is None:
            print(("%s%s %s" % (indent, self.__name, self.__value)).rstrip(), file=f)
        else:
            print(("%s%s [%s] %s" % (indent, self.__name, self.__unit, self.__value)).rstrip(), file=f)

    def load(self, line):
        """Load this keyword from a file-like object"""
        words = line.split()
        try:
            float(words[0])
            self.__name = ""
            self.__value = " ".join(words)
        except ValueError:
            self.__name = words[0].upper()
            if len(words) > 2 and words[1][0]=="[" and words[1][-1]=="]":
                self.unit = words[1][1:-1]
                self.__value = " ".join(words[2:])
            else:
                self.__value = " ".join(words[1:])


    def set_value(self, value):
        """Set the value associated with the keyword"""
        if not isinstance(value, str):
            raise TypeError("A value must be a string, got %s." % value)
        self.__value = value

    def get_name(self):
        """Get the name of the keyword"""
        return self.__name

    def get_value(self):
        """Get the value of the keyword"""
        return self.__value

    def get_unit(self):
        """Get the name of the unit"""
        return self.__unit

    value = property(get_value, set_value)
    name = property(get_name)
    unit = property(get_unit)


class CP2KInputFile(CP2KSection):
    """Data structure representing an entire CP2K input file"""

    @staticmethod
    def read_from_file(filename):
        """
           Arguments:
            | ``filename``  --  the filename of the input file

           Use as follows::

             >>> if = CP2KInputFile.read_from_file("somefile.inp")
             >>> for section in if:
             ...     print section.name
        """
        with open(filename) as f:
            result = CP2KInputFile()
            try:
                while True:
                    result.load_children(f)
            except EOFError:
                pass
        return result

    def __init__(self, children=None):
        """Initialize a new (empty) CP2KInputFile object"""
        CP2KSection.__init__(self, "__ROOT__", children)

    def write_to_file(self, filename):
        """Write the CP2KInput data structure to a file"""
        with open(filename, "w") as f:
            self.dump(f)

    def dump(self, f, indent=''):
        """Dump the CP2KInput data structure to a file-like object"""
        self.dump_children(f, indent)
