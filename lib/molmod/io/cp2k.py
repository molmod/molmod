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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


from molmod.units import angstrom

import numpy


__all__ = ["CellReader", "FileFormatError", "Section", "Keyword", "InputFile"]


class CellReader(object):
    def __init__(self, filename):
        self.f = file(filename)

    def __del__(self):
        self.f.close()

    def __iter__(self):
        return self

    def next(self):
        # skip the blablabla
        self.f.next()
        data = []
        data.append([float(word) for word in self.f.next().split()])
        data.append([float(word) for word in self.f.next().split()])
        data.append([float(word) for word in self.f.next().split()])
        return numpy.array(data)*angstrom


class FileFormatError(Exception):
    pass


class Section(object):
    def __init__(self, name="", children=None, section_parameters=""):
        if not isinstance(name, str):
            raise TypeError("A name must be a string, got %s." % name)
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
        if len(self.__order) != sum(len(values) for values in self.__index.itervalues()):
            return False
        import copy
        tmp = copy.copy(self.__order)
        for key, values in self.__index.iteritems():
            for value in values:
                if value.name != key:
                    return False
                if value in tmp:
                    tmp.remove(value)
                else:
                    return False
                if isinstance(value, Section):
                    if not value._consistent():
                        return False
        return True

    def __iter__(self):
        return iter(self.order)

    def __eq__(self, other):
        if not (isinstance(other, Section) and self.name == other.name):
            return False
        if len(self.__index) != len(other.__index):
            #print "len(self.__index) != len(other.__index)"
            return False
        for key, lself in self.__index.iteritems():
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
                    if not (isinstance(item, Section) or isinstance(item, Keyword)):
                        raise TypeError("The value must be an Section or a Keyword, got: %s." % value)
                    if item.name != key:
                        raise KeyError("The name of the child must correspond to the given key, %s!=%s" % (key, item.name))
                self.__index[key] = value
                # update the ordered list
                index = delete_all(key)
                for item in value[::-1]:
                    self.__order.insert(index, item)
            elif isinstance(value, Section) or isinstance(value, Keyword):
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
        return self.__name

    name = property(get_name)

    def append(self, child):
        if not (isinstance(child, Section) or isinstance(child, Keyword)):
            raise TypeError("The child must be an Section or a Keyword, got: %s." % value)
        l = self.__index.setdefault(child.name, [])
        l.append(child)
        self.__order.append(child)

    def insert(self, index, child):
        if not (isinstance(child, Section) or isinstance(child, Keyword)):
            raise TypeError("The child must be an Section or a Keyword, got: %s." % value)
        l = self.__index.setdefault(child.name, [])
        l.append(child)
        self.__order.insert(index, child)

    def dump_children(self, f, indent=''):
        for child in self.__order:
            child.dump(f, indent+'  ')

    def dump(self, f, indent=''):
        print >> f, ("%s&%s %s" % (indent, self.__name, self.section_parameters)).rstrip()
        self.dump_children(f, indent)
        print >> f, "%s&END %s" % (indent, self.__name)

    def readline(self, f):
        while True:
            line = f.readline()
            if len(line) == 0:
                raise EOFError
            line = line[:line.find('#')]
            line = line.strip()
            if len(line) > 0:
                return line

    def load_children(self, f):
        while True:
            line = self.readline(f)
            if line[0] == '&':
                if line[1:].startswith("END"):
                    check_name = line[4:].strip().upper()
                    if check_name != self.__name:
                        raise FileFormatError("Section end mismatch, pos=%s", f.tell())
                    break
                else:
                    section = Section()
                    section.load(f, line)
                    self.append(section)
            else:
                keyword = Keyword()
                keyword.load(line)
                self.append(keyword)


    def load(self, f, line):
        words = line[1:].split()
        self.__name = words[0].upper()
        self.section_parameters = " ".join(words[1:])
        try:
            self.load_children(f)
        except EOFError:
            raise FileFormatError("Unexpected end of file, section '%s' not ended." % self.__name)


class Keyword(object):
    def __init__(self, name="", value=""):
        self.__name = name.upper()
        self.__value = value

    def __eq__(self, other):
        #print (self.name, other.name), (self.value, other.value)
        return (
            isinstance(other, Keyword) and
            self.name == other.name and
            self.value == other.value
        )

    def dump(self, f, indent=''):
        print >> f, ("%s%s %s" % (indent, self.__name, self.__value)).rstrip()

    def load(self, line):
        words = line.split()
        try:
            float(words[0])
            self.__name == ""
            self.__value = " ".join(words)
        except ValueError:
            self.__name = words[0].upper()
            self.__value = " ".join(words[1:])


    def set_value(self, value):
        if not isinstance(value, str):
            raise TypeError("A value must be a string, got %s." % value)
        self.__value = value

    def get_name(self):
        return self.__name

    def get_value(self):
        return self.__value

    value = property(get_value, set_value)
    name = property(get_name)


class InputFile(Section):
    @staticmethod
    def read_from_file(filename):
        f = file(filename)
        result = InputFile()
        try:
            while True:
                result.load_children(f)
        except EOFError:
            pass
        f.close()
        return result

    def __init__(self, children=None):
        Section.__init__(self, "__ROOT__", children)

    def write_to_file(self, filename):
        f = file(filename, "w")
        self.dump(f)
        f.close()

    def dump(self, f, indent=''):
        self.dump_children(f, indent)



