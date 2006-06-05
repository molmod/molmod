# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of MolMod.
# 
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# 
# --


class KeyValError(Exception):
    pass


class KeyValWriter(object):
    def __init__(self, f, root, indent_step='  '):
        self.f = f
        self.root = root
        self.indent_step = indent_step
        self.locations = {} # map of objects to location strings
        
        # First decide where each object should be written and where references
        # should be used.
        for name, kvo in root.yield_attributes():
            assert isinstance(kvo, KeyValObject)
            self.locate_object(name, kvo, "$")
        
        # then write the output file.
        for name, kvo in root.yield_attributes():
            assert isinstance(kvo, KeyValObject)
            self.write_dispatch(name, kvo, "$", indent='')
    
    def locate_object(self, name, kvo, parent_location):
        def count_colons(s):
            return sum(1 for c in s if c==":")
            
        location = "%s:%s" % (parent_location, name)
        existing_location = self.locations.get(kvo)
        if existing_location == None or \
           count_colons(existing_location) > count_colons(location):
            self.locations[kvo] = location
        
        for child_name, child_kvo in kvo.yield_attributes():
            if isinstance(child_kvo, KeyValObject):
                self.locate_object(child_name, child_kvo, location)
    
    def write_object(self, name, kvo, parent_location, indent):
        assert kvo.class_name != '' or name != ''

        reference_location = self.locations.get(kvo)
        if parent_location != None and name != '':
            location = "%s:%s" % (parent_location, name)
        else:
            location = None

        if reference_location == None or reference_location == location:
            if kvo.class_name == '':
                self.f.write("%s%s:(\n" % (indent, name))
            else:
                self.f.write("%s%s<%s>:(\n" % (indent, name, kvo.class_name))
            for name, value in kvo.yield_attributes():
                self.write_dispatch(name, value, location, indent+self.indent_step)
            self.f.write("%s)\n" % indent)
        else:
            self.f.write("%s%s = %s\n" % (indent, name, reference_location))
            
    def write_dispatch(self, name, value, parent_location, indent):
        if isinstance(value, list):
            self.write_list(name, value, parent_location, indent)
        elif isinstance(value, KeyValObject):
            self.write_object(name, value, parent_location, indent)
        else:
            self.write_atom(name, value, parent_location, indent)
            
    def write_list(self, name, sequence, parent_location, indent):
        contains_kvo = reduce(lambda x,y: x or y, [isinstance(value, KeyValObject) for value in sequence], False)
        only_atoms = reduce(lambda x,y: x and y, [not isinstance(value, (KeyValObject, list)) for value in sequence], True)
        if name=='':
            if only_atoms:
              self.f.write("%s[ " % (indent))
            else:
              self.f.write("%s[\n" % (indent))
        else:
            if contains_kvo:
                symbol = ":"
            else:
                symbol = " ="
            if only_atoms:
                self.f.write("%s%s%s [ " % (indent, name, symbol))
            else:
                self.f.write("%s%s%s [\n" % (indent, name, symbol))
        for value in sequence:
            self.write_dispatch('', value, None, indent+self.indent_step)
        if only_atoms:
            self.f.write("]\n")
        else:
            self.f.write("%s]\n" % indent)

    def write_atom(self, name, value, parent_location, indent):
        if isinstance(value, basestring):
            formatted_value = '"%s"' % value
        else:
            formatted_value = str(value)
    
        if name=='':
            self.f.write("%s " % formatted_value)
        else:
            self.f.write("%s%s = %s\n" % (indent, name, formatted_value))


class KeyValObject(object):
    def __init__(self, class_name='', attributes={}):
        self.__dict__.update(attributes)
        self.class_name = class_name
        
    def yield_attributes(self):
        for name, value in self.__dict__.iteritems():
            if name != 'class_name':
                yield name, value

