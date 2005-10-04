# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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


class KeyVal(object):
    def __init__(self, named_items=[]):
        self.names = [pair[0] for pair in named_items]
        self.dict_items = dict(named_items)
        
    def __getitem__(self, name):
        return self.dict_items[name]
        
    def __setitem__(self, name, item):
        if name not in self.dict_items:
            self.names.append(name)
        self.dict_items[name] = item
        
    def write_list(self, f, l, indent='', noindent=False):
        if not isinstance(l, list):
            f.write("%s\n" % l)
        elif reduce(lambda x, y: x or y, (isinstance(item, list) for item in l), False):
            if noindent:
                f.write("[\n")
            else:
                f.write(indent+"[\n")
            for item in l:
                self.write_list(f, item, indent + '  ')
            f.write(indent+"]\n")
        else:
            if noindent:
                f.write('[%s' % l[0])
            else:
                f.write(indent+'[%s' % l[0])
            for item in l[1:]:
                f.write(" %s" % item)
            f.write("]\n")
    
    def write_stream(self, f, indent='', location='$'):
        if location == '$':
            self.clear_locations()
        for name in self.names:
            item = self.dict_items[name]
            if isinstance(item, int) or isinstance(item, float) or isinstance(item, str) or isinstance(item, list):
                f.write("%s%s=" % (indent, name))
                self.write_list(f, item, indent, True)
            elif isinstance(item, KeyValObject):
                f.write("%s%s" % (indent, name)) 
                item.write_stream(f, indent, location+':'+name)
            elif item!=None:
                raise KeyValError, "Object not supported %s=%s" % (name, item)

    def clear_locations(self):
        for item in self.dict_items.itervalues():
            if isinstance(item, KeyValObject):
                item.clear_locations()


class KeyValObject(KeyVal):
    def __init__(self, class_name='', named_items=[]):
        KeyVal.__init__(self, named_items)
        self.class_name = class_name
        self.location = None

    def write_stream(self, f, indent='', location='$'):
        if self.location == None:
            if self.class_name != '':
                f.write("<%s>:(\n" % self.class_name)
            else:
                f.write(":(\n")
            KeyVal.write_stream(self, f, indent+'  ', location)
            f.write(indent+")\n")
            self.location = location
        else:
            f.write("=%s\n" % self.location)
            
    def clear_locations(self):
        self.location = None
        KeyVal.clear_locations(self)
