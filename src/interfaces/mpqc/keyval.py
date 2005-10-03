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
    def __init__(self, items={}):
        self.items = items
    
    def write_stream(self, f, indent='', location='$:'):
        if location == '$:':
            self.clear_locations()
        for name, item in self.items.iteritems():
            if isinstance(item, int) or isinstance(item, float) or isinstance(item, str) or isinstance(item, list):
                print >> f, "%s%s=%s" % (indent, name, item)
            elif isinstance(item, KeyValObject):
                print >> f, "%s%s" % (indent, name), 
                item.write_stream(f, indent, location+':'+name)
            elif item!=None:
                raise KeyValError, "Object not supported %s=%s" % (name, item)

    def clear_locations(self):
        for item in self.items.itervalues():
            if isinstance(item, KeyValObject):
                item.clear_locations()


class KeyValObject(KeyVal):
    def __init__(self, class_name='', items={}):
        KeyVal.__init__(self, items)
        self.class_name = class_name
        self.location = None

    def write_stream(self, f, indent='', location='$:'):
        if location == None:
            if self.class_name != '':
                print >> f, "<%s>:(" % self.class_name
            else:
                print >> f, ":("
            KeyVal.write_stream(self, f, indent+'  ')
            print >> f, ")"
            self.location = location
        else:
            print >> f, "=%s" % self.location
            
    def clear_locations(self):
        self.location = None
        KeyVal.clear_locations(self)
