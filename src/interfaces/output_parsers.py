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


from os.path import isfile

class OutputParser(object):
    def __init__(self, file_parsers=[]):
        self.file_parsers = {}
        self.add_parsers(file_parsers)
        
    def add_parsers(self, file_parsers):
        for file_parser in file_parsers:
            existing_file_parsers = self.file_parsers.get(file_parser.extension)
            if existing_file_parsers == None:
                self.file_parsers[file_parser.extension] = [file_parser]
            else:
                existing_file_parsers.append(file_parser)
            
    def parse(self, prefix):
        result = {}
        for extension, file_parsers in self.file_parsers.iteritems():
            for file_parser in file_parsers:
                file_parser.reset()
            filename = "%s.%s" % (prefix, extension)
            if isfile(filename):
                f = file(filename, 'r')
                for line in f:
                    for file_parser in file_parsers:
                        file_parser.conditioned_parse(line)
                f.close()
            for file_parser in file_parsers:
                result[file_parser.label] = file_parser.result()
        return result


class FileParser(object):
    extension = None

    def __init__(self, label, condition=None):
        self.label = label
        self.condition = condition
        self.reset()
    
    def reset(self):
        raise NotImplementedError
    
    def conditioned_parse(self, line):
        if (self.condition == None) or self.condition():
            self.parse(line)
    
    def parse(self, line):
        raise NotImplementedError

    def result(self):
        raise NotImplementedError
