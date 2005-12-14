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
        self.clear()
        self.add_parsers(file_parsers)
        
    def clear(self):
        self.file_parsers = {}
        
    def add_parsers(self, file_parsers):
        for file_parser in file_parsers:
            tag = (file_parser.filename, file_parser.extension)
            existing_file_parsers = self.file_parsers.get(tag)
            if existing_file_parsers == None:
                self.file_parsers[tag] = [file_parser]
            else:
                existing_file_parsers.append(file_parser)
            
    def parse(self, directory, prefix):
        result = {}
        for (filename, extension), file_parsers in self.file_parsers.iteritems():
            for file_parser in file_parsers:
                file_parser.reset()
            if extension:
                path = "%s/%s%s" % (directory, prefix, filename)
            else:
                path = "%s/%s" % (directory, filename)
            if isfile(path):
                f = file(path, 'r')
                for line in f:
                    #print line[:-1]
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


class MultiLineParser(FileParser):
    def __init__(self, label, activator, deactivator, condition=None):
        FileParser.__init__(self, label, condition)
        self.activator = activator
        self.deactivator = deactivator

    def reset(self):
        self.active = False

    def parse(self, line):
        if self.active:
            if self.deactivator != None and self.deactivator.search(line) != None:
                self.active = False
                self.stop_collecting()
            else:
                self.collect(line)
        elif self.activator != None and self.activator.search(line) != None:
            self.active = True
            self.start_collecting()
        elif self.activator == None and self.deactivator == None:
            self.collect(line)

    def start_collecting(self):
        raise NotImplementedError

    def collect(self, line):
        raise NotImplementedError

    def stop_collecting(self):
        raise NotImplementedError


