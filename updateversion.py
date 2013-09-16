#!/usr/bin/env python
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
#--
#!/usr/bin/env python

import re, sys



rules = [
    ('setup.py', '^    version=\'(...+)\',$'),
    ('molmod/__init__.py', '^__version__ = \'(...+)\'$'),
    ('doc/conf.py', '^version = \'(...+)\'$'),
    ('doc/conf.py', '^release = \'(...+)\'$'),
    ('doc/tutorial/install.rst', '^    http://users.ugent.be/~tovrstra/molmod/molmod-(...+).tar.gz.$'),
    ('doc/tutorial/install.rst', '^    wget http://users.ugent.be/~tovrstra/molmod/molmod-(...+).tar.gz$'),
    ('doc/tutorial/install.rst', '^    tar -xvzf molmod-(...+).tar.gz$'),
    ('doc/tutorial/install.rst', '^    cd molmod-(...+)$'),
]


if __name__ == '__main__':
    newversion = sys.argv[1]

    for fn, regex in rules:
        r = re.compile(regex)
        with open(fn) as f:
            lines = f.readlines()
        for i in xrange(len(lines)):
            line = lines[i]
            m = r.match(line)
            if m is not None:
                lines[i] = line[:m.start(1)] + newversion + line[m.end(1):]
        with open(fn, 'w') as f:
            f.writelines(lines)
