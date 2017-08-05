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


from __future__ import print_function

import numpy as np


__all__ = ['load_chk', 'dump_chk']


def load_chk(filename):
    '''Load a checkpoint file

       Argument:
        | filename  --  the file to load from

       The return value is a dictionary whose keys are field labels and the
       values can be None, string, integer, float, boolean or an array of
       strings, integers, booleans or floats.

       The file format is similar to the Gaussian fchk format, but has the extra
       feature that the shapes of the arrays are also stored.
    '''
    with open(filename) as f:
        result = {}
        while True:
            line = f.readline()
            if line == '':
                break
            if len(line) < 54:
                raise IOError('Header lines must be at least 54 characters long.')
            key = line[:40].strip()
            kind = line[47:52].strip()
            value = line[53:-1] # discard newline
            if kind == 'str':
                result[key] = value
            elif kind == 'int':
                result[key] = int(value)
            elif kind == 'bln':
                result[key] = value == 'True'
            elif kind == 'flt':
                result[key] = float(value)
            elif kind[3:5] == 'ar':
                if kind[:3] == 'str':
                    dtype = np.dtype('U')
                elif kind[:3] == 'int':
                    dtype = int
                elif kind[:3] == 'bln':
                    dtype = bool
                elif kind[:3] == 'flt':
                    dtype = float
                else:
                    raise IOError('Unsupported kind: %s' % kind)
                shape = tuple(int(i) for i in value.split(','))
                array = np.zeros(shape, dtype)
                if array.size > 0:
                    work = array.ravel()
                    counter = 0
                    while True:
                        short = f.readline().split()
                        if len(short) == 0:
                            raise IOError('Insufficient data')
                        for s in short:
                            if dtype is bool and s.lower() in ['True', '1', 'Yes']:
                                work[counter] = True
                            elif dtype == np.dtype('U'):
                                work[counter] = s
                            else:
                                work[counter] = dtype(s)
                            counter += 1
                            if counter == array.size:
                                break
                        if counter == array.size:
                            break
                result[key] = array
            elif kind == 'none':
                result[key] = None
            else:
                raise IOError('Unsupported kind: %s' % kind)
    return result


def dump_chk(filename, data):
    '''Dump a checkpoint file

       Argument:
        | filename  --  the file to write to
        | data  -- a dictionary whose keys are field labels and the values can
                   be None, string, integer, float, boolean, an array/list of
                   strings, integers, floats or booleans.

       The file format is similar to the Gaussian fchk format, but has the extra
       feature that the shapes of the arrays are also stored.
    '''
    with open(filename, 'w') as f:
        for key, value in sorted(data.items()):
            if not isinstance(key, str):
                raise TypeError('The keys must be strings.')
            if len(key) > 40:
                raise ValueError('Key strings can not be longer than 40 characters.')
            if isinstance(value, str):
                if len(value) > 256:
                    raise TypeError('Only small strings are supported (256 chars).')
                if '\n' in value:
                    raise ValueError('The string can not contain new lines.')
                print('%40s  kind=str   %s' % (key.ljust(40), value), file=f)
            elif isinstance(value, bool):
                print('%40s  kind=bln   %s' % (key.ljust(40), value), file=f)
            elif isinstance(value, (int, np.integer)):
                print('%40s  kind=int   %i' % (key.ljust(40), value), file=f)
            elif isinstance(value, float):
                print('%40s  kind=flt   %22.15e' % (key.ljust(40), value), file=f)
            elif isinstance(value, np.ndarray) or isinstance(value, list) or \
                 isinstance(value, tuple):
                if isinstance(value, list) or isinstance(value, tuple):
                    value = np.array(value)
                if value.dtype.fields is not None:
                    raise TypeError('Arrays with fields are not supported.')
                shape_str = ','.join(str(i) for i in value.shape)
                if issubclass(value.dtype.type, (str, np.unicode)):
                    for cell in value.flat:
                        if len(cell) >= 22:
                            raise ValueError('In case of string arrays, a string may contain at most 21 characters.')
                        if ' ' in cell or '\n' in cell:
                            raise ValueError('In case of string arrays, a string may not contain spaces or new lines.')
                    print('%40s  kind=strar %s' % (key.ljust(40), shape_str), file=f)
                    format_str = '%22s'
                elif issubclass(value.dtype.type, np.integer):
                    print('%40s  kind=intar %s' % (key.ljust(40), shape_str), file=f)
                    format_str = '%22i'
                elif issubclass(value.dtype.type, np.bool_):
                    print('%40s  kind=blnar %s' % (key.ljust(40), shape_str), file=f)
                    format_str = '%22s'
                elif issubclass(value.dtype.type, float):
                    print('%40s  kind=fltar %s' % (key.ljust(40), shape_str), file=f)
                    format_str = '%22.15e'
                else:
                    raise TypeError('Numpy array type %s not supported.' % value.dtype.type)
                short_len = 4
                short = []
                for x in value.ravel():
                    short.append(x)
                    if len(short) == 4:
                        print(' '.join(format_str  % s for s in short), file=f)
                        short = []
                if len(short) > 0:
                    print(' '.join(format_str  % s for s in short), file=f)
            elif value is None:
                print('%40s  kind=none   None' % key.ljust(40), file=f)
            else:
                raise TypeError('Type %s not supported.' % type(value))
