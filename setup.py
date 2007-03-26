#!/usr/bin/env python
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

if __name__ == "__main__":
    import glob
    from distutils.core import setup


    version = '0.1.2'

    setup(name='MolMod',
        version=version,
        description='MolMod is a collection of molecular modelling tools for python.',
        author='Toon Verstraelen',
        author_email='Toon.Verstraelen@UGent.be',
        url='http://molmod.ugent.be/projects/molmod',
        package_dir = {'molmod': 'src'},
        data_files=[
            ('share/molmod/%s' % version, glob.glob('share/*.csv')),
        ],
        packages=[
            'molmod',
            'molmod.data',
        ],
        scripts=['scripts/mrandomize'],
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
            'Topic :: Science/Engineering :: Molecular Science'
        ],
    )


