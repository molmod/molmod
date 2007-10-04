#!/usr/bin/env python
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
# --

if __name__ == "__main__":
    import glob
    from numpy.distutils.core import setup
    from numpy.distutils.extension import Extension

    helpers = Extension(
        name='molmod.helpers',
        sources=['src/helpers.f90']
    )
    setup(
        name='MolMod',
        version='0.001',
        description='MolMod is a collection of molecular modelling tools for python.',
        author='Toon Verstraelen',
        author_email='Toon.Verstraelen@UGent.be',
        url='https://molmod.ugent.be/zeobuilder/',
        package_dir = {'molmod': 'src'},
        ext_modules=[helpers],
        packages=[
            'molmod',
            'molmod.data',
        ],
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

    from distutils.core import setup
    setup(
        name='MolMod',
        version='0.001',
        description='MolMod is a collection of molecular modelling tools for python.',
        author='Toon Verstraelen',
        author_email='Toon.Verstraelen@UGent.be',
        url='https://molmod.ugent.be/zeobuilder/',
        data_files=[
            ('share/molmod', glob.glob('share/*.csv')),
        ],
    )


