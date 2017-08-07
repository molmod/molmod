..
    : MolMod is a collection of molecular modelling tools for python.
    : Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
    : for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
    : reserved unless otherwise stated.
    :
    : This file is part of MolMod.
    :
    : MolMod is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : MolMod is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Release history
###############

**Version 1.4.1** August 7, 2017

- Fix windows compatibility bug.

**Version 1.4.0** August 6, 2017

- Testing on Windows instances with AppVeyor, with deployment of Windows packages to
  anaconda.
- Testing on Travis with OSX instances, with deployment of OSX packages to
  anaconda.
- Fix tests failing due to unfortunate random numbers.
- Skip tests on known weaknesses of binning.

**Version 1.3.5** August 5, 2017

- Python 3 bug fix

**Version 1.3.4** August 5, 2017

- Python 3 bug fixes

**Version 1.3.2** August 3, 2017

- Specify versions of dependencies in setup.py.

**Version 1.3.1** August 3, 2017

- Fix parallel testing issue with ``tmpdir`` contextwrapper.

**Version 1.3.0** August 3, 2017

- Python 3 support.

**Version 1.2.1** August 2, 2017

- Switch to setuptools in setup.py.
- Use Cython to compile the extension instead of f2py.
- Use Conda for testing on Travis-CI.
- Automatic deployment of new version on PyPI, Github and anaconda.org
- Simplified installation.
- Many cleanups.
- Fix mistake in the Kabsch algorithm.

**Version 1.1** September 9, 2014

- DLPoly history reader also works for restarted calculation.
- Testing on Travis-CI
- Various small improvements

**Version 1.0** June 5, 2013

- First stable release of MolMod on Github.
