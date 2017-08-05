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

Databases
=========

The modules :mod:`molmod.periodic`, :mod:`molmod.bonds` and
:mod:`molmod.isotopes` are not loaded automatically with the import statement
``from molmod import *``, but have to be imported explicitely due to their long
loading time, i.e. ::

  >>> from molmod.periodic import periodic
  >>> from molmod.bonds import bonds
  >>> from molmod.isotopes import ame2003, nubtab03

:mod:`molmod.periodic` -- The periodic system
---------------------------------------------

.. automodule:: molmod.periodic
   :members:

:mod:`molmod.bonds` -- The bond length database
-----------------------------------------------

.. automodule:: molmod.bonds
   :members:

:mod:`molmod.isotopes` -- The isotope database
----------------------------------------------

.. automodule:: molmod.isotopes
   :members:
