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
