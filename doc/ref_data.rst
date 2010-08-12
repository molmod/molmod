Databases
=========

The modules :mod:`molmod.periodic`, :mod:`molmod.bonds` and
:mod:`molmod.isotopes` are not loaded automatically with the import statement
``from molmod import *``, but have to be imported explicitely due to their long
loading time, i.e. ::

  >>> from molmod.periodic import periodic
  >>> from molmod.bonds import bonds
  >>> from molmod.isotopes import ame2003, nubtab03

The periodic system
-------------------

.. automodule:: molmod.periodic
   :members:

The bond length database
------------------------

.. automodule:: molmod.bonds
   :members:

The isotope database
--------------------

This is an extension of the Graph object with molecular features

.. automodule:: molmod.isotopes
   :members:
