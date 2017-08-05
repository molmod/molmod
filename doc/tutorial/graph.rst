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

Working with molecular graphs
=============================

Introduction
~~~~~~~~~~~~

Molecular graphs are objects that describe the connectivity of atoms in a
molecule. In graph speech, every atom is a vertex and every bond is an edge.
Some background on `Graph Theory <http://en.wikipedia.org/wiki/Graph_theory>`_
may be helpful. The :class:`molmod.molecular_graphs.MolecularGraph` class in
MolMod can store atomic numbers associated with vertices and bond orders
associated with edges. The class :class:`molmod.graphs.Graph` is the base class
that implements all the abstract graph theory without refering to atoms and
bonds.

.. highlight:: python
   :linenothreshold: 5

Examples
~~~~~~~~

Edges
-----

The first example shows how one creates a molecular graph.. The most important
attribute of a graph object is ``edges``. This is (in most cases) a list of
bonds in the molecule, but more edges may be included if that would be useful.

File: ``molmod/examples/002_graphs/a_graphs.py``

.. literalinclude:: ../../molmod/examples/002_graphs/a_graphs.py

From the edges many related properties can be derived. They are accessible
as attributes of a graph object and are only constructed the first time they
are used. A few examples, e.g. ``neighbors``, ``distances``, or ``symmetries``.

Neighbors
---------

The ``neighbors`` attribute is a dictionary that relates each vertex, i.e atom,
in the graph with its direct neighbors. It can be used to get the direct
chemical environment of an atom, or just to count the number of bonds.

File: ``molmod/examples/002_graphs/b_neighbors.py``

.. literalinclude:: ../../molmod/examples/002_graphs/b_neighbors.py

Distances
---------

The ``distances`` attribute is a square matrix with a row and a column for every
vertex. Each element contains the (minimal) number of bonds between two atoms.
The matrix is constructed with the `Floyd-Warshall algorithm
<http://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm>`_.

File: ``molmod/examples/002_graphs/c_distances.py``

.. literalinclude:: ../../molmod/examples/002_graphs/c_distances.py


Symmetries
----------

All isomorphisms of a molecular graph can be requested. For molecular graphs,
only those isomorphisms are considered that do not alter the atomic numbers.

File: ``molmod/examples/002_graphs/d_symmetries.py``

.. literalinclude:: ../../molmod/examples/002_graphs/d_symmetries.py


Problems
~~~~~~~~

* Write a program that prints the atom indexes of all the C-O bonds in the
  caffeine molecule.

* Write a program that prints the atom indexes of all the Nitrogen atoms that
  are not bonded to a methyl group.

* Write a program that finds all Nitrogen atoms in the caffeine that are at
  least three bonds away from an Oxygen atom.

* Download the `benzene
  <http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=241&disopt=3DSaveSDF>`_
  molecule from Pubchem and write a program that prints out all the graph
  isomorphisms as permutations. The SDF files from pubchem already contain the
  bond graph, so there is no need for ``mol.set_default_graph()`` after loading
  the molecule. How is the number of isomorphisms related to the rotational
  symmetry number of benzene?
