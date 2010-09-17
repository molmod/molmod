Working with molecular graphs
=============================

Molecular graphs are objects that describe the connectivity of atoms in a
molecule. In graph speech, every vertex is an atom and every edge is a bond.
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

The first example shows how a molecular graph can be created. The most important
attribute of a graph object is ``edges``. This typically is a list of bonds in the
molecule.

File: ``examples/002_graphs/a_graphs.py``

.. literalinclude:: ../examples/002_graphs/a_graphs.py

From the edges many related data structures can be derived. They are accessible
as attributes of a graph object and are only constructed the first time they
are used. A few examples, e.g. ``neighbors``, ``distances``, or ``symmetries``.

Neighbors
---------

The ``neighbors`` attribute is a dictionary that relates each vertex, i.e atom,
in the graph with its direct neighbors. It can be used to get the direct
chemical environment of an atom, or just to count the number of bonds.

File: ``examples/002_graphs/b_neighbors.py``

.. literalinclude:: ../examples/002_graphs/b_neighbors.py

Distances
---------

The ``distances`` attribute is a matrix with the number of bonds between two
atoms. The matrix is constructed with the Floyd-Warshall algorithm, which is
implemented in C to make it fast.

File: ``examples/002_graphs/c_distances.py``

.. literalinclude:: ../examples/002_graphs/c_distances.py


Symmetries
----------

All isomorphisms of a molecular graph can be requested. For molecular graphs,
only those isomorphisms are considered that do not alter the atomic numbers.

File: ``examples/002_graphs/d_symmetries.py``

.. literalinclude:: ../examples/002_graphs/d_symmetries.py


Problems
~~~~~~~~

* Write a program that prints the atom indexes of all the C-O bonds in the
  caffeine molecule.

* Write a program that prints the atom indexes of all the Nitrogen atoms that
  are not bonded to a methyl group.

* Write a program that finds all Nitrogen molecules in the caffeine that are at
  least three bonds away from an Oxygen atom.

* Download the `benzene
  <http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=241&disopt=3DSaveSDF>`_
  molecule from Pubchem and write a program that prints out all the graph
  isomorphisms as permutations. Note that the SDF files from pubchem already
  contain the bond graph, so there is no need to call
  :meth:`Molecule.set_default_graph` after loading the molecule. How is the
  number of isomorphisms related to the rotational symmetry number of benzene?
