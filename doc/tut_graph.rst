Working with molecular graphs
=============================

.. highlight:: python
   :linenothreshold: 5

Edges
~~~~~

The first example shows how a molecular graph can be created. The most important
attribute of a graph object is ``edges``. This typically is a list of bonds in the
molecule.

File: ``examples/002_graphs/a_graphs.py``

.. literalinclude:: ../examples/002_graphs/a_graphs.py

From the edges many related data structures can be derived. They are accessible
as attributes of a graph object and are only constructed the first time they
are used. A few examples, e.g. ``neighbors``, ``distances``, or ``symmetries``.

Neighbors
~~~~~~~~~

The ``neighbors`` attribute is a dictionary that relates each vertex, i.e atom,
in the graph with its direct neighbors. It can be used to get the direct
chemical environment of an atom, or just to count the number of bonds.

File: ``examples/002_graphs/b_neighbors.py``

.. literalinclude:: ../examples/002_graphs/b_neighbors.py

Distances
~~~~~~~~~

The ``distances`` attribute is a matrix with the number of bonds between two
atoms. The matrix is constructed with the Floyd-Warshall algorithm, which is
implemented in C to make it fast.

File: ``examples/002_graphs/c_distances.py``

.. literalinclude:: ../examples/002_graphs/c_distances.py


Symmetries
~~~~~~~~~~

All isomorphisms of a molecular graph can be requested. For molecular graphs,
only those isomorphisms are considered that do not alter the atomic numbers.

File: ``examples/002_graphs/d_symmetries.py``

.. literalinclude:: ../examples/002_graphs/d_symmetries.py
