===============================
Atom Mapping
===============================

Atom maappings can have a large variety of different use cases. The general
goal of an atom mapping is to find an assignment between two sets of atoms
based on a certain motivation. A very common approach is to find an
assignment of atoms, between two molecules, which are considers similar/equal
leading to an MCS estimation.

For finding such atom mappings multiple different approaches were described
in literature. One is to use 2D isomorphic graph problem solving in order to
estimate the MCS. Alternatively one can use the 3D geometry information of
the atoms in order to find such a mapping, like Kartograf does with its atom
mapper. Kartograf assumes, that the input coordinates of the molecules are
highly overlapping due to prior alignment or modelling approaches, before the
mapping stage is perfomed. This allows Kartograf to find very efficiently an
atom mapping with a minimal atom displacement.

In the case of classical atom mappings for hybrid topology free energy
calculations, a core region of the molecules is required. That region is
actually a sub space of the actual atom mapping and it ensures, that all
mapped atoms are connected via covalent bonds in their original molecule, as
this otherwise might lead to problems in the sampling during the simulations
(communicating mapped regions).

Approach:

.. image:: ../_static/img/Kartograf_mapping_approach.png


Atom Mapping Filters
---------------------

Additionally can rules be applied during Kartograf's mapping algorithm,
such that for example hydrogens are only mapped on hydrogens or no ring
breakage occurs in the mappped core. Such rules might be necessary in order to
ensure sampling of physical relevant configuraitons or serve other purposes.



