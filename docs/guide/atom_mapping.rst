Atom Mapping
===============================
The Atom mapper in Kartograf is can be used for hybrid topology approaches in relative free energy calculations.
The focus of this algorithm is, to use the given input coordinates and use them a given truth, that shall only be minimally modified.
This allows users to give pre-modelled ligand alignments to the mapping algorithm and always getting a mapping of the molecules,
that is not changing the input structures significantly.

This approach is very robust, but critically depends on the input structures of course.

Additionally can rules be applied during the mapping algorithm, such that for example hydrogens are only mapped on hydrogens or no ring breakage happens in the mappped core.

(More here soon)

Approach:

.. image:: ../_static/img/Kartograf_mapping_approach.png


Atom Mapping Filters
---------------------