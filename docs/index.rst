Welcome to Kartograf's documentation!
=========================================
Kartograf is a packgage for atom mappings focussing on 3D geometries.
This package can for example be used to generate hybrid topology systems, where an atom mapping is required to determine the core region of the approach.
But of course, there exist also other use cases for this package.
The atom mapper takes two sets of coordinates of molecules as input.
Optionally those sets of coordinates can be aligned onto each other, checkout
the ``atom_aligner`` module functions
of Kartograf that offer a shape alignment implementation and an MCS-skeleton alignment.
The ``atom_mapper`` can be used to generate the 3D geometry-focused atom
mapping, the algorithm is described in the related publication of Kartograf (see reference).
Additionally, rule-based filter functions can be provided to demap atoms,
that do not fulfill the desired criteria, see ``filters``.
Several mapping scoring metrics are provided, that evaluate geometric
properties of your mapping, from ``mapping_metrics``, which might be
useful for checking the quality of your mappings.
Finally, there is a visualization function ``display_mappings_3d`` that can be
used to check out the mappings with a Jupyter Notebook widget.

You can find our Preprint on `Ries, B.; Alibay, I.; Swenson, D. W. H; Baumann, H. M.; Henry, M. M.; Eastwood, J. R. B.; Gowers, R. J. - Kartograf: An Accurate Geometry-Based Atom Mapper for Hybrid Topology Relative Free Energy Calculations, Chemrxiv (2023) <https://doi.org/10.26434/chemrxiv-2023-0n1pq>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorial
   guide
   api
   CHANGELOG


