Molecule Alignment Algorithms
=============================

Kartograf provides two complementary algorithms for aligning small molecules in 3D space.
Both functions accept a molecule to be moved and a reference molecule, and return a new, aligned copy of the input molecule without modifying the original.

.. note::

   Both functions operate on :class:`~gufe.SmallMoleculeComponent` objects and return a deep copy of the input molecule with updated 3D coordinates.
   The reference molecule is never modified.

----

.. _align-mol-skeletons:

Skeleton-Based Alignment (MCS)
-------------------------------

.. autofunction:: kartograf.atom_aligner.align_mol_skeletons
   :no-index:

.. _align-mol-shape:

Shape-Based Alignment (Open3DAlign)
-------------------------------------

.. autofunction:: kartograf.atom_aligner.align_mol_shape
   :no-index:

Choosing Between the Two Methods
---------------------------------

+--------------------------------+---------------------------+---------------------------+
| Criterion                      | Skeleton (MCS)            | Shape (O3A)               |
+================================+===========================+===========================+
| Requires shared substructure   | Yes                       | No                        |
+--------------------------------+---------------------------+---------------------------+
| Sensitive to scaffold changes  | Low                       | High                      |
+--------------------------------+---------------------------+---------------------------+
| Best for congeneric series     | ✓                         |                           |
+--------------------------------+---------------------------+---------------------------+
| Best for diverse scaffolds     |                           | ✓                         |
+--------------------------------+---------------------------+---------------------------+
| Alignment driven by            | Atom topology             | 3D volume & charges       |
+--------------------------------+---------------------------+---------------------------+
| Underlying RDKit module        | ``rdFMCS`` / ``AllChem``  | ``rdMolAlign`` (O3A)      |
+--------------------------------+---------------------------+---------------------------+

For a congeneric series where you want to respect the common core, use :ref:`align_mol_skeletons <align-mol-skeletons>`.
For structurally diverse molecules where binding-volume overlap is the primary concern, use :ref:`align_mol_shape <align-mol-shape>`.
