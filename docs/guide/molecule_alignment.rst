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
 
.. py:function:: align_mol_skeletons(mol, ref_mol)
 
   Aligns a molecule to a reference by superimposing their shared
   **Maximum Common Substructure (MCS)**.
 
   This function uses RDKit's :mod:`rdFMCS` module to identify the largest
   common subgraph between the two molecules, then calls
   :func:`~rdkit.Chem.AllChem.AlignMol` to minimise the RMSD over the matched
   atom pairs. Because the atom-type comparator is set to
   ``CompareAny``, the MCS search matches atoms regardless of their element,
   making it tolerant of scaffold hops and heteroatom substitutions.
 
   :param mol: The molecule to be aligned (moved). A deep copy is made
               internally, so the caller's object is not affected.
   :type mol: :class:`~gufe.SmallMoleculeComponent`
   :param ref_mol: The reference molecule that provides the target coordinates.
   :type ref_mol: :class:`~gufe.SmallMoleculeComponent`
   :returns: An aligned copy of *mol* superimposed on *ref_mol* via the MCS.
   :rtype: :class:`~gufe.SmallMoleculeComponent`
 
   **When to use this algorithm**
 
   Choose skeleton-based alignment when:
 
   * The two molecules share a recognisable common scaffold (e.g. an
     R-group optimisation series).
   * You want the alignment to reflect chemical similarity rather than
     overall 3D shape.
   * The molecules differ primarily in peripheral substituents and the
     core geometry should be conserved.
 
   **Algorithm outline**
 
   1. Find the MCS of the two molecules using ``CompareAny`` atom typing
      (topology-only matching).
   2. Convert the MCS SMARTS pattern to an atom-index mapping between
      the two molecules.
   3. Call :func:`~rdkit.Chem.AllChem.AlignMol` with the explicit atom map to
      minimise the RMSD over the MCS atoms only.
 
   **Example**
 
   .. code-block:: python
 
      from kartograf.alignment import align_mol_skeletons
      from gufe import SmallMoleculeComponent
 
      mol     = SmallMoleculeComponent.from_sdf("ligand.sdf")
      ref_mol = SmallMoleculeComponent.from_sdf("reference.sdf")
 
      aligned_mol = align_mol_skeletons(mol, ref_mol)
 
   .. warning::
 
      If the two molecules share no common substructure the MCS will be
      empty and the alignment will be undefined. Pre-check your molecules
      if a shared scaffold cannot be assumed.
 
----
 
.. _align-mol-shape:
 
Shape-Based Alignment (Open3DAlign)
-------------------------------------
 
.. py:function:: align_mol_shape(mol, ref_mol)
 
   Aligns a molecule to a reference by maximising the **overlap of their
   3D shapes** using the Open3DAlign (O3A) algorithm.
 
   This function wraps RDKit's :func:`~rdkit.Chem.rdMolAlign.GetO3A`, which
   scores alignment quality using a combination of atom-pair distances and
   partial-charge similarities. The alignment is purely geometry-driven and
   does not require any shared substructure.
 
   :param mol: The molecule to be aligned (moved). A deep copy is made
               internally, so the caller's object is not affected.
   :type mol: :class:`~gufe.SmallMoleculeComponent`
   :param ref_mol: The reference molecule that provides the target coordinates.
   :type ref_mol: :class:`~gufe.SmallMoleculeComponent`
   :returns: An aligned copy of *mol* whose 3D shape best overlaps *ref_mol*.
   :rtype: :class:`~gufe.SmallMoleculeComponent`
 
   **When to use this algorithm**
 
   Choose shape-based alignment when:
 
   * The molecules belong to different chemical series (scaffold hops,
     bioisosteric replacements) but are expected to occupy a similar
     binding volume.
   * No obvious MCS exists or the MCS is too small to anchor a reliable
     structural overlay.
   * You wish to compare or cluster molecules by 3D pharmacophoric shape.
 
   **Algorithm outline**
 
   1. Compute the Open3DAlign score between the probe and reference
      molecules using :func:`~rdkit.Chem.rdMolAlign.GetO3A`.
   2. Call :meth:`Align` to apply the optimal rigid-body rotation and
      translation that maximises shape overlap.
   3. The alignment score (a float) is logged at ``DEBUG`` level for
      diagnostic purposes.
 
   **Example**
 
   .. code-block:: python
 
      from kartograf.alignment import align_mol_shape
      from gufe import SmallMoleculeComponent
 
      mol     = SmallMoleculeComponent.from_sdf("ligand.sdf")
      ref_mol = SmallMoleculeComponent.from_sdf("reference.sdf")
 
      aligned_mol = align_mol_shape(mol, ref_mol)
 
   .. note::
 
      Open3DAlign requires that both molecules carry 3D coordinates and
      (optionally) partial charges for optimal scoring. Ensure that conformers
      have been generated before calling this function.
 
----
 
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
 
For a congeneric series where you want to respect the common core, use
:ref:`align_mol_skeletons <align-mol-skeletons>`. For structurally diverse
molecules where binding-volume overlap is the primary concern, use
:ref:`align_mol_shape <align-mol-shape>`.
 
----
 
See Also
--------
 
* `RDKit MCS documentation <https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html>`_
* `RDKit Open3DAlign documentation <https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html>`_
* `Kartograf repository <https://github.com/OpenFreeEnergy/kartograf>`_
