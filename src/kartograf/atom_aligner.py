# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging
from copy import deepcopy

from gufe import SmallMoleculeComponent
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

logger = logging.getLogger(__name__)


def align_mol_skeletons(
    mol: SmallMoleculeComponent,
    ref_mol: SmallMoleculeComponent,
) -> SmallMoleculeComponent:
    """
    Aligns a molecule to a reference by superimposing their shared **Maximum Common Substructure (MCS)**.

    This function uses RDKit's :mod:`rdkit.Chem.rdFMCS` module to identify the largest common subgraph between the two molecules, then calls :func:`~rdkit.Chem.rdMolAlign.AlignMol` to minimise the RMSD over the matched atom pairs.
    Because the atom-type comparator is set to :attr:`~rdkit.Chem.rdFMCS.AtomCompare.CompareAny` the MCS search matches atoms regardless of their element, making it tolerant of scaffold hops and heteroatom substitutions.

    :param mol: The molecule to be aligned (moved).
    :param ref_mol: The reference molecule that provides the target coordinates.
    :returns: An aligned copy of ``mol`` superimposed on ``ref_mol`` via the MCS.

    **When to use this algorithm**

    Choose skeleton-based alignment when:

    * The two molecules share a recognisable common scaffold (e.g. an R-group optimisation series).
    * You want the alignment to reflect chemical similarity rather than overall 3D shape.
    * The molecules differ primarily in peripheral substituents and the core geometry should be conserved.

    **Algorithm outline**

    1. Find the MCS of the two molecules using :attr:`~rdkit.Chem.rdFMCS.AtomCompare.CompareAny` atom typing (topology-only matching).
    2. Convert the MCS SMARTS pattern to an atom-index mapping between the two molecules.
    3. Call :func:`rdkit.Chem.rdMolAlign.AlignMol` with the explicit atom map to minimise the RMSD over the MCS atoms only.

    **Example**

    .. invisible-code-block: python

       # We do this so we can keep our examples easy to read
       # We copy test data into a temp dir so we can run our
       # example in an isolated dir + rename the files
       import shutil, os, tempfile
       _tmp = tempfile.mkdtemp()
       _orig_cwd = os.getcwd()
       shutil.copy("src/kartograf/tests/data/lig_E29.sdf", _tmp + "/ligand.sdf")
       shutil.copy("src/kartograf/tests/data/lig_E6.sdf",  _tmp + "/reference.sdf")
       os.chdir(_tmp)

    .. code-block:: python

       from kartograf import align_mol_skeletons
       from gufe import SmallMoleculeComponent

       mol = SmallMoleculeComponent.from_sdf_file("ligand.sdf")
       ref_mol = SmallMoleculeComponent.from_sdf_file("reference.sdf")

       aligned_mol = align_mol_skeletons(mol, ref_mol)

    .. invisible-code-block: python

       os.chdir(_orig_cwd)  # restore cwd

    .. warning::

       If the two molecules share no common substructure the MCS will be empty and the alignment will be undefined.
       Pre-check your molecules if a shared scaffold cannot be assumed.
    """

    mol = deepcopy(mol)

    mol1b = ref_mol._rdkit
    mol2b = mol._rdkit

    # MCS
    p = rdFMCS.MCSParameters()
    p.AtomTyper = rdFMCS.AtomCompare.CompareAny

    res = rdFMCS.FindMCS([mol1b, mol2b], p)

    # convert match to mapping
    q = Chem.MolFromSmarts(res.smartsString)
    logging.debug(q)

    m1_idx = mol1b.GetSubstructMatch(q)
    m2_idx = mol2b.GetSubstructMatch(q)
    logging.debug(m1_idx, m2_idx)

    idx_mappings = list(zip(m2_idx, m1_idx))

    rms = AllChem.AlignMol(
        prbMol=mol2b,
        refMol=mol1b,
        atomMap=idx_mappings,
    )
    logging.debug(rms)

    mol._rdkit = mol2b
    return mol


def align_mol_shape(mol: SmallMoleculeComponent, ref_mol: SmallMoleculeComponent) -> SmallMoleculeComponent:
    """
    Aligns a molecule to a reference by maximising the **overlap of their 3D shapes** using the Open3DAlign (O3A) algorithm.

    This function wraps RDKit's :func:`~rdkit.Chem.rdMolAlign.GetO3A`, which scores alignment quality using a combination of atom-pair distances and partial-charge similarities.
    The alignment is purely geometry-driven and does not require any shared substructure.

    :param mol: The molecule to be aligned (moved).
    :param ref_mol: The reference molecule that provides the target coordinates.
    :returns: An aligned copy of ``mol`` whose 3D shape best overlaps ``ref_mol``.

    **When to use this algorithm**

    Choose shape-based alignment when:

    * The molecules belong to different chemical series (scaffold hops, bioisosteric replacements) but are expected to occupy a similar binding volume.
    * No obvious MCS exists or the MCS is too small to anchor a reliable structural overlay.
    * You wish to compare or cluster molecules by 3D pharmacophoric shape.

    **Algorithm outline**

    1. Compute the Open3DAlign score between the probe and reference molecules using :func:`~rdkit.Chem.rdMolAlign.GetO3A`.
    2. Call :meth:`Align` to apply the optimal rigid-body rotation and translation that maximises shape overlap.
    3. The alignment score (a float) is logged at ``DEBUG`` level for diagnostic purposes.

    **Example**

    .. invisible-code-block: python

       # We do this so we can keep our examples easy to read
       # We copy test data into a temp dir so we can run our
       # example in an isolated dir + rename the files
       import shutil, os, tempfile
       _tmp = tempfile.mkdtemp()
       _orig_cwd = os.getcwd()
       shutil.copy("src/kartograf/tests/data/lig_E29.sdf", _tmp + "/ligand.sdf")
       shutil.copy("src/kartograf/tests/data/lig_E6.sdf",  _tmp + "/reference.sdf")
       os.chdir(_tmp)

    .. code-block:: python

       from kartograf import align_mol_shape
       from gufe import SmallMoleculeComponent

       mol = SmallMoleculeComponent.from_sdf_file("ligand.sdf")
       ref_mol = SmallMoleculeComponent.from_sdf_file("reference.sdf")

       aligned_mol = align_mol_shape(mol, ref_mol)

    .. invisible-code-block: python

       os.chdir(_orig_cwd )  # restore cwd

    .. note::

       Open3DAlign requires that both molecules carry 3D coordinates and (optionally) partial charges for optimal scoring.
       Ensure that conformers have been generated before calling this function.
    """
    mol = deepcopy(mol)

    mol1b = ref_mol._rdkit
    mol2b = mol._rdkit
    pyO3A = rdMolAlign.GetO3A(
        prbMol=mol2b,
        refMol=mol1b,
    )
    score = pyO3A.Align()
    logging.debug(f"alignment score: {score}")

    mol._rdkit = mol2b
    return mol
