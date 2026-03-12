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
    Alignment based on MCS.

    This is a wrapper for RDKit - MCS align that aligns one molecule to another based on the shared MCS skeleton.

    Parameters
    ----------
    mol : SmallMoleculeComponent
        Molecule to be aligned to `ref_mol` (will be moved).
    ref_mol : SmallMoleculeComponent
        Molecule with the reference positions.

    Returns
    -------
    SmallMoleculeComponent
        Aligned copy of `mol`.
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


def align_mol_shape(mol: SmallMoleculeComponent, ref_mol: SmallMoleculeComponent) -> Chem.Mol:
    """
    Alignment based on shape.

    This is a wrapper for RDKit / Open3DAlign that aligns the shape of one molecule to another.

    Parameters
    ----------
    mol : SmallMoleculeComponent
        Molecule to be aligned to `ref_mol` (will be moved).
    ref_mol : SmallMoleculeComponent
        Molecule with the reference positions.

    Returns
    -------
    Chem.Mol
        Aligned RDKit molecule.
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
