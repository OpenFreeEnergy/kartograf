# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolAlign

from gufe import SmallMoleculeComponent

import logging

logger = logging.getLogger(__name__)


def align_mol_sceletons(
    mol: SmallMoleculeComponent,
    ref_mol: SmallMoleculeComponent,
) -> SmallMoleculeComponent:
    """
        This i a Wrapper for rdkit - MCS align
        Aligns very simply molecule to the reference molecule,
        based on the shared MCS - Sceleton.

    Parameters
    ----------
    mol : SmallMoleculeComponent
        molecule to be aligned to molA (will be moved)
    ref_mol : SmallMoleculeComponent
        molecule with the reference_positions.

    Returns
    -------
    SmallMoleculeComponent
        return an aligned copy of molB
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


def align_mol_shape(
    mol: SmallMoleculeComponent, ref_mol: SmallMoleculeComponent
) -> Chem.Mol:
    """
        This is a Wrapper for rdkit / OPEN3DAlign
        Aligns shape based two SmallMoleculeComponents.

    Parameters
    ----------
    mol : SmallMoleculeComponent
        molecule to be aligned to molA (will be moved)
    ref_mol : SmallMoleculeComponent
        molecule with the reference_positions.

    Returns
    -------
    SmallMoleculeComponent
        return an aligned copy of molB
    """
    mol = deepcopy(mol)

    mol1b = ref_mol._rdkit
    mol2b = mol._rdkit
    pyO3A = rdMolAlign.GetO3A(
        prbMol=mol2b,
        refMol=mol1b,
    )
    score = pyO3A.Align()
    logging.debug("alignment score: " + str(score))

    mol._rdkit = mol2b
    return mol
