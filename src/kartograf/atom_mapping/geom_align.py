from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem

from gufe import SmallMoleculeComponent

import logging

log = logging.getLogger(__name__)


def align_mol_sceletons(
    mol: SmallMoleculeComponent,
    ref_mol: SmallMoleculeComponent,
) -> SmallMoleculeComponent:
    """
        Aligns very simply molecule to the reference molecule, based on the shared MCS.

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

    p = rdFMCS.MCSParameters()
    p.AtomTyper = rdFMCS.AtomCompare.CompareAny

    res = rdFMCS.FindMCS([mol1b, mol2b], p)

    # convert match to mapping'
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
