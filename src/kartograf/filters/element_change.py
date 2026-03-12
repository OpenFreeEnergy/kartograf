# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging

from rdkit import Chem

logger = logging.getLogger(__name__)


def filter_atoms_h_only_h_mapped(molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]) -> dict[int, int]:
    """
    Forces a mapping to only allow hydrogens to map to other hydrogens.

    Parameters
    ----------
    molA : Chem.Mol
        The first molecule.
    molB : Chem.Mol
        The second molecule.
    mapping : dict[int, int]
        The initial mapping between atoms of molA and molB.

    Returns
    -------
    dict[int, int]
        The filtered mapping where only hydrogens are mapped to other hydrogens.
    """

    filtered_mapping = {}
    for atomA_idx, atomB_idx in mapping.items():
        atomA = molA.GetAtomWithIdx(atomA_idx)
        atomB = molB.GetAtomWithIdx(atomB_idx)

        if (atomA.GetAtomicNum() == 1 and atomB.GetAtomicNum() == 1) or (
            atomA.GetAtomicNum() != 1 and atomB.GetAtomicNum() != 1
        ):
            filtered_mapping[atomA_idx] = atomB_idx
            logger.debug(
                f"keep mapping for atomIDs ({atomA_idx}, {atomB_idx}): {atomA.GetAtomicNum()} {atomB.GetAtomicNum()}"
            )
        else:
            logger.debug(
                f"no mapping for atomIDs ({atomA_idx}, {atomB_idx}): {atomA.GetAtomicNum()} {atomB.GetAtomicNum()}"
            )

    return filtered_mapping


def filter_element_changes(molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]) -> dict[int, int]:
    """
    Forces a mapping to exclude any alchemical element changes in the core.

    Parameters
    ----------
    molA : Chem.Mol
        The first molecule.
    molB : Chem.Mol
        The second molecule.
    mapping : dict[int, int]
        The initial mapping between atoms of molA and molB.

    Returns
    -------
    dict[int, int]
        The filtered mapping excluding any alchemical element changes.
    """
    filtered_mapping = {}

    for i, j in mapping.items():
        if molA.GetAtomWithIdx(i).GetAtomicNum() != molB.GetAtomWithIdx(j).GetAtomicNum():
            continue
        filtered_mapping[i] = j

    return filtered_mapping


def filter_hybridization_changes(molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]) -> dict[int, int]:
    """
    Forces a mapping to exclude any alchemical atom hybridization changes.
    Such a change could be for example -C-C >> -C=C.

    Parameters
    ----------
    molA : Chem.Mol
        The first molecule.
    molB : Chem.Mol
        The second molecule.
    mapping : dict[int, int]
        The initial mapping between atoms of molA and molB.

    Returns
    -------
    dict[int, int]
        The filtered mapping excluding any alchemical atom hybridization changes.
    """
    filtered_mapping = {}

    for i, j in mapping.items():
        if molA.GetAtomWithIdx(i).GetHybridization() != molB.GetAtomWithIdx(j).GetHybridization():
            continue
        filtered_mapping[i] = j

    return filtered_mapping
