# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from rdkit import Chem
import logging

logger = logging.getLogger(__name__)

def filter_bond_breaks(mol_a: Chem.Mol, mol_b: Chem.Mol, mapping: dict[int, int]):
    """
    Remove cases where a bond would be broken during the transformation.
    These are detected by checking that all bonds of ligand A which are in the mapping also form a bond in ligand B,
    if they do not the atoms are removed from the mapping.

    See <https://github.com/OpenFreeEnergy/kartograf/issues/88> for more details.

    Parameters
    ----------
    mol_a
        Molecule at state A with atom indices in the keys of the mapping dict
    mol_b
        Molecule at state B with atom indices in the values of the mapping dict
    mapping
        The mapping between molecule A and B

    Returns
    -------
        A filtered mapping removing any entries which map a broken bond.
    """

    to_remove = []
    # generate a list of bonds in molecule b
    mol_b_bonds = [sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in mol_b.GetBonds()]

    # loop over the bonds in molecule A and check they are connected in molecule B if they are in the mapping
    for bond in mol_a.GetBonds():
        atom_a = bond.GetBeginAtomIdx()
        atom_b = bond.GetEndAtomIdx()

        if atom_a in mapping and atom_b in mapping:
            # the bond is in the mapped section
            # make sure the mapped atoms are also bonded
            mapped_bond = sorted((mapping[atom_a], mapping[atom_b]))
            if mapped_bond not in mol_b_bonds:
                # if the bond is not in molecule b remove the bond atoms from the mapping
                to_remove.extend([atom_a, atom_b])

    new_mapping = {k: v for k, v in mapping.items() if k not in to_remove}

    return new_mapping
