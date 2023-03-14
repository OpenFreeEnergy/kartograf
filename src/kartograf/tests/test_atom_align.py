# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest

from rdkit import Chem

from kartograf.atom_align import align_mol_sceletons
from gufe import SmallMoleculeComponent


@pytest.fixture(scope="session")
def stereo_chem_problem():
    smiles = [
        "C[C@H](F)Br",
        "C[C@@H](F)Br",
    ]

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]
    [Chem.rdDistGeom.EmbedMultipleConfs(m, 1) for m in mols]
    mols = [SmallMoleculeComponent(m) for m in mols]

    return mols


def test_stereo_align(stereo_chem_problem):
    """
    Currently a smoke test
    """
    molA, molB = stereo_chem_problem
    aligned_molA = align_mol_sceletons(molA, molB)
    