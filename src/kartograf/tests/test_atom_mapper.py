# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest

from rdkit import Chem

from kartograf import kartograf_atom_mapper
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
    Chem.rdMolAlign.AlignMol(mols[0], mols[1])

    return mols


def test_stereo_mapping(stereo_chem_problem):
    """
    Currently a smoke test
    """
    expected_solution = {}
    geom_mapper = kartograf_atom_mapper(
        atom_max_distance=0.95, atom_map_hydrogens=True
    )  # mapping_algorithm.minimal_spanning_tree
    geom_mapping = next(
        geom_mapper.suggest_mappings(
            SmallMoleculeComponent(stereo_chem_problem[0]),
            SmallMoleculeComponent(stereo_chem_problem[1]),
        )
    )
