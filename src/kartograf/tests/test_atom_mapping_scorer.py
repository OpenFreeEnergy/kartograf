import pytest

from rdkit import Chem
from gufe import SmallMoleculeComponent, LigandAtomMapping

from kartograf.atom_mapping_scorer import (
    mapping_volume_ratio,
    mapping_area_ratio,
    mappings_rmsd,
    norm_mapping_rmsd,
    number_of_mapped_atoms_ratio,
)


@pytest.fixture(scope="session")
def stereo_mapping():
    smiles = [
        "C[C@H](F)Br",
        "C[C@@H](F)Br",
    ]

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]
    [Chem.rdDistGeom.EmbedMultipleConfs(m, 1) for m in mols]
    Chem.rdMolAlign.AlignMol(mols[0], mols[1])

    return LigandAtomMapping(
        SmallMoleculeComponent(mols[0]),
        SmallMoleculeComponent(mols[1]),
        {2: 7, 4: 5, 5: 4, 6: 6, 7: 2, 0: 0, 1: 1, 3: 3},
    )


def test_score_mapping_volume_ratio(stereo_mapping):
    """
    Currently a smoke test
    """
    score = mapping_volume_ratio(stereo_mapping)
    print(score)


def test_score_mappings_rmsd(stereo_mapping):
    """
    Currently a smoke test
    """
    score = mappings_rmsd(stereo_mapping)
    print(score)


def test_score_norm_mapping_rmsd(stereo_mapping):
    """
    Currently a smoke test
    """
    score = norm_mapping_rmsd(stereo_mapping)
    print(score)


def test_score_mapping_area_ratio(stereo_mapping):
    """
    Currently a smoke test
    """
    score = mapping_area_ratio(stereo_mapping)
    print(score)


def test_score_mapping_area_ratio(stereo_mapping):
    """
    Currently a smoke test
    """
    score = number_of_mapped_atoms_ratio(stereo_mapping)
    print(score)
