import pytest
from kartograf.atom_mapping_scorer import (
    mapping_volume_ratio,
    mapping_rmsd,
    mapping_ratio_of_mapped_atoms,
    mapping_shape_distance,
    default_kartograf_scorer
)

from .conf import stereo_chem_mapping, benzene_benzene_mapping


def test_score_mapping_volume_ratio(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = mapping_volume_ratio()
    score = scorer(stereo_chem_mapping)
    print(score)


def test_score_mappings_rmsd(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = mapping_rmsd()
    score = scorer.get_rmsd(stereo_chem_mapping)
    print(score)


def test_score_norm_mapping_rmsd(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = mapping_rmsd()
    score = scorer.get_rmsd_p(stereo_chem_mapping)
    print(score)



def test_score_tanimoto_shape_identical(benzene_benzene_mapping):
    """
    Currently a smoke test
    """
    scorer = mapping_shape_distance()
    score = scorer(benzene_benzene_mapping)

    assert isinstance(score, float)
    assert score == 0, "Score of identical Molecule should be 0"


def test_default_kartograph_scorer_identical(benzene_benzene_mapping):
    """
    Currently a smoke test
    """
    scorer = default_kartograf_scorer()
    score = scorer(benzene_benzene_mapping)

    assert isinstance(score, float)
    assert score == 0, "Score of identical Molecule should be 0"

def test_default_kartograph_scorer_identical(benzene_benzene_mapping):
    """
    Currently a smoke test
    """
    scorer = mapping_ratio_of_mapped_atoms()
    score = scorer(benzene_benzene_mapping)

    assert isinstance(score, float)
    assert score == 0, "Score of identical Molecule should be 0"