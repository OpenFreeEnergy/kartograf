

import pytest
from kartograf.atom_mapping_scorer import (
    mapping_volume_ratio,
    mapping_area_ratio,
    mappings_rmsd,
    norm_mapping_rmsd,
    number_of_mapped_atoms_ratio,
)

from .conf import stereo_mapping




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
