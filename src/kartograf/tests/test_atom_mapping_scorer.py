import pytest
from kartograf.atom_mapping_scorer import DefaultKartografScorer
from kartograf.mapping_metrics import (
    MappingRMSDScorer,
    MappingShapeMismatchScorer,
    MappingShapeOverlapScorer,
    MappingVolumeRatioScorer,
    MappingRatioMappedAtomsScorer,

)

from kartograf.mapping_metrics.metric_shape_difference import _MappingShapeDistanceScorer
from .conf import stereo_chem_mapping, benzene_benzene_mapping, benzene_benzene_empty_mapping


def test_score_mappings_rmsd(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = MappingRMSDScorer()
    score = scorer.get_rmsd(stereo_chem_mapping)
    print(score)


def test_score_norm_mapping_rmsd(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = MappingRMSDScorer()
    score = scorer.get_rmsd_p(stereo_chem_mapping)
    print(score)

def test_score_mapping_volume_ratio(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = MappingVolumeRatioScorer()
    score = scorer(stereo_chem_mapping)
    print(score)


def test_score_shape_dist(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = _MappingShapeDistanceScorer()
    score = scorer(stereo_chem_mapping)
    print(score)

def test_score_shape_overlap(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = MappingShapeOverlapScorer()
    score = scorer(stereo_chem_mapping)
    print(score)

def test_score_shape_mismatch(stereo_chem_mapping):
    """
    Currently a smoke test
    """
    scorer = MappingShapeMismatchScorer()
    score = scorer(stereo_chem_mapping)
    print(score)


@pytest.mark.parametrize("scorer_class", [ MappingRMSDScorer,
    _MappingShapeDistanceScorer,
    MappingShapeMismatchScorer,
    MappingShapeOverlapScorer,
    MappingVolumeRatioScorer,
    MappingRatioMappedAtomsScorer,
    DefaultKartografScorer])
def test_scorer_identical_molecules(scorer_class, benzene_benzene_mapping):
    """
    Currently a smoke test
    """
    scorer = scorer_class()
    score = scorer(benzene_benzene_mapping)

    assert isinstance(score, float)
    assert score == 0, "Score of identical Molecule should be 0"


@pytest.mark.parametrize("scorer_class, exp", [(MappingRMSDScorer, 0),
                                               (_MappingShapeDistanceScorer, 1),
                                               (MappingShapeMismatchScorer,1),
                                               (MappingShapeOverlapScorer, 1),
                                               (MappingRatioMappedAtomsScorer,1)])
def test_scorer_empty_mapping(scorer_class, exp:float, benzene_benzene_empty_mapping):
    """
    Currently a smoke test
    """
    scorer = scorer_class()
    score = scorer(benzene_benzene_empty_mapping)

    assert isinstance(score, float)
    assert score == exp, "Score of non-Mapping should be "+str(exp)+" for "+str(scorer_class.__name__)+" but is: "+str(score)


@pytest.mark.parametrize("scorer_class, exp", [(MappingVolumeRatioScorer, 1),
                                               (DefaultKartografScorer,1)])
def test_scorer_empty_mapping_err(scorer_class, exp:float, benzene_benzene_empty_mapping):
    """
    Currently a smoke test
    """
    with pytest.raises(ValueError) as exc:
        scorer = scorer_class()
        score = scorer(benzene_benzene_empty_mapping)

    assert "Mapping is too small to calculate convex hull" in str(exc.value)