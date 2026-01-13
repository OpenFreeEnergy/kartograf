# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest
from gufe import LigandAtomMapping

from kartograf.mapping_metrics.metrics_mapping_comparisons import jaccard_score


def test_mapping_comparison_jcs_identical(benzene_benzene_mapping):
    """
    Check identical mappings
    """
    mapping = benzene_benzene_mapping

    score = jaccard_score(mappingA=mapping, mappingB=mapping)

    assert score == 1


def test_mapping_comparison_jcs_diverse(benzene_benzene_mapping):
    """
    Check completly different mappings
    """
    # Totallydifferent indices
    mapping = benzene_benzene_mapping
    mappingA = LigandAtomMapping(
        componentA=mapping.componentA,
        componentB=mapping.componentB,
        componentA_to_componentB={v: v for v in range(0, 6)},
    )
    mappingB = LigandAtomMapping(
        componentA=mapping.componentA,
        componentB=mapping.componentB,
        componentA_to_componentB={v: v for v in range(6, 12)},
    )

    score = jaccard_score(mappingA=mappingA, mappingB=mappingB)

    assert score == 0

    # Check index order change
    r1 = list(range(0, 6))
    r2 = list(range(6, 12))
    index_mixing = list(zip(r1, r2))

    mappingA = LigandAtomMapping(
        componentA=mapping.componentA,
        componentB=mapping.componentB,
        componentA_to_componentB={k: v for k, v in index_mixing},
    )
    mappingB = LigandAtomMapping(
        componentA=mapping.componentA,
        componentB=mapping.componentB,
        componentA_to_componentB={v: k for k, v in index_mixing},
    )
    print(mappingA.componentA_to_componentB, mappingB.componentA_to_componentB)
    score = jaccard_score(mappingA=mappingA, mappingB=mappingB)

    assert score == 0


def test_mapping_comparison_jcs_empty_mapping(benzene_benzene_mapping, benzene_benzene_empty_mapping):
    """
    Check empty mappings
    """
    mappingA = benzene_benzene_empty_mapping
    mappingB = benzene_benzene_mapping
    with pytest.raises(ValueError) as exc:
        jaccard_score(mappingA=mappingA, mappingB=mappingB)

    assert "Mapping A does not contain any mapped atoms: set()" in str(exc.value)

    mappingA = benzene_benzene_mapping
    mappingB = benzene_benzene_empty_mapping
    with pytest.raises(ValueError) as exc:
        jaccard_score(mappingA=mappingA, mappingB=mappingB)

    assert "Mapping B does not contain any mapped atoms: set()" in str(exc.value)
