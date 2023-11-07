import numpy as np

from gufe import AtomMapping


def jaccard_score(mappingA:AtomMapping, mappingB:AtomMapping) -> float:
    """
        The Jaccard score is a normalized score ([0,1]) , that gives insight
        on the selected atom pair diversity of two compared mappings.
        Diversity is expressed here as the change in atom assignments,
        if the score is 1 the mappings are identical. If the score is 0,
        the two mapping consist of entirely different atom mapping pairs.

    Parameters
    ----------
    mappingA: AtomMapping
    mappingB: AtomMapping

    Returns
    -------
    float
        the calculated Jaccard Score [0, 1]
    """

    mappingA_pairs = set(mappingA.componentA_to_componentB.items())
    mappingB_pairs = set(mappingB.componentA_to_componentB.items())

    intersection = mappingA_pairs.intersection(mappingB_pairs)
    union = mappingA_pairs.union(mappingB_pairs)

    return np.round(len(intersection) / len(union), 2)
