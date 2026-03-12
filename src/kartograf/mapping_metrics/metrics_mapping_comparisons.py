from gufe import AtomMapping


def jaccard_score(mappingA: AtomMapping, mappingB: AtomMapping) -> float:
    """
    Calculate the Jaccard score for mapping diversity.

    The Jaccard score is a normalized score ([0,1]) that gives insight
    into the selected atom pair diversity of two compared mappings.
    Diversity is expressed here as the change in atom assignments.
    If the score is 1, the mappings are identical. If the score is 0,
    the two mappings consist of entirely different atom mapping pairs.

    Parameters
    ----------
    mappingA : AtomMapping
        The first atom mapping.
    mappingB : AtomMapping
        The second atom mapping.

    Returns
    -------
    float
        The calculated Jaccard Score [0, 1].
    """

    mappingA_pairs = set(mappingA.componentA_to_componentB.items())
    mappingB_pairs = set(mappingB.componentA_to_componentB.items())

    if len(mappingA_pairs) == 0:
        raise ValueError(f"Mapping A does not contain any mapped atoms: {mappingA_pairs}")
    if len(mappingB_pairs) == 0:
        raise ValueError(f"Mapping B does not contain any mapped atoms: {mappingB_pairs}")

    intersection = mappingA_pairs.intersection(mappingB_pairs)
    union = mappingA_pairs.union(mappingB_pairs)

    return len(intersection) / len(union)
