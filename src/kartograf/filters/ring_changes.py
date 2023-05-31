import logging

logger = logging.getLogger(__name__)


def filter_ringsize_changes(molA, molB,
                            mapping: dict[int, int]) -> dict[int, int]:
    """Prevents mutating the size of rings in the mapping"""
    riA = molA.GetRingInfo()
    riB = molB.GetRingInfo()

    filtered_mapping = {}

    for i, j in mapping.items():
        at_i = molA.GetAtomWithIdx(i)
        at_j = molB.GetAtomWithIdx(j)

        # AtomRingSizes gives as tuple of ring sizes that this atom is part of
        rs_i = set(riA.AtomRingSizes(at_i))
        rs_j = set(riB.AtomRingSizes(at_j))

        # if there's any intersection in ring size, we're ok
        if not (rs_i & rs_j):
            # e.g. {5} and {6} we've mutated a ring, don't include
            continue
        filtered_mapping[i] = j

    return filtered_mapping


def filter_ringbreak_changes(molA, molB,
                             mapping: dict[int, int]) -> dict[int, int]:
    """Prevent any ring cleaving transformations in the mapping"""
    filtered_mapping = {}

    for i, j in mapping.items():
        at_i = molA.GetAtomWithIdx(i)
        at_j = molB.GetAtomWithIdx(j)
        if at_i.IsInRing() ^ at_j.IsInRing():
            continue
        filtered_mapping[i] = j

    return filtered_mapping
