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

        if not (at_i.IsInRing() and at_j.IsInRing()):
            # only applicable for atoms in rings
            # other cases are handled with other filters
            filtered_mapping[i] = j
            continue

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
    """Prevent any ring cleaving transformations in the mapping

    This filter prevents any non-ring atom turning into a ring atom (or
    vice versa)
    """
    filtered_mapping = {}

    for i, j in mapping.items():
        at_i = molA.GetAtomWithIdx(i)
        at_j = molB.GetAtomWithIdx(j)
        if at_i.IsInRing() ^ at_j.IsInRing():
            continue
        filtered_mapping[i] = j

    return filtered_mapping


def filter_whole_rings_only(molA, molB,
                            mapping: dict[int, int]) -> dict[int, int]:
    """Ensure that any mapped rings are wholly mapped"""
    proposed_mapping = {**mapping}

    for mol in [molA, molB]:  # loop over A->B and B->A directions
        ri = mol.GetRingInfo().AtomRings()  # gives list of tuple of atom indices
        ok: dict[int, bool] = {}  # if a ring atom, is this ok to keep?
        for ring in ri:
            # for each ring, atoms are ok if all the atoms of this ring are in
            # the mapping
            all_in_ring = all(atom in proposed_mapping for atom in ring)
            for atom in ring:
                ok[atom] = all_in_ring

        filtered_mapping = {}
        for i, j in proposed_mapping.items():
            if ok.get(i, True):  # if not present, isn't in ring so keep
                filtered_mapping[i] = j

        # reverse the mapping to check B->A (then reverse again)
        proposed_mapping = {v: k for k, v in filtered_mapping}

    return proposed_mapping
