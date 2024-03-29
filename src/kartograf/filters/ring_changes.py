# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from collections import defaultdict
import logging
from rdkit import Chem

logger = logging.getLogger(__name__)


def filter_ringsize_changes(
    molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
) -> dict[int, int]:
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
        rs_i = set(riA.AtomRingSizes(i))
        rs_j = set(riB.AtomRingSizes(j))

        # if there's any intersection in ring size, we're ok
        if not (rs_i & rs_j):
            # e.g. {5} and {6} we've mutated a ring, don't include
            continue
        filtered_mapping[i] = j

    return filtered_mapping


def filter_ringbreak_changes(
    molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
) -> dict[int, int]:
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


def filter_whole_rings_only(
    molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
) -> dict[int, int]:
    """Ensure that any mapped rings are wholly mapped"""
    proposed_mapping = {**mapping}

    for mol in [molA, molB]:  # loop over A->B and B->A directions
        ri = (
            mol.GetRingInfo().AtomRings()
        )  # gives list of tuple of atom indices
        # for each ring, are we fully present?
        ringok: dict[frozenset[int], bool] = {}
        # for each atom, which (maybe unmapped) rings is it present in?
        atom2ring: dict[int, list[frozenset[int]]] = defaultdict(list)
        for ring in ri:
            ringset = frozenset(ring)

            ringok[ringset] = all(atom in proposed_mapping for atom in ring)
            for atom in ring:
                atom2ring[atom].append(ringset)

        filtered_mapping = {}
        for i, j in proposed_mapping.items():
            # if not in any rings, we're ok
            if i not in atom2ring:
                filtered_mapping[i] = j
                continue

            # if in any rings, at least one must be ok
            # e.g. if on edge of fused rings, one ring being completely mapped is ok
            if any(ringok[r] for r in atom2ring[i]):
                filtered_mapping[i] = j

        # reverse the mapping to check B->A (then reverse again)
        proposed_mapping = {v: k for k, v in filtered_mapping.items()}

    return proposed_mapping
