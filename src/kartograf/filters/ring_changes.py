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


def filter_hybridization_rings(
        molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
) -> dict[int, int]:
    """Ensure that any mapped rings are either both aromatic or aliphatic

    e.g. this filter would unmap hexane to benzene type transformations
    """

    def get_atom_ring_hybridization_map(rdmol: Chem.Mol) -> dict[int, bool]:
        """
        For each atom, determines information about ring hybridization

        Parameters
        ----------
        rdmol: Chem.Mol

        Returns
        -------
        dict[int, bool]:
            returns a dict, that maps each atom's index to if it is always in aromatic rings
            (If True this atom exists entirely in aromatic rings, if False it is
            not entirely aromatic and therefore not necessarily sterically restraint
            by a pi-orbital-system.)
        """
        riInf = rdmol.GetRingInfo()

        # get_ring_hybridization
        # for each ring, the ring is aromatic if all atoms within are aromatic
        # maps ring index to aromaticity as bool
        is_ring_aromatic = {}
        for i, ring in enumerate(riInf.BondRings()):
            is_ring_aromatic[i] = all(
                rdmol.GetBondWithIdx(atomI).GetIsAromatic()
                for atomI in ring)

        # first iterate over all rings and determine if they are aromatic
        # map atoms to ring aromaticities
        atom_ring_map = defaultdict(list)
        for ri, r in enumerate(riInf.AtomRings()):
            for a in r:
                atom_ring_map[a].append(is_ring_aromatic[ri])

        # then with all rings traversed, crush this information down to a single bool per atom
        # maps atom index to all ring aromaticity
        atom_aromatic = {}
        for a, v in atom_ring_map.items():
            atom_aromatic[a] = all(v)

        return atom_aromatic

    atomA_ring_hyb_map = get_atom_ring_hybridization_map(molA)
    atomB_ring_hyb_map = get_atom_ring_hybridization_map(molB)

    # Filtering Mapping
    filtered_mapping = {}
    for ai, aj in mapping.items():
        ai_only_arom_sys = atomA_ring_hyb_map[ai]
        aj_only_arom_sys = atomB_ring_hyb_map[aj]

        if ai_only_arom_sys == aj_only_arom_sys:
            filtered_mapping[ai] = aj

    return filtered_mapping
