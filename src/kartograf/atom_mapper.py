# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import copy
import dill
import numpy as np
from enum import Enum

from collections import OrderedDict
from collections.abc import Iterator

from rdkit import Chem

from scipy.sparse import csr_matrix
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import connected_components

from typing import Callable, Iterable, Optional, Union

from gufe import SmallMoleculeComponent, ProteinComponent
from gufe import AtomMapping, AtomMapper, LigandAtomMapping

from numpy.typing import NDArray

import logging

from .filters import (
    filter_atoms_h_only_h_mapped,
    filter_ringbreak_changes,
    filter_ringsize_changes,
    filter_whole_rings_only,
)

logger = logging.getLogger(__name__)


# Enums:
class mapping_algorithm(Enum):
    """
    This enum helps selecting the optimization algorithm
    for the distance graph based atom mapping.
    """

    linear_sum_assignment = "LSA"
    minimal_spanning_tree = "MST"


_mapping_alg_type = Callable[[NDArray, float], dict[int, int]]


# Helper:
vector_eucledean_dist = calculate_edge_weight = lambda x, y: np.sqrt(
    np.sum(np.square(y - x), axis=1)
)


# Implementation of Mapper:
class KartografAtomMapper(AtomMapper):
    atom_max_distance: float
    _map_exact_ring_matches_only: bool
    atom_map_hydrogens: bool
    mapping_algorithm: _mapping_alg_type

    _filter_funcs: list[
        Callable[[Chem.Mol, Chem.Mol, dict[int, int]], dict[int, int]]
    ]

    def __init__(
            self,
            *,
            atom_max_distance: float = 0.95,
            atom_map_hydrogens: bool = True,
            map_hydrogens_on_hydrogens_only: bool = False,
            map_exact_ring_matches_only: bool = True,
            additional_mapping_filter_functions: Optional[Iterable[Callable[[
                Chem.Mol, Chem.Mol, dict[int, int]], dict[int, int]]]] = None,
            _mapping_algorithm: str = mapping_algorithm.linear_sum_assignment,
    ):
        """ Geometry Based Atom Mapper
        This mapper is a homebrew, that utilises rdkit in order
        to generate an atom-mapping based on the coordinates of two molecules.

        Parameters
        ----------
        atom_max_distance : float, optional
            geometric criteria for two atoms, how far their distance
            can be maximal (in Angstrom). Default 0.95
        map_hydrogens_on_hydrogens_only : bool, optional
            map hydrogens only on hydrogens. Default False
        map_exact_ring_matches_only : bool, optional
            if true, only rings with matching ringsize and same bond-orders
            will be mapped. Additionally no ring-breaking is permitted. default
            False
        additional_mapping_filter_functions : Iterable[Callable[[Chem.Mol,
        Chem.Mol, Dict[int, int]], Dict[int, int]]], optional
            with this optional parameter you can further filter the distance
            based mappings with your own custom filters, provided as iterables.
            as default we suggest to avoid ring size/breaking changes and only
            allow whole rings to be mapped
        _mapping_algorithm : str, optional
            mapping_algorithm.linear_sum_assignment - this allows to swap the
            optimization algorithm. Not recommended.

        """
        self.atom_max_distance = atom_max_distance
        self.atom_map_hydrogens = atom_map_hydrogens
        self._map_exact_ring_matches_only = map_exact_ring_matches_only
        self._map_hydrogens_on_hydrogens_only = map_hydrogens_on_hydrogens_only

        self._filter_funcs = []
        if map_hydrogens_on_hydrogens_only:
            self._filter_funcs.append(filter_atoms_h_only_h_mapped)
        if self._map_exact_ring_matches_only:
            self._filter_funcs.extend(
                [
                    filter_ringbreak_changes,
                    filter_ringsize_changes,
                    filter_whole_rings_only,
                ]
            )
        if additional_mapping_filter_functions is not None:
            self._filter_funcs.extend(additional_mapping_filter_functions)

        if _mapping_algorithm == mapping_algorithm.linear_sum_assignment:
            self._map_hydrogens_on_hydrogens_only = True
            self.mapping_algorithm = self._linearSumAlgorithm_map
        elif _mapping_algorithm == mapping_algorithm.minimal_spanning_tree:
            self.mapping_algorithm = self._minimalSpanningTree_map
        else:
            raise ValueError(
                f"Mapping algorithm not implemented or unknown (options: MST "
                f"or LSA). got key: {_mapping_algorithm}"
            )

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self) -> dict:
        built_in_filters = {
            filter_atoms_h_only_h_mapped,
            filter_whole_rings_only,
            filter_ringsize_changes,
            filter_ringbreak_changes,
        }
        additional_filters = [
            dill.dumps(f) for f in self._filter_funcs
            if f not in built_in_filters
        ]

        # rather than serialise _filter_funcs, we serialise the arguments
        # that lead to the correct _filter_funcs being added
        #
        # then reverse engineer the _mapping_algorithm argument
        # this avoids serialising the function directly
        map_arg = {
            self._linearSumAlgorithm_map: 'LSA',
            self._minimalSpanningTree_map: 'MST',
        }[self.mapping_algorithm]

        return {
            'atom_max_distance': self.atom_max_distance,
            'atom_map_hydrogens': self.atom_map_hydrogens,
            'map_hydrogens_on_hydrogens_only': self._map_hydrogens_on_hydrogens_only,
            'map_exact_ring_matches_only': self._map_exact_ring_matches_only,
            '_mapping_algorithm': map_arg,
            'filters': additional_filters,
        }

    @classmethod
    def _from_dict(cls, d: dict):
        # replace _mapping_algorithm key to enum
        map_arg = d.pop('_mapping_algorithm', 'LSA')

        map_alg = {
            'LSA': mapping_algorithm.linear_sum_assignment,
            'MSR': mapping_algorithm.minimal_spanning_tree,
        }[map_arg]

        d['_mapping_algorithm'] = map_alg

        d['additional_mapping_filter_functions'] = [
            dill.loads(f) for f in d.pop('filters', [])
        ]

        return cls(**d)

    @property
    def map_hydrogens_on_hydrogens_only(self) -> bool:
        """this property is a shortcut for setting hydrogen shall be mapped
        only on hydrogen filter."""
        return self._map_hydrogens_on_hydrogens_only

    @map_hydrogens_on_hydrogens_only.setter
    def map_hydrogens_on_hydrogens_only(self, s: bool):
        self._map_hydrogens_on_hydrogens_only = s
        if s and filter_atoms_h_only_h_mapped not in self._filter_funcs:
            self._filter_funcs.insert(0, filter_atoms_h_only_h_mapped)
        elif filter_atoms_h_only_h_mapped in self._filter_funcs:
            self._filter_funcs.remove(filter_atoms_h_only_h_mapped)

    @property
    def map_exact_ring_matches_only(self) -> bool:
        """this property is a shortcut for setting exact ring matches only."""
        return self._map_exact_ring_matches_only

    @map_exact_ring_matches_only.setter
    def map_exact_ring_matches_only(self, s: bool):
        self._map_exact_ring_matches_only = s

        for f in [filter_ringbreak_changes, filter_ringsize_changes,
                  filter_whole_rings_only, ]:
            if s and f not in self._filter_funcs:
                self._filter_funcs.append(f)
            elif f in self._filter_funcs:
                self._filter_funcs.remove(f)

    """
       Private - Set Operations
    """

    @classmethod
    def _filter_mapping_for_max_overlapping_connected_atom_set(
            cls,
            moleculeA: Chem.Mol,
            moleculeB: Chem.Mol,
            atom_mapping: dict[int, int],
    ) -> dict[int, int]:
        """ Find connected core region from raw mapping
            This algorithm finds the maximal overlapping connected  set of
            two molecules and a given mapping. In order to accomplish this
            the following steps are executed:
            1. find all connnected sets for molecule A and molecule B
            2. get the maximal overlap between a set of molecule A and a
            set of molecule B
            3. reduce the atom_mapping to the found maximal overlapping sets.

        Parameters
        ----------
        moleculeA : Chem.Mol
            molecule A which mapped atoms are represented by atom_mapping keys
        moleculeB : Chem.Mol
            molecule B which mapped atoms are represented by atom_mapping values
        atom_mapping : Dict[int, int]
            input atom_mapping, to be filtered.

        Returns
        -------
        Dict[int, int]
            a filtered mapping, containing the maximal overlapping filter.

        """
        # get connected sets from mappings
        sets_a = cls._get_connected_atom_subsets(
            moleculeA, list(atom_mapping.keys())
        )
        sets_b = cls._get_connected_atom_subsets(
            moleculeB, list(atom_mapping.values())
        )

        log_var = "\t".join(map(str, [sets_a, sets_b, atom_mapping]))
        logger.debug(f"Found connected sets: {log_var}")

        # get maximally overlapping largest sets
        ((max_set_a_id, max_set_b_id), max_set,
         ) = cls._get_maximal_mapping_set_overlap(sets_a, sets_b, atom_mapping)

        # filter mapping for atoms in maximally overlapping sets
        atom_mapping = cls._filter_mapping_for_set_overlap(
            sets_a[max_set_a_id], sets_b[max_set_b_id], atom_mapping
        )

        set_overlaps = "\t".join(map(str, [max_set_a_id, max_set_b_id,
                                           max_set, atom_mapping]))

        logger.debug(f"Found connected set overlaps: {set_overlaps}")

        return atom_mapping

    @staticmethod
    def _get_connected_atom_subsets(
            mol: Chem.Mol, to_be_searched: list[int]
    ) -> list[set[int]]:
        """ find connected sets in mappings
        Get the connected sets of all to_be_searched atom indices in mol.
        Connected means the atoms in a resulting connected set are connected
        by covalent bonds.

        Parameters
        ----------
        mol : Chem.Mol
            molecule graph containing the bond information
        to_be_searched : List[int]
            atom indices, that shall be grouped by connections.

        Returns
        -------
        List[Set[int]]
            return list of connected sets from the to_be_searched atom index
            list found in mol.
        """
        # Get bonded atoms sets
        bonded_atom_sets = []
        for aid in to_be_searched:
            a = mol.GetAtomWithIdx(int(aid))
            bond_atoms = [b.GetOtherAtomIdx(aid) for b in a.GetBonds()]
            connected_atoms = set([aid] + bond_atoms)

            # filter to not connect non mapped atoms
            filtered_connected_atoms = set(
                filter(lambda x: x in to_be_searched, connected_atoms)
            )

            # NO mapped atom connected?
            if len(filtered_connected_atoms) == 1:
                continue
            else:
                bonded_atom_sets.append(filtered_connected_atoms)
        bonded_atom_sets = list(sorted(bonded_atom_sets))

        logger.debug(f"Bonded atom Sets: {bonded_atom_sets}")

        # add new set or merge with existing set
        atom_mappings = {i: a for i, a in enumerate(to_be_searched)}
        r_atom_mappings = {a: i for i, a in enumerate(to_be_searched)}
        max_n = len(to_be_searched)

        connection_matrix = np.zeros(shape=(max_n, max_n))
        for set_i in bonded_atom_sets:
            for ai in set_i:
                for aj in set_i:
                    connection_matrix[r_atom_mappings[ai]][
                        r_atom_mappings[aj]
                    ] = 1

        logger.debug(f"Adjacency Matrix: {connection_matrix}")

        # get connected components
        matrix = csr_matrix(connection_matrix)
        n_components, labels = connected_components(
            csgraph=matrix, directed=False
        )

        logger.debug(f"Connected Components: {n_components}\t{labels} ")

        # Translate back
        merged_connnected_atom_sets = []
        labels_sizes = [
            label
            for label, _ in sorted(
                np.vstack(np.unique(labels, return_counts=True)).T,
                key=lambda x: x[1],
                reverse=True,
            )
        ]
        for selected_label in labels_sizes:
            connected_atoms_ids = list(
                sorted(
                    [
                        atom_mappings[i]
                        for i, label in enumerate(labels)
                        if label == selected_label
                    ]
                )
            )
            merged_connnected_atom_sets.append(connected_atoms_ids)

        logger.debug(f"Merged Sets: {merged_connnected_atom_sets}")

        return merged_connnected_atom_sets

    @staticmethod
    def _get_maximal_mapping_set_overlap(
            sets_a: Iterable[set], sets_b: Iterable[set],
            mapping: dict[int, int]
    ) -> tuple[set, set]:
        """get the largest set overlaps in the mapping of set_a and set_b.

        Parameters
        ----------
        sets_a : Iterable[Set]
            connected sets resulting from the mapping in state A
        sets_b : Iterable[Set]
            connected sets resulting from the mapping in state B
        mapping : Dict[int, int]
            proposed atom mapping from state A to state B

        Returns
        -------
        Tuple[Set, Set]
            returns the found maximal and maximal overlapping sets of a and
            the set of b
        """
        # Calculate overlaps
        logger.debug(f"Input: {sets_a} {sets_b} {mapping}")

        max_set_combi = {}
        for ida, set_a in enumerate(sets_a):
            mapped_set_a = {mapping[a] for a in set_a}
            for idb, set_b in enumerate(sets_b):
                set_intersection_size = len(mapped_set_a.intersection(set_b))
                if set_intersection_size:
                    max_set_combi[(ida, idb)] = {
                        "overlap_count": set_intersection_size,
                        "set_a_size": len(set_a),
                        "set_b_size": len(set_b),
                    }

        # sort overlap by size
        overlap_sorted_sets = OrderedDict(
            sorted(
                max_set_combi.items(),
                key=lambda x: x[1]["overlap_count"],
                reverse=False,
            )
        )

        logger.debug(f"Set overlap properties {overlap_sorted_sets}")

        (
            max_overlap_set_a_id,
            max_overlap_set_b_id,
        ) = overlap_sorted_sets.popitem()
        return max_overlap_set_a_id, max_overlap_set_b_id

    @staticmethod
    def _filter_mapping_for_set_overlap(
            set_a: set[int], set_b: set[int], mapping: dict[int, int]
    ) -> dict[int, int]:
        """This filter reduces the mapping dict to only in the sets contained
        atom IDs

        Parameters
        ----------
        set_a : Set[int]
            this set contains the allowed keys
        set_b : Set[int]
            this set contains the allowed values
        mapping : Dict[int, int]
            this mapping is to be filtered by the set values.

        Returns
        -------
        Dict[int, int]
            filtered mapping

        """
        filtered_mapping = {}
        for atom_a, atom_b in mapping.items():
            if atom_a in set_a and atom_b in set_b:
                filtered_mapping[atom_a] = atom_b
        return filtered_mapping

    """
        Utils
    """

    @staticmethod
    def _get_full_distance_matrix(
            atomA_pos: NDArray,
            atomB_pos: NDArray,
            metric: Callable[
                [Union[float, Iterable], Union[float, Iterable]],
                Union[float, Iterable],
            ] = vector_eucledean_dist,
    ) -> np.array:
        """calculates a full distance matrix between the two given input
        position matrixes.


        Parameters
        ----------
        atomA_pos : NDArray
            position matrix A
        atomB_pos : NDArray
            position matrix B
        metric : Callable[[Union[float, Iterable], Union[float, Iterable]],
        Union[float, Iterable]], optional
            the applied metric to calculate the distance matrix. default
            metric: eucledean distance.

        Returns
        -------
        np.array
            returns a distance matrix.

        """
        distance_matrix = []
        for atomPosA in atomA_pos:
            atomPos_distances = metric(atomPosA, atomB_pos)
            distance_matrix.append(atomPos_distances)
        distance_matrix = np.array(distance_matrix)
        return distance_matrix

    @staticmethod
    def _mask_atoms(
            mol, mol_pos, map_hydrogens: bool, masked_atoms: list[int],
    ) -> tuple[dict, list]:
        """Mask atoms such they are not considered during the mapping.

        Parameters
        ----------
        mol: Chem.Mol
            molecule to be masked
        mol_pos: List[List[int]]
            positions of mol.
        map_hydrogens: bool
            should hydrogens be mapped?
        masked_atoms: List[int]
            already masked atoms.

        Returns
        -------
        Tuple[Dict, List]
            the return values are a dictionary, mapping atom ids onto the new
            reduced dimensions and reduced positions.
        """
        pos = []
        masked_atomMapping = {}
        for atom in mol.GetAtoms():
            atom_id = atom.GetIdx()
            if (
                    map_hydrogens or atom.GetSymbol() != "H"
            ) and atom_id not in masked_atoms:
                masked_atomMapping[len(masked_atomMapping)] = atom_id
                pos.append(mol_pos[atom_id])
        pos = np.array(pos)

        return masked_atomMapping, pos

    def _minimalSpanningTree_map(
            self, distance_matrix: NDArray, max_dist: float
    ) -> dict[int, int]:
        """MST Mapping
        This function is a numpy graph based implementation to build up an
        Atom Mapping purely on 3D criteria.

        Parameters
        ----------
        distance_matrix: NDArray
            distances of atoms to each other.
        max_dist: float
            maximal distance of a atoms in a mapping.

        Returns
        -------
        Dict[int, int]
            mapping of atoms.
        """

        # distance matrix:  - full graph
        logger.debug(f"Got Distance Matrix: \n {distance_matrix}")
        edges = []
        for i, distance_row in enumerate(distance_matrix):
            for j, dist in enumerate(distance_row):
                edges.append([float(dist), int(i), int(j)])

        edges = np.array(edges)
        logger.debug(f"Edges: \n{edges}")

        # priority queue based on edge weight
        sorted_edges = edges[edges[:, 0].argsort()]
        logger.debug(f"Sorted Edges: \n{sorted_edges}")

        # MST like algorithm
        mapping = {}
        for w, x, y in sorted_edges:
            if w >= max_dist:  # filter for max dist
                break
            else:
                if x not in mapping.keys() and y not in mapping.values():
                    mapping[int(x)] = int(y)

        return mapping

    @staticmethod
    def _linearSumAlgorithm_map(
            distance_matrix: NDArray, max_dist: float
    ) -> dict[int, int]:
        """ LSA mapping
        This function is a LSA based implementation to build up an Atom
        Mapping purely on 3D criteria.

        Parameters
        ----------
        distance_matrix: NDArray
            distances of atoms to each other.
        max_dist: float
            maximal distance of a atoms in a mapping.

        Returns
        -------
        Dict[int, int]
            mapping of atoms.
        """
        row_ind, col_ind = linear_sum_assignment(distance_matrix)
        raw_mapping = list(zip(map(int, row_ind), map(int, col_ind)))
        # filter for mask value
        mapping = dict(
            filter(lambda x: distance_matrix[x] < max_dist, raw_mapping)
        )

        return mapping

    def _additional_filter_rules(
            self, molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
    ) -> dict[int, int]:
        """apply additional filter rules to the given mapping.

        Parameters
        ----------
        molA : Chem.Mol
            mol, with atoms contained in mapping keys.
        molB : Chem.Mol
            mol, with atoms contained in mapping values.
        mapping : Dict[int, int]
            mapping to be filtered

        Returns
        -------
        Dict[int, int]
            filtered mapping
        """
        logger.debug(f"Before filters mapping is {mapping}")
        for filter_rule in self._filter_funcs:
            mapping = filter_rule(molA, molB, mapping)
            logger.debug(f"After {filter_rule} mapping is {mapping}")
        return mapping

    def _calculate_mapping(
            self,
            molA: Chem.Mol,
            molB: Chem.Mol,
            max_d: float = 0.95,
            masked_atoms_molA: Optional[list[int]] = None,
            masked_atoms_molB: Optional[list[int]] = None,
            pre_mapped_atoms: Optional[dict[int, int]] = None,
            map_hydrogens: bool = True,
    ) -> dict[int, int]:
        """ main mapping function - private
            find a mapping between two molecules based on 3D coordinates.
            This is a helper function for the suggest mapping functions.

        Parameters
        ----------
        molA : Chem.Mol
            rdkit Molecule A
        molB : Chem.Mol
            rdkit Molecule B
        max_d : float, optional
            maximal allowed distance for atoms to be mapped, by default 0.95
        masked_atoms_molA : list, optional
            remove categorical atoms by index of MolA from the mapping,
            by default []
        masked_atoms_molB : list, optional
            remove categorical atoms by index of MolB from the mapping,
            by default []
        pre_mapped_atoms : dict, optional
            pre_mapped atoms, that need to be part of the mapping, by default {}
        map_hydrogens : bool, optional
            if True map hydrogens as well, will hapen after heavy atom
            mapping, by default True

        Returns
        -------
        AtomMapping
            _description_
        """
        if masked_atoms_molA is None:
            masked_atoms_molA = []
        if masked_atoms_molB is None:
            masked_atoms_molB = []
        if pre_mapped_atoms is None:
            pre_mapped_atoms = dict()

        molA_pos = molA.GetConformer().GetPositions()
        molB_pos = molB.GetConformer().GetPositions()
        masked_atoms_molA = copy.deepcopy(masked_atoms_molA)
        masked_atoms_molB = copy.deepcopy(masked_atoms_molB)
        pre_mapped_atoms = copy.deepcopy(pre_mapped_atoms)

        if len(pre_mapped_atoms) > 0:
            masked_atoms_molA.extend(pre_mapped_atoms.keys())
            masked_atoms_molB.extend(pre_mapped_atoms.values())

        logger.debug(
            f"Positions lengths: {len(molA_pos)} {len(molB_pos)} ")

        pos = "\n\t".join(map(str, [molA_pos, molB_pos]))
        logger.debug(f"Positions: \n{pos}")
        logger.debug(f"masked_atoms_molA: {masked_atoms_molA}")
        logger.debug(f"masked_atoms_molB: {masked_atoms_molB}")
        logger.debug(f"pre_mapped_atoms: {pre_mapped_atoms}")
        logger.debug(f"map_hydrogens: {map_hydrogens}")

        logger.info("Masking Atoms")

        molA_masked_atomMapping, molA_pos = self._mask_atoms(
            mol=molA,
            mol_pos=molA_pos,
            masked_atoms=masked_atoms_molA,
            map_hydrogens=map_hydrogens,
        )
        molB_masked_atomMapping, molB_pos = self._mask_atoms(
            mol=molB,
            mol_pos=molB_pos,
            masked_atoms=masked_atoms_molB,
            map_hydrogens=map_hydrogens,
        )

        masked_atom_ids_str = "\n\t".join(map(str, [molA_masked_atomMapping,
                                                    molB_masked_atomMapping]))
        unmasked_atom_ids_str = "\n\t".join(map(str, [molA_pos, molB_pos]))
        logger.debug(f"Remove Masked Atoms:\n{masked_atom_ids_str}\n"
                     f"{unmasked_atom_ids_str}")
        logger.debug(f"filtered Masked Positions lengths:{len(molA_pos)}\t"
                     f"{len(molB_pos)}")

        if len(molA_pos) == 0 or len(molB_pos) == 0:
            if len(pre_mapped_atoms) == 0:
                logger.warning("no mappable Atoms were found!")
            return pre_mapped_atoms

        # Calculate mapping
        logger.info("Build Distance Matrix")
        # distance matrix:  - full graph
        distance_matrix = self._get_full_distance_matrix(molA_pos, molB_pos)
        logger.debug(f"Distance Matrix: \n{distance_matrix}")

        # Mask distance matrix with max_d
        # np.inf is considererd as not possible in lsa implementation -
        # therefore use a high value
        self.mask_dist_val = max_d * 10 ** 6
        masked_dmatrix = np.array(
            np.ma.where(
                distance_matrix < max_d, distance_matrix, self.mask_dist_val
            )
        )
        logger.debug(f"masked Distance Matrix: {max_d}\n\t{masked_dmatrix}")

        # solve atom mappings
        logger.info("Calculate Mapping")
        mapping = self.mapping_algorithm(
            distance_matrix=masked_dmatrix, max_dist=self.mask_dist_val
        )
        logger.debug(f"Raw Mapping: {mapping}")

        if len(mapping) == 0:  # TODO: check if this is correct
            if len(pre_mapped_atoms) == 0:
                logger.warning("no mapping could be found!")
            return pre_mapped_atoms

        # reverse any prior masking:
        mapping = {
            molA_masked_atomMapping[k]: molB_masked_atomMapping[v]
            for k, v in mapping.items()
        }

        # filter mapping for rules:
        if self._filter_funcs is not None:
            mapping = self._additional_filter_rules(molA, molB, mapping)

        if len(pre_mapped_atoms) > 0:
            mapping.update(pre_mapped_atoms)
        logger.debug(f"reverse Masking Mapping: {mapping}")

        if len(mapping) == 0:
            if len(pre_mapped_atoms) == 0:
                logger.warning("no mapping could be found, after applying "
                               "filters!")
            return pre_mapped_atoms

        # Reduce mapping to maximally overlapping two connected sets
        logger.info("Find Maximal overlapping connected sets of mapped atoms")
        mapping = self._filter_mapping_for_max_overlapping_connected_atom_set(
            moleculeA=molA, moleculeB=molB, atom_mapping=mapping
        )
        logger.debug(f"Set overlap Mapping: {mapping}")

        return mapping

    """
        Mapping functions
    """

    def suggest_mapping_from_rdmols(
            self,
            molA: Chem.Mol,
            molB: Chem.Mol,
            masked_atoms_molA: Optional[list] = None,
            masked_atoms_molB: Optional[list] = None,
            pre_mapped_atoms: Optional[dict] = None,
    ) -> dict[int, int]:
        """ Mapping Function with RDkit
        The function effectivly maps the two molecules on to each other and
        applies the given settings by the obj.

        Parameters
        ----------
        molA, molB : rdkit.Chem.Mol
            two rdkit molecules that should be mapped onto each other
        masked_atoms_molA : list, optional
            remove categorical atoms by index of MolA from the mapping,
            by default []
        masked_atoms_molB : list, optional
            remove categorical atoms by index of MolB from the mapping,
            by default []
        pre_mapped_atoms : dict, optional
            pre_mapped atoms, that need to be part of the mapping, by default {}
        """
        if masked_atoms_molA is None:
            masked_atoms_molA = []
        if masked_atoms_molB is None:
            masked_atoms_molB = []
        if pre_mapped_atoms is None:
            pre_mapped_atoms = {}

        molA = Chem.Mol(molA)
        molB = Chem.Mol(molB)
        masked_atoms_molA = copy.deepcopy(masked_atoms_molA)
        masked_atoms_molB = copy.deepcopy(masked_atoms_molB)
        pre_mapped_atoms = copy.deepcopy(pre_mapped_atoms)

        logger.info("#################################")
        logger.info("Map Heavy Atoms ")
        logger.info("#################################")

        mapping = self._calculate_mapping(
            molA,
            molB,
            max_d=self.atom_max_distance,
            masked_atoms_molA=masked_atoms_molA,
            masked_atoms_molB=masked_atoms_molB,
            pre_mapped_atoms=pre_mapped_atoms,
            map_hydrogens=False,
        )

        if self.atom_map_hydrogens:
            logger.info("#################################")
            logger.info("Map Hydrogen Atoms: ")
            logger.info("#################################")

            pre_mapped_atoms.update(mapping)

            mapping = self._calculate_mapping(
                molA,
                molB,
                max_d=self.atom_max_distance,
                masked_atoms_molA=masked_atoms_molA,
                masked_atoms_molB=masked_atoms_molB,
                pre_mapped_atoms=pre_mapped_atoms,
                map_hydrogens=self.atom_map_hydrogens,
            )

        return mapping

    def _split_protein_components_molecules(self, protein: ProteinComponent) -> list[ProteinComponent]:
        """
        Aims at splitting a protein component into different protein components based on the
        connectivity of the molecules that compose it. Useful for mapping multimer components
        or proteins with structural waters or similarly.

        This returns a ``ProteinComponent`` object with the name attribute indicating the starting
        index in the original component.
        """
        from rdkit.Chem.rdmolops import GetMolFrags
        rdmol = protein.to_rdkit()
        fragments_indices = GetMolFrags(rdmol)
        components = []
        for index_tuple in fragments_indices:
            edit_rdmol_frag = Chem.EditableMol(rdmol)
            remove_indices = []
            for atom in rdmol.GetAtoms():
                atom_index = atom.GetIdx()
                if not (atom_index in index_tuple):
                    remove_indices.append(atom_index)
            # Need to remove separately https://github.com/rdkit/rdkit/issues/1366
            for atom_idx in sorted(remove_indices, reverse=True):
                edit_rdmol_frag.RemoveAtom(atom_idx)
            #  Create component with the remaining molecule
            frag_rdmol = edit_rdmol_frag.GetMol()
            Chem.SanitizeMol(frag_rdmol)
            # FIXME: We are storing the starting index in name, this is super hacky :/
            protein_comp = ProteinComponent(frag_rdmol, name=f"frag_{index_tuple[0]}")
            components.append(protein_comp)

        return components

    def suggest_mappings(
            self, A: SmallMoleculeComponent, B: SmallMoleculeComponent
    ) -> Iterator[AtomMapping]:
        """ Mapping generator - Gufe
        return a generator for atom mappings.

        Parameters
        ----------
        A : SmallMoleculeComponent
            molecule A to be mapped.
        B : SmallMoleculeComponent
            molecule B to be mapped.

        Returns
        -------
        Iterator[AtomMapping]
            returns an interator of possible atom mappings.
        """
        # TODO: We probably want to modularize the following in methods
        if isinstance(A, ProteinComponent) or isinstance(B, ProteinComponent):
            # 1. identify Component Chains
            componentA_chains = self._split_protein_components_molecules(A)
            componentB_chains = self._split_protein_components_molecules(B)

            # 2. calculate all possible mappings
            largest_mappings = []
            for A_chain in componentA_chains:
                largest_overlap_map = {}  # Initialize to empty map
                largest_overlap_component = componentB_chains[0]  # Initialization
                for B_chain in componentB_chains:
                    # This reinitializes indices, that's why we need the index information from
                    #  split components.
                    current_map = self.suggest_mapping_from_rdmols(
                        molA=A_chain.to_rdkit(), molB=B_chain.to_rdkit()
                    )
                    if len(largest_overlap_map) < len(current_map):
                        largest_overlap_component = B_chain
                        largest_overlap_map = current_map
                # TODO: Do we need a better suited object here instead of LigandAtomMapping?
                mapping_obj = LigandAtomMapping(A_chain, largest_overlap_component, largest_overlap_map)
                # At the end of the loop mapping_obj should have the largest map overlap
                largest_mappings.append(mapping_obj)

            # Merge all the largest mappings for each component into a single mapping
            merged_map = {}
            for mapping_obj in largest_mappings:
                start_a = int(mapping_obj.componentA.name.split("_")[-1])
                start_b = int(mapping_obj.componentB.name.split("_")[-1])
                shifted_map = {a_idx + start_a: b_idx + start_b for a_idx, b_idx in
                               mapping_obj.componentA_to_componentB.items()}
                merged_map.update(shifted_map)
            yield LigandAtomMapping(A, B, merged_map)

        else:  # SmallMoleculeComponent case
            yield LigandAtomMapping(
                A,
                B,
                self.suggest_mapping_from_rdmols(
                    molA=A.to_rdkit(), molB=B.to_rdkit()
                ),
            )
