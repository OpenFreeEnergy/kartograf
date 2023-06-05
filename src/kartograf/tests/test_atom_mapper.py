# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest
from kartograf import KartografAtomMapper
from kartograf.atom_mapper import filter_atoms_h_only_h_mapped

from .conf import (
    naphtalene_benzene_molecules,
    naphtalene_benzene_mapping,
    stereco_chem_molecules,
    stereo_chem_mapping,
)
from copy import deepcopy

def check_mapping_vs_expected(mapping, expected_mapping):
    assert len(expected_mapping) == len(mapping.componentA_to_componentB)

    diff = []
    for exp_k, exp_v in expected_mapping.items():
        if exp_k not in mapping.componentA_to_componentB:
            diff.append((exp_k, exp_v, None, "missing mapping for atom"))
        elif exp_v != mapping.componentA_to_componentB[exp_k]:
            diff.append(
                (
                    exp_k,
                    exp_v,
                    mapping.componentA_to_componentB[exp_k],
                    "-",
                    "differently mapped atoms",
                )
            )

    if len(diff) > 0:
        print("Differences: expected key, expected value, actual value")
        print("\n".join(map(lambda x: "\t".join(map(str, x)), diff)))
        raise ValueError("mapping did not match expected mapping!")


def test_mapping_naphtalene_benzene(
    naphtalene_benzene_molecules, naphtalene_benzene_mapping
):
    """
    Test mapping of naphtalene to benzene.
    """
    expected_mapping = naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_naphtalene_benzene_noHs(naphtalene_benzene_molecules):
    """
    Test mapping of naphtalene to benzene without H-atoms.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=True,
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_naphtalene_benzene_noHs_add_filter(naphtalene_benzene_molecules):
    """
    Test mapping of naphtalene to benzene without H-atoms added as additional filter.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
        additional_mapping_filter_functions=[filter_atoms_h_only_h_mapped],
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_stereo_mapping(stereco_chem_molecules, stereo_chem_mapping):
    """
    Test weird Stereochemistry mapping.
    """
    expected_mapping = stereo_chem_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95, atom_map_hydrogens=True
    )  # mapping_algorithm.minimal_spanning_tree
    geom_mapping = next(
        geom_mapper.suggest_mappings(
            stereco_chem_molecules[0],
            stereco_chem_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_to_From_dict():
    mapper = KartografAtomMapper()
    d1 = mapper._to_dict()
    mapper2 = KartografAtomMapper._from_dict(d1)
    d2 = mapper2._to_dict()

    for key, val1 in d1.items():
        val2 = d2[key]

        if(val1 != val2):
            raise ValueError("they need to be identical.")

    mapper2.atom_max_distance = 10
    d3 = mapper2._to_dict()
    print(d3)
    for key, val1 in d1.items():
        val2 = d3[key]

        if(val1 != val2 and key != "atom_max_distance"):
            raise ValueError("they need to be identical.")
        if(key == "atom_max_distance" and val1 == val2 ):
            raise ValueError("they must not be identical.")


def test_filter_property():
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=False)
    first_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = True
    second_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = False
    third_filters = deepcopy(mapper._filter_funcs)

    assert len(first_filters) == len(third_filters)
    assert len(first_filters) == len(second_filters)-1

    assert filter_atoms_h_only_h_mapped in second_filters
    assert filter_atoms_h_only_h_mapped not in first_filters
    assert filter_atoms_h_only_h_mapped not in first_filters

